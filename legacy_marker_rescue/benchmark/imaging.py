"""Image-only intensity/scale adapter; reference masks have a separate decoder."""
from __future__ import annotations
import io
import math
import numpy as np
from PIL import Image
from scipy.ndimage import label, find_objects

def decode_image(data: bytes, preprocessing: dict) -> tuple[np.ndarray, dict]:
    with Image.open(io.BytesIO(data)) as im:
        if getattr(im, 'n_frames', 1) != 1:
            raise ValueError('MULTIFRAME_IMAGE_UNSUPPORTED: frames are not silently collapsed')
        if im.width * im.height > preprocessing['max_image_pixels']:
            raise ValueError('IMAGE_PIXEL_SAFETY_LIMIT')
        mode = im.mode
        if mode == 'RGBA':
            if np.asarray(im)[..., 3].min() != 255:
                raise ValueError('TRANSPARENT_IMAGE_UNSUPPORTED')
            arr = np.asarray(im.convert('L'), dtype=float)
        elif mode in ('RGB', 'P'):
            arr = np.asarray(im.convert('L'), dtype=float)
        elif mode in ('L', 'I', 'F', 'I;16', 'I;16L', 'I;16B', '1'):
            arr = np.array(im, dtype=float)
        else:
            raise ValueError('UNSUPPORTED_IMAGE_MODE: ' + mode)
    if arr.ndim != 2 or not np.isfinite(arr).all():
        raise ValueError('INVALID_IMAGE_ARRAY')
    lo, hi = np.percentile(arr, [preprocessing['lower_percentile'], preprocessing['upper_percentile']])
    fallback = False
    if hi <= lo:
        lo = float(arr.min())
        hi = float(arr.max())
        fallback = True
    flat = hi <= lo
    normalized = np.zeros_like(arr) if flat else np.clip((arr - lo) * (255.0 / (hi - lo)), 0, 255)
    inverted = bool(not flat and np.median(normalized) > 127.5)
    if inverted:
        normalized = 255.0 - normalized
    return (normalized, {'original_mode': mode, 'width': arr.shape[1], 'height': arr.shape[0], 'original_min': float(arr.min()), 'original_max': float(arr.max()), 'normalization_low': float(lo), 'normalization_high': float(hi), 'normalization_fallback_minmax': fallback, 'constant_image': bool(flat), 'inverted_from_image_median': inverted, 'geometric_transform': 'none', 'note': 'Analysis working copy only. Original source bytes unchanged; no reference labels used.'})

def rescale_height(gray: np.ndarray, height: int) -> tuple[np.ndarray, float, float]:
    h, w = gray.shape
    nw = max(1, int(round(w * height / h)))
    im = Image.fromarray(gray.astype(np.float32))
    transformed = np.asarray(im.resize((nw, height), resample=Image.Resampling.BILINEAR), dtype=float)
    return (np.clip(transformed, 0, 255), nw / w, height / h)

def decode_mask(data: bytes, shape: tuple[int, int]) -> tuple[np.ndarray, dict]:
    with Image.open(io.BytesIO(data)) as im:
        if getattr(im, 'n_frames', 1) != 1:
            raise ValueError('MULTIFRAME_MASK_UNSUPPORTED')
        if im.size != (shape[1], shape[0]):
            raise ValueError('MASK_IMAGE_DIMENSION_MISMATCH')
        arr = np.asarray(im)
        if arr.ndim != 2:
            raise ValueError('REFERENCE_MASK_MUST_BE_SINGLE_CHANNEL')
        values = np.unique(arr)
    if not np.isfinite(values).all():
        raise ValueError('NONFINITE_MASK')
    allowed01 = np.isin(values, [0, 1]).all()
    allowed255 = np.isin(values, [0, 255]).all()
    if not (allowed01 or allowed255):
        raise ValueError('UNRECOGNIZED_MASK_ENCODING: ' + repr(values[:20].tolist()))
    return (arr != 0, {'original_values': values.tolist(), 'foreground_value': 1 if allowed01 else 255, 'background_value': 0})

def component_reference(mask: np.ndarray) -> tuple[np.ndarray, list[dict]]:
    """All 8-connected foreground objects, including single pixels and edge objects.

    Touching band annotations remain one component; no guessed splitting, area
    pruning, or biological interpretation is added.
    """
    if mask.ndim != 2:
        raise ValueError('Mask must be 2D')
    labels, n = label(mask, structure=np.ones((3, 3), dtype=np.uint8))
    objects = find_objects(labels)
    counts = np.bincount(labels.ravel(), minlength=n + 1)
    rows = []
    for i, sl in enumerate(objects, 1):
        if sl is None:
            raise ValueError('Noncontiguous component identifiers')
        yy, xx = np.nonzero(labels[sl] == i)
        ys = yy + sl[0].start + 0.5
        xs = xx + sl[1].start + 0.5
        rows.append({'component_id': i, 'x0': sl[1].start, 'x1': sl[1].stop, 'y0': sl[0].start, 'y1': sl[0].stop, 'x_center': float(np.mean(xs)), 'y_center': float(np.mean(ys)), 'area_px': int(counts[i]), 'touches_image_edge': bool(sl[0].start == 0 or sl[1].start == 0 or sl[0].stop == mask.shape[0] or (sl[1].stop == mask.shape[1]))})
    return (labels, rows)
