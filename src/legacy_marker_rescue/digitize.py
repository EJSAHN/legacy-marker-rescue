from __future__ import annotations
from pathlib import Path
from PIL import Image
import fitz


def render_pdf_page(pdf_path: str | Path, page_index: int, dpi: int = 300) -> Image.Image:
    doc = fitz.open(str(pdf_path))
    page = doc[page_index]
    scale = dpi / 72.0
    pix = page.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=False)
    return Image.frombytes("RGB", [pix.width, pix.height], pix.samples)


def crop_image(image: Image.Image, crop_box: list[int]) -> Image.Image:
    return image.crop(tuple(crop_box))


def detect_vertical_segments(image: Image.Image, dark_threshold: int = 80, min_height: int = 8, min_width: int = 1) -> list[dict]:
    gray = image.convert("L")
    w, h = gray.size
    px = gray.load()
    visited = set()
    segments = []
    for y in range(h):
        for x in range(w):
            if (x, y) in visited or px[x, y] > dark_threshold:
                continue
            stack = [(x, y)]
            visited.add((x, y))
            xs, ys = [], []
            while stack:
                cx, cy = stack.pop()
                xs.append(cx); ys.append(cy)
                for nx in (cx-1, cx, cx+1):
                    for ny in (cy-1, cy, cy+1):
                        if 0 <= nx < w and 0 <= ny < h and (nx, ny) not in visited and px[nx, ny] <= dark_threshold:
                            visited.add((nx, ny))
                            stack.append((nx, ny))
            xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)
            if (ymax - ymin + 1) >= min_height and (xmax - xmin + 1) >= min_width:
                segments.append({
                    "x_min": xmin, "x_max": xmax, "y_min": ymin, "y_max": ymax,
                    "x_mid": (xmin + xmax) / 2, "y_mid": (ymin + ymax) / 2,
                    "height": ymax - ymin + 1, "width": xmax - xmin + 1,
                })
    return segments


def x_to_kb(x: float, x0: float, x25: float, min_kb: float = 0.0, max_kb: float = 2.5) -> float:
    return min_kb + (x - x0) * (max_kb - min_kb) / (x25 - x0)
