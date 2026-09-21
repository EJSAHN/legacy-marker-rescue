"""Extract human lane geometry without band marks or exclusion masks."""
from __future__ import annotations
import math
from typing import Any

def finite(value: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or (not math.isfinite(float(value))):
        raise ValueError('Coordinates must be finite numbers')
    return float(value)

def intersect(a: list[float], b: list[float]) -> bool:
    return min(a[2], b[2]) > max(a[0], b[0]) and min(a[3], b[3]) > max(a[1], b[1])

def strip_geometry(review: dict) -> dict:
    """Only allowed geometry/category fields cross the detector boundary.

    In particular, annotations, exclusions, reviewed flags, counts, and events
    are not included. A change to the reader's clicks cannot change a profile.
    """
    geometry = {'schema_version': '1.0', 'coordinate_system': 'native_pixel_edges', 'image_sha256': review['image_sha256'], 'width_px': review['width_px'], 'height_px': review['height_px'], 'panels': {}}
    for pid, panel in review['panels'].items():
        geometry['panels'][pid] = {'lanes': [{k: lane[k] for k in ('id', 'box', 'label', 'label_source', 'role', 'assessment')} for lane in panel['lanes']]}
    validate_geometry(geometry)
    return geometry

def validate_geometry(geometry: dict) -> None:
    allowed = {'schema_version', 'coordinate_system', 'image_sha256', 'width_px', 'height_px', 'panels'}
    if set(geometry) != allowed or geometry['schema_version'] != '1.0' or geometry['coordinate_system'] != 'native_pixel_edges':
        raise ValueError('Unknown or non-geometry detector input fields')
    w, h = (geometry['width_px'], geometry['height_px'])
    if isinstance(w, bool) or isinstance(h, bool) or (not isinstance(w, int)) or (not isinstance(h, int)) or (min(w, h) < 1):
        raise ValueError('Invalid native image dimensions')
    if not isinstance(geometry['image_sha256'], str) or not geometry['image_sha256']:
        raise ValueError('Image hash missing')
    if not isinstance(geometry['panels'], dict) or not geometry['panels']:
        raise ValueError('No panel geometry supplied')
    ids = set()
    for pid, panel in geometry['panels'].items():
        if not isinstance(pid, str) or not pid or set(panel) != {'lanes'} or (not isinstance(panel['lanes'], list)):
            raise ValueError('Invalid panel geometry')
        prior = []
        for lane in panel['lanes']:
            if set(lane) != {'id', 'box', 'label', 'label_source', 'role', 'assessment'}:
                raise ValueError('Lane input contains unknown fields; band annotations are not detector inputs')
            lid = lane['id']
            if not isinstance(lid, str) or not lid or lid in ids:
                raise ValueError('Missing/duplicate lane identifier')
            ids.add(lid)
            if not isinstance(lane['box'], list) or len(lane['box']) != 4:
                raise ValueError('Invalid lane box')
            box = [finite(x) for x in lane['box']]
            if not (0 <= box[0] < box[2] <= w and 0 <= box[1] < box[3] <= h):
                raise ValueError('Lane box outside image or reversed')
            if any((intersect(box, other) for other in prior)):
                raise ValueError('Overlapping lanes in one panel')
            prior.append(box)
            if lane['role'] not in {'sample', 'marker', 'gap', 'unresolved'}:
                raise ValueError('Unknown lane role')
            if lane['assessment'] not in {'readable', 'unreadable'}:
                raise ValueError('Unknown readability')
            if not isinstance(lane['label'], str) or lane['label_source'] not in {'printed_here', 'shared_header', 'unresolved'}:
                raise ValueError('Invalid label provenance')

def pixel_range(lo: float, hi: float) -> tuple[int, int]:
    """Integer array indices whose pixel CENTERS are in half-open [lo, hi)."""
    return (math.ceil(lo - 0.5), math.ceil(hi - 0.5))
