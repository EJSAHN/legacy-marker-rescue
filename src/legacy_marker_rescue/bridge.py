from __future__ import annotations


def match_legacy_to_modern(legacy_bands: list[dict], modern_bands: list[dict], tolerance_kb: float) -> list[dict]:
    rows = []
    for legacy in legacy_bands:
        lk = float(legacy["kb"])
        for modern in modern_bands:
            mk = float(modern["bp"]) / 1000.0
            delta = abs(lk - mk)
            if delta <= tolerance_kb:
                out = dict(legacy)
                out.update({
                    "modern_band": modern["band"],
                    "modern_kb": mk,
                    "delta_kb": delta,
                    "modern_prevalence": modern.get("prevalence", ""),
                    "modern_class": modern.get("band_class", ""),
                })
                rows.append(out)
    return rows
