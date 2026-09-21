#!/usr/bin/env python3
"""Run marker-image and genome analyses."""
import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_name, '1')
if __name__ == '__main__':
    from legacy_marker_rescue.cli import main
    raise SystemExit(main())
