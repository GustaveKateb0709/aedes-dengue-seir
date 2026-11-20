"""Recompute the published static PRCC (R0 at 28 degC) for side-by-side
comparison with the dynamic analyses.

Nothing in ../src is modified and no published number is altered: this simply
re-runs ``sensitivity.run()`` with its own defaults (n=2000, seed 20260816,
T=28 degC) and stores the output as JSON next to the new results, so the three
estimators (static / dynamic / equal-width) can be tabulated together.
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'src'))
sys.path.insert(0, HERE)

import sensitivity as S
import dv_common as C


def main():
    r = S.run()
    payload = {
        'source': 'sensitivity.run() from ../src (recomputed, unchanged defaults)',
        'estimator': 'PRCC vs algebraic R0 at fixed T = 28 degC',
        'param_names': list(r['param_names']),
        'param_labels': list(r['param_labels']),
        'prcc': [float(x) for x in r['prcc']],
        'p_value': [float(x) for x in r['p_value']],
        'r0_min': float(r['r0_min']),
        'r0_max': float(r['r0_max']),
        'r0_median': float(r['r0_median']),
        'n_samples': int(r['n_samples']),
        'T': float(r['T']),
        'rank': list(range(1, len(r['prcc']) + 1)),
        'meta': C.make_meta(script='make_static_reference.py', seed=20260816,
                            n_samples=int(r['n_samples']), wall_seconds=0.0,
                            per_run_seconds=0.0,
                            extra={'note': 'no ODE solves; algebraic R0 only'}),
    }
    print(f"{'rank':<5}{'param':<12}{'PRCC':>9}{'p':>12}")
    for i, name in enumerate(payload['param_names']):
        print(f"{i + 1:<5}{name:<12}{payload['prcc'][i]:>9.3f}"
              f"{C.fmt_p(payload['p_value'][i]):>12}")
    print(f"R0 range {payload['r0_min']:.2f}-{payload['r0_max']:.2f} "
          f"(median {payload['r0_median']:.2f})")
    C.write_json(os.path.join(HERE, 'static_prcc_reference.json'), payload)


if __name__ == '__main__':
    main()
