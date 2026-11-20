"""Stamp provenance onto the existing result files.

`../src` is shared with other work in this project and was touched by another
process while these analyses were running, so each result file gets

  * the SHA-256 of the source files, and
  * a numeric re-check that the baseline run still reproduces the values the
    analyses were computed with.

No result number is modified here: only `meta` fields are added.
"""

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'src'))
sys.path.insert(0, HERE)

import dv_common as C
import model as M
import r0 as R0

TARGETS = ['dynamic_prcc.json', 'nv_scaling.json', 'clamp_sensitivity.json',
           'equalwidth_prcc.json', 'static_prcc_reference.json',
           'constant_param_prcc.json']

EXPECTED = {
    'peak_ih_n1096': 1650.756836,
    'peak_day_n1096': 233.848208,
    'cum_infected_n1096': 0.968944,
    'peak_ih_default_grid_n1000': 1649.777450,
    'R0_at_28C': 2.928116,
}


def main():
    C.reset_rates()
    m = C.run_case(M.baseline_params())
    t, y = M.solve(M.baseline_params())          # default 1000-point grid
    check = {
        'peak_ih_n1096': m['peak_ih'],
        'peak_day_n1096': m['peak_day'],
        'cum_infected_n1096': m['cum_infected'],
        'peak_ih_default_grid_n1000': float(np.max(y[2])),
        'R0_at_28C': R0.compute_r0(M.baseline_params(), 28.0),
    }
    ok = (abs(check['peak_ih_n1096'] - EXPECTED['peak_ih_n1096']) < 1e-5
          and abs(check['peak_ih_default_grid_n1000'] - EXPECTED['peak_ih_default_grid_n1000']) < 1e-5
          and abs(check['cum_infected_n1096'] - EXPECTED['cum_infected_n1096']) < 1e-6
          and abs(check['R0_at_28C'] - EXPECTED['R0_at_28C']) < 1e-5)
    print(f"[verify] baseline reproduction against the values used in the runs: {ok}")
    for k, v in check.items():
        print(f"   {k:<28} {v:.9f}  (expected {EXPECTED[k]})")

    stamp = {
        'src_sha256': C.src_fingerprint(),
        'src_numeric_recheck': check,
        'src_numeric_recheck_matches_run_values': bool(ok),
    }

    for name in TARGETS:
        path = os.path.join(HERE, name)
        with open(path, encoding='utf-8') as fh:
            data = json.load(fh)
        data.setdefault('meta', {}).update(stamp)
        with open(path, 'w', encoding='utf-8') as fh:
            json.dump(data, fh, indent=2, ensure_ascii=False)
        print(f"[stamped] {name}")


if __name__ == '__main__':
    main()
