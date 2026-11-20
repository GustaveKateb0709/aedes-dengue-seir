"""Runtime configuration for the dengue modelling package.

The only tunable global is the pseudo-random seed used by the Latin hypercube
sampling in ``sensitivity.py``.  The default (20260816) reproduces every number
reported in the manuscript; override it without editing any file:

    DENGUE_SEED=42 python run_all.py
    python run_all.py --seed 42

Both forms resolve to the same value.  The resolution order is:

1. the ``--seed`` command-line flag of ``run_all.py`` (highest priority),
2. the ``DENGUE_SEED`` environment variable,
3. ``DEFAULT_SEED``.
"""

import os

DEFAULT_SEED = 20260816
ENV_VAR = "DENGUE_SEED"


def get_seed(default=DEFAULT_SEED):
    """Return the effective seed, honouring the DENGUE_SEED environment variable."""
    raw = os.environ.get(ENV_VAR)
    if raw is None or raw.strip() == "":
        return default
    try:
        return int(raw)
    except ValueError:
        raise ValueError(
            f"{ENV_VAR} must be an integer, got {raw!r}"
        ) from None


#: Effective seed, resolved once at import time.
SEED = get_seed()
