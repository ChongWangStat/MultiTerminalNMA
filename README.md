# MultiTerminalNMA — code and reproducibility materials

Companion reproducibility archive for **Multi-terminal study resistances in network meta-analysis: correlated arms, two random-walk conjectures, and non-unique evidence routes** (Wang, Jin, O’Connor; *Research Synthesis Methods*, Research in Brief, manuscript).

This is a **methodological paper**, not a clinical treatment-effect evaluation. The trial-inspired five-study network is illustrative.

## Contents

- `brief/code/`: GLS/current identities, route audit and sharp bounds, correlated-arm checks, depression-network diagnostics, and numerical consistency checks.
- `brief/figures/` and `brief/data/`: scripts and selected fixed inputs from the Research-in-Brief source.
- `psoriasis/`: the Python and R scripts used for the Cochrane-derived psoriasis network, plus summary results independently re-run in November 2026.
- `simulation/`: seeded structural-simulation generator and verified aggregate outputs (18,000 analyses).
- `REPRODUCIBILITY.md`: dependency and execution instructions, and known input-data restrictions.

## Quick start

Use Python 3 with NumPy, SciPy, pandas and matplotlib.

```bash
python -m pip install -r requirements.txt
python brief/code/verify_math.py
python brief/code/correlated_arm_examples.py
python simulation/simulation_structural_final.py --output-dir simulation/reproduced
```

The psoriasis study-level file `pso-psa-data.csv` is **not redistributed** because the project has not established permission to publish that underlying extraction. To rerun the psoriasis analysis, obtain the original extraction through its authorized source and place it in `psoriasis/` with that exact filename. See `REPRODUCIBILITY.md`. The provided result summaries can be inspected without the raw input.

## Earlier related work

- [Contrast-Space Projection for Network Meta-Analysis — source code](https://github.com/ChongWangStat/CSP-NMA).
- [CSP preprint](https://arxiv.org/abs/2604.21994).

## Reproducibility scope

The numerical checks validate the algebra and stored illustrative examples. They do not validate trial-level clinical inference. Prior third-party published example inputs retain their original attribution. See `PROVENANCE.md`.

License for this companion package should be confirmed by all authors before announcing a reusable distribution license. Previously released CSP code included here retains its original copyright and license; see `CSP_LICENSE_NOTICE.md`.
