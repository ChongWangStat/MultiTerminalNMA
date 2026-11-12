# MultiTerminalNMA — code and reproducibility materials

Companion reproducibility archive for **Multi-terminal study resistance and electrical evidence flow in network meta-analysis** (Wang, Jin, O’Connor; *Research Synthesis Methods*, Research in Brief, manuscript).

This is a **methodological paper**, not a clinical treatment-effect evaluation. The trial-inspired five-study network is illustrative.

## Multi-terminal resistance framework

For a multi-arm study with a full-rank contrast matrix `C_k` and positive-definite covariance `Sigma_k` of its within-study contrasts, `Sigma_k` is a generalized resistance matrix in contrast coordinates. The corresponding terminal conductance is `K_k = C_k^T (Sigma_k)^(-1) C_k`. The electrical terminal law `i_k = K_k v_k` reproduces the fitted generalized least-squares coefficients. A study with correlated arms need not be realizable by ordinary positive resistors; the multi-terminal law still holds. The code is not a recommendation to allocate evidence by random walks.


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

Original software code is released under the **MIT License**; see [`LICENSE`](LICENSE) and [`LICENSE_SCOPE.md`](LICENSE_SCOPE.md). Reused CSP code retains its original authorship and MIT notice (see [`CSP_LICENSE_NOTICE.md`](CSP_LICENSE_NOTICE.md)). Externally sourced data and publications are not relicensed, and the restricted psoriasis extraction is excluded.

## Manuscript revision v19.7 (November 12, 2026)

The manuscript now makes three distinct NMA contributions explicit: (1) a covariance-aware multi-terminal electrical representation and criterion for ordinary positive-resistor realization; (2) proofs of two conjectures stated by Davies (2026) for independent arms; (3) identification conditions and sharp bounds for study-labeled routes. It distinguishes arbitrary within-study pairing masses from the uniquely defined maximal direct component of the CSP method, and states that complete-route decompositions require an acyclic study-terminal flow graph.

The five-study numerical illustration now uses anonymous study names **S1-S5** in its figures. Source CSVs retain their historical identifiers for exact numerical audit; no valid trial-specific clinical inference is claimed. The paired figure-rendering scripts in `brief/figures/` produce the revised labels. The supplementary counterexample is reproducible with `python3 brief/code/positive_resistor_cycle.py`; it demonstrates that study-terminal cycles can arise despite exclusively positive mesh conductances.

## Focused Research-in-Brief revision v19.8

The editorial review makes the three *NMA-specific* contributions visible: exact proofs of both Davies (2026) conjectures for independent-arm models; a covariance-based multi-terminal electrical interpretation for correlated multi-arm studies, with an exact criterion for realization by ordinary positive resistors; and route-identification conditions and sharp route bounds. Established GLS information matrices, electrical current identities, and contrast invariance are credited to earlier work rather than claimed as mathematical inventions.

A new numerical finding shows why the usual variance-triangle check alone is insufficient for a four-arm study. See [four_arm_triangle_counterexample.py](brief/code/four_arm_triangle_counterexample.py) and its [verified output](brief/audit_output/four_arm_triangle_counterexample.txt): every variance triangle inequality holds while one off-diagonal terminal conductance is **+0.05**, so an ordinary positive-resistor representation is impossible even though the covariance-based multi-terminal law remains valid.

The manuscript's five-study network is explicitly illustrative, not a clinical-effect analysis. The seeded 18,000-analysis structural simulation and psoriasis route-audit scripts remain in the repository; the externally sourced psoriasis raw arm data are not redistributed. Source and accompanying PDF files for v19.8 are provided in the author's submission package.
