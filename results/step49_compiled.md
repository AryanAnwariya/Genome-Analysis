# GM12878 results (n_init=100, lowest-inertia seed of 3, 26,226 Rao A1-B3 bins)

## INTERKR graph

| Method | Seed | Silh | DB | ARI(Rao) | NMI(Rao) | ARI(SCI) | NMI(SCI) | NMI(chr) | Largest cluster |
|---|---|---|---|---|---|---|---|---|---|
| LINE paper settings (25M, SCI merge) | 2 | 0.005 | 8.34 | 0.049 | 0.060 | 0.041 | 0.048 | 0.006 | 26% |
| LINE 1000M, per-order L2 | 3 | 0.076 | 3.33 | 0.443 | 0.487 | 0.268 | 0.311 | 0.192 | 30% |
| Weighted DeepWalk | 3 | 0.202 | 1.49 | 0.341 | 0.448 | 0.258 | 0.334 | 0.029 | 27% |
| HOPE (Katz) | 3 | 0.093 | 2.86 | 0.550 | 0.543 | 0.305 | 0.334 | 0.199 | 30% |
| DeepWalk + HOPE (L2 each) | 3 | 0.093 | 2.70 | 0.526 | 0.544 | 0.306 | 0.357 | 0.154 | 30% |

## Pair-O/E graph

| Method | Seed | Silh | DB | ARI(Rao) | NMI(Rao) | ARI(SCI) | NMI(SCI) | NMI(chr) | Largest cluster |
|---|---|---|---|---|---|---|---|---|---|
| LINE paper settings (25M, SCI merge) | 1 | -0.000 | 8.77 | 0.195 | 0.290 | 0.148 | 0.204 | 0.023 | 26% |
| LINE order 1 only (1000M) | 2 | 0.089 | 2.24 | 0.099 | 0.347 | 0.042 | 0.198 | 0.048 | 67% |
| LINE 1000M, per-order L2 | 3 | 0.058 | 3.46 | 0.551 | 0.560 | 0.352 | 0.388 | 0.060 | 37% |
| Weighted DeepWalk | 1 | 0.274 | 1.30 | 0.428 | 0.527 | 0.324 | 0.380 | 0.028 | 37% |
| HOPE (Katz) | 1 | 0.091 | 2.69 | 0.616 | 0.624 | 0.350 | 0.407 | 0.052 | 26% |
| DeepWalk + HOPE (L2 each) | 1 | 0.104 | 2.57 | 0.633 | 0.646 | 0.362 | 0.421 | 0.049 | 27% |

Reference: published SCI vs Rao ARI 0.503, NMI 0.567. NMI(chr) reference: Rao 0.100, published SCI 0.035.

Caveats:
1. Agreement with Rao is not the same as biologically better; biological validation is pending.
2. DW+HOPE and pair-O/E were informed by inspecting earlier results; each was pre-stated and tested once.
3. LINE sample count and merge were tuned after seeing scores; paper-settings LINE is the main LINE result.
4. Own choices: INTERKR normalisation, HOPE beta and dimension, DeepWalk parameters.
5. On INTERKR, HOPE, DW+HOPE and LINE 1000M have chromosome NMI 1.5-2x Rao's 0.100.
6. HOPE and DW+HOPE are roughly tied; no single best method is claimed.
7. LINE paper settings on INTERKR contains a 16-node outlier cluster; LINE order-1-only on pair-O/E has one cluster holding 67% of bins.
