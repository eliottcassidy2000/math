# Pre-registration, mac-mini-2026-11-09-genesis

Committed before the runs it describes. Already done before this file (exploratory, not pre-registered):
- the extension search over the seven 3-arc-strong order-11 non-creators (`ext12.py`): 12 escaper classes;
- two order-12 annealing runs (`sa12.c`, seeds 11 and 22);
- band certificates for band_3(13) and band_4(16) (`wit.c`, `verify_bands.py`);
- the first-creator class shares at orders 3-11 (`fc_share.py`).

Definitions:
- An **escaper** is an n-tournament that is none of the following:
  - a first creator (H not in S(n - 1));
  - a first creator plus a vertex;
  - non-strong;
  - non-prime;
  - lambda <= 2.

  This is THM-4942's five-operation grammar.
- The **decomposability defect** dec(T) is the fewest arc reversals reaching a non-strong or non-prime tournament,
  that is, min(lambda, module defect).

## Q1-Q2: order 13, annealing restricted to minimum semidegree >= 4 (`sa13.c`)

The energy uses S(11) exactly for the deletions and the proxy H(T) <= 545,835 = max H(12) for T. Every hit is
certified afterwards: H(T) in S(12) by an explicit witness (`wit.c`), and every deletion by klein's S(11).
Two seeds, at most 1.5M steps each.

- **Q1.** At least one order-13 escaper with lambda = 4 is found and certified. Prediction: YES (60%).
- **Q2.** At least one certified order-13 escaper has dec >= 4, so every module costs >= 4 reversals. Prediction:
  YES (45%).

## Q3: the spanning-tree contrast (`tau_share.c`)

Consider all labelled graphs on m vertices, m = 4..8. Call a connected graph **new** if tau(G) is not attained by any
connected graph on m - 1 vertices.

- **Q3.** The labelled share of new graphs among connected graphs decreases from m = 6 to m = 8. Prediction: YES
  (65%).

  For tournaments the class share of first creators increases at every order from 7 to 11. A decrease would make the
  two worlds opposite in this respect.

Outcomes will be appended below without editing the predictions.

## Outcomes

(to be appended)
