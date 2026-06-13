---
id: HYP-1901
status: OPEN
source: codex-2026-05-31-S450
related:
  - THM-357
  - THM-360
  - THM-369
  - HYP-1794
  - HYP-1802
  - HYP-1813
  - HYP-1836
  - HYP-1868
  - HYP-1894
  - HYP-1890
  - HYP-1891
  - HYP-1900
  - HYP-1902
---

# HYP-1901: LRC analogues factor through protected anti-Bohr boundary incidence

## Statement

Every exact formulation of the Lonely Runner Conjecture factors through the
same finite boundary-incidence object.

Start with a one-parameter subgroup

```text
t -> t*v in (R/Z)^k
```

and forbidden neighborhoods of the coordinate subtori

```text
||v_i*t|| < 1/(k+1).
```

The exact proof problem is not primarily about the speeds, the volume of the
bad set, or the chromatic number of a distance graph.  It is about the finite
protected boundary of this anti-Bohr cover:

```text
endpoints = equality states ||v_i*t|| = 1/(k+1)
protection = strict coverage of such an endpoint by another forbidden arc
counterexample = full open cover with every endpoint protected
```

Distance-graph regular circular colorings, subtorus/cube avoidance,
view-obstruction, zonotope covering radius, lonely-runner spectra, and the
repo's endpoint/IP/Bruhat-Tits ledgers are useful precisely to the extent that
they preserve this boundary incidence.  This refines HYP-1900's labelled
incidence-core claim by making the protected anti-Bohr boundary the object
that any exact analogue must transport.

## Evidence

`lrc_analogy_atlas_s450.py` checks the shared finite-boundary object on exact
examples:

```text
initial n=14:
  boundary_only, forbidden=1, unprotected=6, first_unprotected=1/14

n14 seven-ladder:
  positive_gap, forbidden=142/143, max_gap/th=5/924,
  unprotected=84, first_unprotected=9/98

n14 S380 gate ladder:
  positive_gap, forbidden=142/143, max_gap/th=5/1848,
  unprotected=168, first_unprotected=15/196
```

The same examples can be read as:

```text
failed multiplier colorings in a distance graph,
line/cube intersections in a torus,
view-obstruction residues,
zonotope covering debt,
IP endpoint rows,
or p-adic frontier mass.
```

But the invariant that remains visible across all these readings is the
boundary debt: positive gap or unprotected endpoint.  The near-counterexample
ladders shrink real gaps while increasing exposed endpoint debt.

The 2026 web pass also changes the strategic frontier.  Current
computer-assisted preprints prove the conjecture through the title range of
thirteen runners, making `n=14` the first frontier denominator if those
preprints are accepted as the working state.  That lines up with the repo's
long emphasis on `n=14` as the first mixed product-building case.

## Predictions

1. Any analogy that cannot identify the image of an unprotected endpoint will
   not produce a durable proof invariant.
2. The useful functor from zonotopes to the repo should send endpoint debt to
   covering-radius row debt or facet-normal debt.
3. The useful functor from distance graphs to the repo should send failed
   multiplier colorings to the same endpoint-protection hypergraph, not merely
   to ordinary chromatic number bounds.
4. The useful functor from spectra to the repo should send the recursion
   `accumulation(S(n)) = S(n-1)` to runner deletion / endpoint deletion.
5. For `n=14`, a serious certificate should be readable in at least three
   equivalent languages: endpoint peel, IP dual rows, and zonotope covering
   debt.
6. A Zeckendorf/Ostrowski analogue is useful only if it gives a normal form for
   endpoint debt, not merely a Fibonacci label for the final numbers.

## Sources

Repo:

- `04-computation/lrc_analogy_atlas_s450.py`
- `04-computation/lrc_zeckendorf_bridge_s451.py`
- `05-knowledge/results/lrc_analogy_atlas_s450.out`
- `05-knowledge/results/lrc_zeckendorf_bridge_s451.out`
- `07-reflections/lrc-analogy-atlas-s450.md`
- `07-reflections/lrc-zeckendorf-bridge-s451.md`
- THM-357.
- THM-360.
- HYP-1836.
- HYP-1890.
- HYP-1891.
- HYP-1894.
- HYP-1900.
- HYP-1902.

Web:

- Rosenfeld, eight runners: https://arxiv.org/abs/2509.14111
- Trakulthongchai, nine and ten: https://arxiv.org/abs/2511.22427
- Sungkawichai-Trakulthongchai, eleven through thirteen: https://arxiv.org/abs/2604.23906
- Barajas-Serra, seven runners and regular chromatic number: https://arxiv.org/abs/0710.4495
- Henze-Malikiosis, zonotope covering radius: https://arxiv.org/abs/1609.01939
- Cambridge finite-checking/zonotope restatement: https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/linearly-exponential-checking-is-enough-for-the-lonely-runner-conjecture-and-some-of-its-variants/A51A991DE89B8C9C2E2FF13FBD4501DA
- Giri-Kravitz, spectra: https://arxiv.org/abs/2304.01462
