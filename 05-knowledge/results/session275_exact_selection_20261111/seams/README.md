# Exact seam selection lane

- [Proof](PROOF.md): bifan curvature, exact divisor inverse, hidden root menu,
  renewable +2 sibling bridge, and **balanced-port transport of an actual
  single-free-edge 61→63 portal to all scales**.
- [Controls](controls.json): declared finite checks, not a spectrum extension.

Reproduce from any working directory:

```sh
python3 /path/to/repo/05-knowledge/results/session275_exact_selection_20261111/seams/check.py
```

`explore.py` supplies direct subset Hamilton counts and the physical graph
constructor. `bifan.py` supplies the root-packet computation; `check.py`
recomputes special packets with permutation enumeration, tests actual
extended graphs, and checks the divisor inverse.

The new all-order claims follow from the explicit source/sink laws, the
bilinear factorization, and the linear conserved charges. Controls cover
1,140 actual bifan graphs, all 5,405 rooted labelled inputs through order5
for the *false* conjectural curvature floor, and independent explicit
counterexamples and portal graphs up to order14. A separate optional random
search for additional one-sided twin roots was interrupted after a suitable
exact balanced-fan bridge was found; no completeness claim rests on it.

The low-jump portal source is supplied by this session's landscape lane.
No all-order connection between successive translated corridors has been
proved. The same-root fan language remains sparse, despite its renewable
interfaces and exact inverse.
