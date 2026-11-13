# Marked cubic structure, polynomial exchange and Sierpiński precision

See [PROOF.md](PROOF.md) for exact hypotheses and lineage.

1. **Complete relative boundary classifier:** with boundary incidences fixed, an interior tournament's full marked continuation matrix is equivalent to its augmented labelled cyclic-triple hypergraph. The new interface argument uses the corrected inherited THM-4803 C3/path-cover theorem. Cubic internal observations plus quadratic port-start observations suffice at every size; they do not reconstruct the tournament isomorphism class or give a fast evaluator.
2. **Polynomial pivot atlas:** destination occupancy annihilates the individual Schur denominator. Every two-pivot exchange is exactly five positive terms.
3. **Sharp prime-power inverse cutoff:** geometric powers truncate at the least K with v_p(K!)≥s modulo p^s, sharply even for tournament cycle denominators. Divided-power multiplication has the Lucas/Sierpiński carry pattern.
4. **Physical synthesis correction:** frozen source/sink guards forbid not only a single changed saturated seam but *any* boundary change. The marked reciprocity identity first forced a reverse response; the stronger module theorem then ruled out the overconservative paired repair. A constructive move must explicitly change and compensate a subset-count or low-degree port response.
5. **Triangle-group boundary:** (2,3,r) with r odd coprime3 acts trivially by global preserve/reverse symmetries of one tournament, explaining why (2,3,5)/(2,3,7) analogies require another carrier.

Reproduction from this directory:

```bash
python3 verify.py
python3 seam_probe.py
python3 paired_seam_probe.py
python3 rigidity_probe.py
```

All four also pass under `python3 -O`. The three physical probes each exhaust all labelled tournaments of orders3..6 (33,864 total), with fixed ports0,1. The rigidity probe checks both the degree3+degree2 observer and the augmented-C3 equivalence. `verify.py` independently checks path series against DFS, polynomial exchange, marked reciprocity, prime-power cutoffs and formal rank-one seam identities. The formal seam tests are not physical witnesses. Chan–Pak remains open.
