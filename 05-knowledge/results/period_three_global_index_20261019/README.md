# Period-three global index, resultant, and local delay

19 October 2026. A concise continuation of the parabolic/cotangent research.

The already known period-three iterative residue **51/49** follows without expanding the degree-eight return germ: the two ordinary fixed points contribute `6/49` to the global finite residue sum, leaving `−2/49` at each of the three parabolic cycle points.

The repeated **49** is connected by exact identities, rather than a numerical resemblance:

- collision resultant `Res(f−x, Φ₃)=16c²+4c+7`;
- at `c=−7/4`, `Φ₃=P²` and `Res(f−x,P)=7`;
- the cyclic gap product squared is `disc(P)=49`;
- the parameter unfolding product `a_r α_r=49` yields a proved local delay `N_return~π/(7√ε)` for `ε=c+7/4>0`.

The cyclotomic model, discriminant, and residue value are explicitly credited to preceding work. The contribution here is their global-index/resultant synthesis and the local unfolding-scale connection. An independent paired-cycle coalescence calculation checks the index. The note includes multiplicity, nonpolynomial, and coordinate-scaling boundaries. It makes no global critical-orbit escape, analytic embedding, or major-conjecture claim.

Read [the proof](period_three_global_index.md) and [independent review](review/REVIEW.md).

## Reproduce

Python 3 and SymPy are required. From this directory:

```sh
python3 verify_period_three_index.py
python3 -O verify_period_three_index.py
python3 review/audit_period_three_index.py
python3 -O review/audit_period_three_index.py
sha256sum -c SHA256SUMS
```

The primary verifier has 62 always-active exact checks. Normal and optimized outputs are byte-identical. The independent review supplies 32 separately derived exact controls and audits the local first-passage proof; finite computations are not substitutes for that limit proof.
