# Independent threshold audit — first checkpoint gate

**Verdict: PASS.** No blocking error found in `ROOT_THRESHOLDS.md` G1–G3,
section ladders, the surface threshold dictionary, or the local-cone sidecar.
Audited 2026-10-25.

## Proof-level checks

- **Degree and integrality:** the maximal root cover has degree `r^(n-1)`,
  while the complete-intersection degree is `r^(n-d-1)`. Divisorial valuations
  show full Kummer independence even for composite `r`; finite flatness then
  excludes additional vertical components or nilpotents.
- **Smoothness and canonical balance:** general position makes every local
  branch collection part of a parameter system. Root replacement is smooth;
  adjunction gives `a=(r-1)n-r(d+1)`. The `a<0` case is not asserted rational.
- **Sections:** `x_j=0` restricts to `X(d-1,n-1;r)` under general position,
  preserving codimension and raising `a` by one. The warning about special
  arrangement sections is necessary and correctly present.
- **Surface lists:** `(r-1)(n-3)=3` or `4` gives exactly the stated two K3
  and three canonically embedded `a=1` uniform-degree cases.
- **Local cone:** write `A=r(m-2)-m=(r-1)(m-2)-2`. For `m>=3,r>=2`,
  `A>=-1`, with `A=-1` only at `(m,r)=(3,2)`. Hence the unique genus-zero
  singular cone is the quadratic conic cone (A1). The boundary `A=0` occurs
  at `(m,r)=(4,2)` or `(3,3)`, both elliptic, discrepancy `-1`. All other
  cases have positive genus and negative discrepancy. In particular a triple
  point at exponent three must NOT be treated as a rational double point.
- **Normality:** distinct spanning branch lines leave only isolated cone
  singularities. A complete intersection is Cohen–Macaulay, so S2 plus R1
  gives normality. Spanning ensures `m<=n-1`, so the stated vertex count
  `r^(n-m-1)` is meaningful. The blow-up resolves each cone and has no
  exceptional rational `(-1)` curve.
- **All-incidence obstruction:** projection to concurrent root coordinates
  is dominant onto the irreducible pencil root curve. Its power-map image
  is the pencil `P^1`, so it cannot be constant. A rational surface cannot
  dominate a positive-genus curve. If no such point exists, only A1 nodes
  occur (`r=2`), or the surface is smooth (`r>=3`), and nonzero geometric
  genus follows from `a>=0`. This proves G3 for all stated exponents.
- **Torification:** a finite torus stratification of an irreducible variety
  has a full-dimensional locally closed torus, necessarily dense open.
  Thus nonrationality blocks the stratification and every nonempty open.
- **Hostile cones and numerical sidecar:** the negative discrepancy correctly
  prevents promoting ample singular dualizing class to general type of a
  resolution. The square-root Euler and chi corrections agree algebraically
  with the topological stratification; chi is not silently identified with
  geometric genus.

## Independent executable checks

`check_threshold_independent.py` does not import the source checker. It compares:

1. The character-eigensheaf sum for `H^d(O_X)` against the complete-intersection
   Hilbert coefficient in **350 cases** (`1<=d<=5`, `2<=r<=8`,
   `d+1<=n<=2d+7`).
2. The Euler characteristic counted over open branch strata against the
   complete-intersection Chern class in the same cases.
3. Riemann–Hurwitz genus and exceptional adjunction discrepancy for **342
   local cones** (`2<=r<=20`, `3<=m<=20`).
4. Both surface Diophantine lists through exponent/root count 100.

All pass in normal and optimized Python; audit gates use always-active explicit
checks rather than assertions. The source `check_thresholds.py` was also rerun successfully
(1,500 Hilbert comparisons plus matrix/section/cone controls). These finite
checks supplement, not replace, the proof above.

## Primary-source verification

- [Hidalgo–Hughes–Leyton-Álvarez, v8](https://arxiv.org/html/2010.04628v8),
  Proposition 1 and §2.4 Proposition 2: the smooth complete intersection and
  canonical coefficient agree after `k=r`, literature `n=our n-1`.
  Their displayed plurigenus asymptotic needs the positivity qualification;
  the packet explicitly avoids importing it outside that range.
- [Billi–Palomino–Sarti, v1](https://arxiv.org/html/2610.21752v1),
  Proposition 3.6: the CY condition is exactly the claimed balance. The
  date 2026-10-15 and preprint identifier were verified, not guessed.
- [Catanese, §4.3, printed p.47](https://www.mathe8.uni-bayreuth.de/de/team/prof-fabrizio-catanese-old/pdf/156.pdf):
  the plane square-root sequence and distinction from one double cover are
  explicitly present.
- [Hirzebruch, primary scan](https://hirzebruch.mpim-bonn.mpg.de/244/1/hirzebrucharr.pdf):
  the introduction uses these full Kummer function-field extensions and
  excludes all lines in one pencil, consistent with the spanning hypothesis.

No change to the mathematical claims is required for the first checkpoint.
