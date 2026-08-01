---
source: kind-pasteur-2026-07-24-S161 (Opus 4.8)
status: RESULT (exact, closes the isolated gap at D=2). Searched leak-Jacobian rank drops via the discriminant.
  Key clarification: rank drops (det J=0) locate counterexample FAMILIES (which generic transversality already
  disfavors); the ISOLATED counterexample (the KZ gap) has FULL rank and is NOT on the rank-drop locus -- it is
  found only by testing emptiness of the full leak IDEAL. In the DFT coordinates the FC(3) leaks are
  integer-coefficient polynomials; for D=2, Groebner({leak1,leak2,leak3}) = (1), so {all leaks=0} is EMPTY: NO
  FC(3) counterexample at D=2, isolated OR family. This is an EXACT no-counterexample proof (cyclic-weight family,
  D=2), strictly stronger than the generic rank scan. Explicit rank-drop discriminant det J recorded. D=3 ideal
  test attempted.
tags: [factorial-conjecture, discriminant, groebner, transversality, exact-proof, roots-of-unity]
related: [kps-S157, kps-S159, kps-S160]
---

# FC(3): the discriminant, and an exact no-counterexample proof at D=2

## 1. Geometry -- what the discriminant actually finds (honest)
Leak map `Phi(c)=(L(f^3),L(f^6),...)`, `P` params, generic rank `P` (transversal, D<=5, kps-S159).
- **Rank-drop locus** `{det J_{PxP}(c)=0}`: where `Phi` fails to be an immersion -- where a **positive-dimensional
  family** of counterexamples could live. Generic transversality already says this is a proper subvariety.
- **Isolated counterexample** (the KZ gap): has **FULL** rank at `c*` (it is a transversal isolated zero of the
  overdetermined leak system), so it is **NOT on the rank-drop locus.** Rank-drop search cannot find it.
> **The definitive test is emptiness of the leak IDEAL** `I=(leak1,...,leak_{P+1})`: `V(I)=varnothing` iff `1 in I`,
> which rules out **all** counterexamples (isolated and family). That -- not the discriminant -- closes the gap.

## 2. The leaks are integer polynomials; D=2 discriminant
In DFT/`(A,P,Q)` coordinates the moment table `T(A,P,Q)=A!P!Q!*g(A,P,Q)` is an integer and the params are free, so
each `L(f^{3j})` is a polynomial in `c` over `Z`. For D=2 (`f=L1+c1 Lbar1^2+c2 s L1`):
```
leak1 = 120 c1^3+576 c1^2 c2+72 c1^2+1008 c1 c2^2+252 c1 c2+18 c1+336 c2^3+126 c2^2+18 c2+1   (deg 3)
leak2 = 665280 c1^6 + ... + 1                                                                  (deg 6)
```
**Rank-drop discriminant** (det of the 2x2 leak-Jacobian on leaks 1,2):
```
det J = -209952 * c1 * ( 2240 c1^6 + 23520 c1^5 c2 + ... + c2^2 )   [c1 * irreducible sextic]
```
So rank drops occur exactly on `{c1=0} cup {sextic=0}` -- a curve. (Where a family *would* live, if one existed.)

## 3. EXACT no-counterexample at D=2 (the real result)
Groebner basis (grevlex, exact over `Q`):
> **`Groebner( leak1, leak2, leak3 ) = (1)`.**
Hence `V(leak1,leak2,leak3)=varnothing`: **no `c` kills even the first three leaks.** A fortiori no `c` kills all
leaks. **There is NO FC(3) counterexample in the cyclic-weight family at D=2 -- isolated or family.** This is an
*exact proof*, and it **closes the isolated (KZ) gap for D=2** -- which the generic rank scan could not do.

Note the rank-drop locus `{det J=0}` is a proper curve and does **not** contain a common zero of the leaks (since
there is none at all); the counterexample, had it existed, would have been the isolated full-rank point the
Groebner test rules out.

## 4. Scope and D>=3 (honest)
- **Family scope:** the cyclic-weight (`C3`-eigenvector) family is the **most favorable** place for a counterexample
  -- there the `3 !| k` moments vanish for free, so only the `3|k` conditions remain. Emptiness *there* is strong
  evidence; a fully general degree-2 `f` (all `C3`-weights, all-`k` conditions -- strictly harder) is a larger
  ideal, a further computation.
- **D>=3:** the ideal test is the definitive tool but the leaks have degree `3j` in `P` variables; the Groebner
  cost grows fast. D=3 (`P=5`, leaks 1..6, degrees up to 18) attempted via **modular** Groebner (mod `2^31-1`);
  result pending / cost-bounded. If it returns `(1)`, D=3 is closed mod `p` (isolated gap too). The generic rank
  scan already gives D<=5 no-family + no-generic-obstruction; the ideal test upgrades chosen degrees to
  no-isolated-either.

## 5. Status
D=2: **exact, definitive** -- no FC(3) counterexample (cyclic-weight), isolated or family; discriminant explicit.
This is the first result to close the isolated-point gap at any degree. The method (integer leaks -> ideal
emptiness, computable mod `p`) is the right definitive tool; scaling it to D>=3 is the concrete frontier.

Files: `/tmp/{disc_d2,gb_d3}.py`. Builds on kps-S157/S159/S160.
