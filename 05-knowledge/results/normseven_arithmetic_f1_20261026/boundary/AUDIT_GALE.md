# Independent audit of the outer Gale twist classification

2026-10-26. Auditor: parallel `boundary_arithmetic` agent, independent of
the Gale producer. **PASS** for the statements of Theorem GA in
[`gale/RESULT.md`](../gale/RESULT.md), with the dependencies and limits below.
No theorem-level defect found. This is a written proof audit and an exact
independent arithmetic path, not a formal theorem-prover certification.

## Proof review

1. **Correct surface.** The equation `w^2=E(E^2-x^2)(E^2-z^2)(E+x+z)`
   is precisely the double-plane K3 called `Z` in the inherited exact-rank
   corollary. Its rank19 is not being transferred without proof from the
   different degree-eight root-intersection model. I read the full
   `outer-picard/EXACT_RANK_COROLLARY.md` and the branch-resolution and
   point-count parts of `picard-upper/INDEPENDENT_UPPER_BOUND.md`.
2. **Eighteen rational invariant classes.** The nine A1 and two D4
   exceptional configurations plus a plane pullback have nondegenerate
   rank18. At a branch node the tangent exceptional conic is
   `w^2=cuv`, so it has the rational points at `u=0` or `v=0`, regardless
   of the constant twist. This is unlike the anisotropic exceptional
   conics in the sibling central-boundary result. At a triple the central
   exceptional component is in the branch; each arm is attached to a
   different rational crossing. Thus every class is defined over Q and
   individually fixed by the deck involution. Counting components does not
   silently assert that the involution fixes their points.
3. **The nineteenth class.** The rational section `P=(dt,d^2)` is inherited
   nontorsion. The class `[P]-[-P]` maps to `2P` in `NS/Triv`, so is
   nonzero. It is anti-invariant under deck negation. It is independent of
   the invariant rank18 subspace; geometric rank19 then exhausts NS over Q.
   Primitivity of P and saturation of the displayed lattice are unnecessary.
4. **Descent and abstract nonisomorphism.** The constant-field comparison
   `w -> w/sqrt(delta)` has cocycle the deck involution. It twists only the
   one-dimensional anti-invariant class, yielding `1^18 + chi_delta`.
   Any abstract k-isomorphism is Galois-equivariant on NS; character traces
   force equality of the quadratic characters. This proves necessity even
   when an isomorphism need not preserve the fibration or plane map.
   Scaling `w` gives the converse. The representation is a rational
   representation, not an integral lattice classification.
5. **Sections and torsion.** The same character action applies to the
   one-dimensional geometric Mordell-Weil free space. A nonsquare constant
   twist therefore has no rational free class. The inherited I0* fibre
   bounds geometric torsion by `(C2)^2`, and the three visible nonzero
   two-torsion sections attain the bound. This keeps k(t) fixed and makes
   no claim about specialized fibre ranks or all rational points on Z.
6. **Counts and Gale sign.** The raw count is
   `q^2+q+1+chi(delta)S_q`; nine node corrections and two D4 corrections
   add `17q`. The `-4` Gale multiplier gives the actual minus-one twist
   after dividing the cover coordinate by2. Equality over a split prime
   or an even extension is only a control, not a rational isomorphism test.

I independently opened [Schuett--Shioda, Elliptic Surfaces, v3](https://arxiv.org/pdf/0907.0298)
and checked Theorem6.3 (printed p23), Lemma7.3 (p32) and Lemma7.8 (p35).
They provide the divisor/section quotient and additive-fibre torsion bound
used here. The required singular fibre and characteristic-zero hypotheses
hold. The inherited rank19 proof remains a named dependency; this audit did
not rerun its Tate/Artin--Tate proof from scratch.

## Independent exact path

`audit_gale.py` imports no producer code. It literally enumerates `(x,z,w)`
over F5,F7,F11,F13,F17,F19 for **every** nonzero twist, separately counts
every exceptional node conic, and checks the twist mean. It recovers
`180/172` at7 and `408/408` at13 without Legendre-symbol sums. It also
computes eight rational duplications on the inherited t=2 specialization,
recovering valuations `-2,-4,...,-16`, and replays the producer checker
normally and with `-O` in a temporary directory. Its independent path is
the literal cubic enumeration; replay is labelled separately.

Run from repository root:

```sh
python3 -B 05-knowledge/results/normseven_arithmetic_f1_20261026/boundary/audit_gale.py
python3 -B -O 05-knowledge/results/normseven_arithmetic_f1_20261026/boundary/audit_gale.py
```

The saved normal and optimized audit receipts are byte-identical,
SHA-256 `d421c6c3ed87c837a9f87c247f983df5cf4f53f7ce2e54b6ab100105511f8f2b`.
The computations corroborate the arithmetic consequences. They do not
compute the NS representation independently of the written divisor proof.
