# Independent audit: conductor-eight classification and quartic/sextic exclusion

20 October 2026. **PASS.**

## Exact reviewed scope

Reviewed `CONDUCTOR_EIGHT_CLASSIFICATION.md`, Sections 1–7, and
`verify_frontier.py`: the three symmetric genus-four degree semigroups,
all nonprincipal characteristic patterns, explicit immersed polynomial
realizations, the complete monic degree-(4,6) delta-four normal form,
the actual monic Keller exclusion for the `<4,5,6>` sector, and the
actual-quadratic-coordinate lemma.

This audit does not review the separate all-degree hyperelliptic obstruction
or the later cubic/quintic addendum. In particular it does not promote a
complete conductor-eight Keller exclusion.

Reviewed source hashes:

```
e632f2e8bb25c43ce491bd09caf2503070175cf908e00b67e390627ec9dcfa70  CONDUCTOR_EIGHT_CLASSIFICATION.md
61b3ea3091cecdf7f76bf786ac98ae5ae1dcc1d10d2f49587f55acdddf5e077f  verify_frontier.py
```

The proof's earlier hash was
`5d6ef5d29e0857b353d4c1f45ae57a33b8cf9e81c181e849734b13fc1fac26af`.
The sole audit correction replaced the misleading description that four
Jacobian coefficients “successively” force four derivatives, in their listed
order, by the correct statement that they form a linear system. The actual
system is nonsingular with determinant `-192b`. No theorem or formula changed.

## Reproduction and independence

The producer checker passes all **57** tests in ordinary Python and with `-O`;
the outputs are byte-identical and match the staged receipts.

The reviewer-owned `independent_checks.py` does not import the producer checker.
It passes **42** checks ordinarily and with `-O`, with byte-identical outputs.
It independently:

- enumerates all 256 subsets below eight and imposes the semigroup/symmetry
  conditions;
- enumerates characteristic tuples from the entire bounded integer universe,
  rather than filtering first to the proposed semigroups;
- derives the successive remainder coefficients by reading and canceling
  leading terms;
- differentiates actual symbolic functions of x, builds the derivative
  coefficient matrix, computes its determinant, and solves it;
- verifies the whole Jacobian, source-centering chain rule, actual example
  immersions, algebra-recovery identities, and quadratic inverse controls.

Reproduction:

```
python3 verify_frontier.py > review/producer_replay.json
python3 -O verify_frontier.py > review/producer_replay_optimized.json
cmp review/producer_replay.json review/producer_replay_optimized.json
python3 review/independent_checks.py > review/independent_results.json
python3 -O review/independent_checks.py > review/independent_results_optimized.json
cmp review/independent_results.json review/independent_results_optimized.json
```

The exact computations corroborate the proofs below. No finite scan of
polynomial coefficients is used to prove an unrestricted Keller statement.

## 1. Symmetry, field extension, and characteristic applicability

The frozen low-delta packet proves symmetry by an explicit trace/residue
argument and identifies codimension with the number of missing degrees.
Its hypotheses exactly match a primitive two-generated subalgebra of k[t]
in characteristic zero. In genus four its Frobenius number is seven.
The paired membership relations under n ↦ 7−n, together with closure under
addition, leave precisely the gap sets

```
{1,3,5,7}, {1,2,4,7}, {1,2,3,7}.
```

The corresponding semigroups are `<2,9>`, `<3,5>`, `<4,5,6>`.
For example, once 2 is present, 4 is present, 3 and 5 are absent, and all
integers at least eight are present; if 2 is absent, 5 is present and
exactly one of 3,4 is present. This is an exhaustive proof.

The constant-field extension step is valid. For a field extension L/k,
flatness preserves the kernel of

```
A → k[t]/k[t]≤n.
```

Consequently `(A ∩ k[t]≤n) ⊗ L` is exactly
`L[p,q] ∩ L[t]≤n`, and every filtered dimension and occurring degree is
unchanged. Geometric integrality is automatic from the inclusion in L[t].
The identity expressing t rationally in p,q also survives. Normalization
is L[t] and the sole infinity point of its projective completion is rational.
Thus applying characteristic-zero one-place theory after algebraic closure
is legitimate; all triangular reductions used here remain over k.

The characteristic numerical axioms and bounds were freshly checked in
Sathaye–Stenerson, *On Plane Polynomial Curves*, Section 3.1, printed
pages 8–9, author-hosted text dated 10 September 2014:
https://www.ms.uky.edu/~sohum/ma561/affgeom/papnew.pdf .
Section 2.2 gives the one-place and residual-rationality convention.
The first two entries' identification with the actual coordinate pole
orders is inherited from the independently audited fixed cubic/quartic
packet, whose `PRIMARY_SOURCES.md` pins Abhyankar–Sathaye (1996),
Section 2, printed p.1062, DOI
https://doi.org/10.1090/S0002-9939-96-03254-6 .
This review does not claim a new retrieval of that second source.

## 2. Exhaustive patterns and polynomial realizability

For a nonprincipal initial pair the cited bounds are nine for the initial
entries and seven for later entries. Strict subsequent gcd drops and the
requirement to contain every primitive generator yield exactly:

```
(2,9), (9,2), (6,9,2), (9,6,2);
(3,5), (5,3);
(4,6,5), (6,4,5).
```

The written case split accounts for all nondividing pairs. A coprime pair
ends the sequence immediately; the noncoprime alternatives have gcd two or
three, so a later primitive entry cannot be skipped without ending the
sequence too soon or failing a strict drop. The first axiom's initial
semigroup membership condition is automatic because
`(δ0/gcd(δ0,δ1))δ1` is a multiple of δ0. The checker applies all remaining
conditions, not merely the gcd test. Reversing the first two entries is
individually checked.

The polynomial-realizability claim is not inferred from abstract
characteristic admissibility. The ordinary examples have coprime degrees.
The quartic/sextic example has a degree-five remainder and the full
four-element basis below. The sextic/nonic identity recovers t² and gives
the exact algebra

```
k[t², t((t²)^4+6(t²)^2+6)].
```

The even/odd decomposition proves its semigroup and codimension exactly.
Derivative gcds are one for all four examples. The nonic example's quartic
in z is squarefree and nonzero at zero, giving four ordinary nodes over C:
at a root α the two tangent slopes are opposite nonzero multiples of
`sqrt(α) p'(α)`. Swapping either polynomial pair realizes the opposite
orientation. No special extension of k is required for the displayed
polynomials or for the immersion identities.

## 3. Full-module argument and quartic/sextic normal form

If p has degree four and k(p,q)=k(t), then q has a monic degree-four
minimal polynomial over k(p) with coefficients in k[p]. Integrality of
q and integral closedness of k[p] justify the coefficient statement.
Therefore `1,q,q²,q³` is a free k[p]-basis of the actual algebra.

For every remainder of the form `r=q²-u(p)q-v(p)`, the transition to
`1,q,r,qr` is triangular with diagonal entries one. This is an actual
unimodular basis transformation, not just four independent elements.
If r has odd degree n, the four basis degrees are `0,6,n,n+6`, with
distinct residues modulo four. Thus no leading cancellation occurs across
the four polynomial-module summands. The missing-degree count is exactly

```
(0+6+n+n+6 - (0+1+2+3))/4 = (n+3)/2.
```

This proves the successive contradictions at n=11,9,7 and the exact
codimension four at n=5. The producer's cancellations have been independently
expanded and are correct. At D=0 with b≠0 the remaining cubic is
`b³ t(3a+2t²)/16`, giving codimension three. Its degree is genuinely three,
since characteristic is zero and b≠0. At b=0 both coordinates are even,
so the fraction field cannot be k(t). Conversely bD≠0 produces a degree-five
element, so the intermediate field degree divides both four and five and
is one. The same module argument then proves sufficiency.

Subtracting the scalar coefficients d,f is harmless for a k-algebra
comparison, but not for an actual (x,t)-Jacobian. The proof explicitly
restores and retains them. No x-dependent target subtraction is used in
the Keller argument.

## 4. The actual Keller argument

The constant-leading gate can be verified directly, without a further
classification input. If P is monic of positive y-degree m and Q has
positive degree n and leading coefficient l(x), the coefficient of
`y^(m+n−1)` in J(P,Q) is `−m l'(x)`. Thus l is constant. Scale by nonzero
constants to make both monic. A divisibility reduction uses a constant
target coefficient; after cancellation the untouched coordinate remains
monic, so the gate applies again to the new positive-degree coordinate.
Every step preserves the generic algebra and a nonzero constant Jacobian.
If a terminal coordinate is y-constant, the Jacobian product forces the
other y-derivative to be a nonzero constant, and the generic algebra is
k[y], contrary to delta four.

Centering the terminal monic quartic is a polynomial source shear of
determinant one. Its coefficients lie in C[x], and it leaves the actual
Jacobian unchanged. Applying the necessary generic normal-form identities
in C(x) is valid because they are equalities of already polynomial
coefficients. No generic generator is asserted to globalize.

At t=0, the whole Jacobian equals

```
b((c−3a/4)d'−f').
```

For nonzero constant J, this makes b a unit in C[x], so b∈C* and b'=0.
With unknowns `(c',d',e',f')`, the coefficients of t⁷,t⁵,t³,t² give matrix

```
[-4       0       0      0]
[-2a      6      -4      0]
[ 0       4c     -2a    -4]
[-2ab     9b/2   -b      0].
```

Its determinant is `−192b`. It is invertible even at `3a=4c`; no illicit
division by that possibly zero expression occurs. Its unique solution is
the displayed equation (E). Substitution into the entire Jacobian yields

```
J = −a'b(3a²−4ac+6bt+4e)/8.
```

The linear t coefficient forces a'=0; the whole expression then vanishes.
This is the required contradiction. All d',f' terms have been accounted
for, and constant-fiber delta or conductor monicity is nowhere assumed.

## 5. Actual quadratic coordinate

After completing the square by a source shear, p=t²+a(x). Evaluating J at
t=0 gives `a' q_t(x,0)=j`, so a'=α∈C*. Hence `(p,t)` is a global
polynomial coordinate system, with `x=(p−t²−β)/α`.
The equation at fixed p is `α q_t=j`, whose polynomial solutions are
`q=(j/α)t+F(p)`. This gives an explicit polynomial inverse and proves the
automorphism claim. The final target reduction is a C-coefficient shear.

This argument applies only to an actual quadratic member of the Keller
pair. The sextic/nonic example correctly prevents confusing it with the
existence of some quadratic element in the generic algebra.

