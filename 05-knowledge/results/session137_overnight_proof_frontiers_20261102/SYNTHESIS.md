# Session137 synthesis: closure, local universality, and exact resource barriers

## Result map

The session closes the global determinant lower-bound problem and strengthens it quantitatively. On the spectrum side, it establishes universality at every modulus while exposing why congruences, convex relaxations, finite edit atlases, and bounded construction resources do not settle exact realization. The harmonic lane proves a complete closed-boundary theorem for a specified family without claiming LRC(14).

Here **PROVED** means an all-order argument within the stated domain. **FINITE-EXACT** means a complete finite universe with its reduction and exact computation recorded. Independent review and replay scopes are listed in [VERIFICATION.md](VERIFICATION.md).

## 1. The determinant target is closed

For a nonempty tournament T, set

```text
u=det(I+A),  v=(-1)^n det A,
D=disc(T)=det(I+A-A^T)/2^(n-1),
s=1^T(I+A-A^T)^(-1)1.
```

Skew-matrix positivity gives `u+v=D>0` and `u-v=sD>0`. Exact vertex-cycle recurrences, together with induction on the weaker assertion `H>=u`, then prove

```text
H >= u+|v|+C_proper >= max(1,s)D >= D.
```

The weaker induction has a valid empty base; the displayed strong bound does not. This is an all-order proof, independent of any strong-tournament census or feedback-number cutoff.

The equality classifications are exact:

- `H=u+|v|` iff T is transitive or the directed triangle.
- `H=D` iff T is transitive.
- `H=u+|v|+C_proper` iff no two directed triangles are vertex-disjoint.
- Two disjoint directed triangles force the proper-cycle surplus to be at least two.

For weighted tournament matrices the bound `H>=u+|v|` remains valid. At orders one and two equality always holds; at higher orders equality requires a deterministic transitive tournament or deterministic directed triangle. The integer proper-cycle refinement is not exported to fractional matrices.

The companion counting theorem places precisely these four functions in #P on nonempty tournament inputs:

```text
H-u,                    H-u-C_all,
H-u-|v|,                H-u-|v|-C_proper.
```

Canonical certificates have length `O(n^2 log n)` and polynomial-time verification. Extending the last two functions to empty input requires an explicit value zero, not evaluation of the subtraction formula. No closure of arbitrary nonnegative GapP functions is used.

[Global proof](discriminant/GLOBAL_PROOF.md) · [Counting theorem](global_referee/COUNTING_GAPS.md)

### Corrections retained

The exact SCC discriminant product law remains valid, but unconditional discriminant supermultiplicativity was false: a strong factor with `disc=6`, `s=2/3` has ordinal-square discriminant 26, below 36. Reducibility also does not force `s<1`. The new proof uses neither error. The superseded pre-closure draft `discriminant/THEORY.md` is not a current proof source.

## 2. Quantitative weighted stability

Let `V=sum_(i<j) A_ij A_ji`, `delta=H-u`, and `Delta=H-u-|v|`. Augmented-skew principal minors yield the cofactor bounds; averaging deletion-root inequalities gives

```text
delta >= c_n V,       n>=2,
c_n = [24-2^(3-n)(n^2+2n+3)]/[n(n-1)],

Delta >= beta_n V,    n>=3,
beta_n = [24-2^(4-n)(n^2-n+3)]/[n(n-1)]
       >= 6/[n(n-1)].
```

Both coefficients are asymptotic to `24/n^2`. At order three `Delta=V` identically. Two symmetric half-weighted pairs in ordinal sum give `c_4=7/8` and `beta_4=3/4`, sharply. This session makes no higher-order sharpness claim. Order two is a necessary exclusion from the strong stability domain: there `Delta=0` even with fractional edges.

For the unordered-edge rounding distance `d=sum_(i<j) min(A_ij,1-A_ij)`,

```text
d <= 2V <= 2Delta/beta_n <= n(n-1)Delta/3,    n>=3.
```

Full off-diagonal matrix L1 distance is twice d. This bounds distance to a deterministic tournament only; it does not identify a transitive or optimal isomorphism class.

[Stability manuscript](discriminant/WEIGHTED_STABILITY.md) · [Polynomial audit](independent/WEIGHTED_POLYNOMIAL_AUDIT.md) · [Earlier cofactor audit](independent/WEIGHTED_STABILITY_AUDIT.md)

## 3. Modular information is exact, but exact values need more

For `p_T(c)`, the number of unordered covers by c nonempty directed paths, every homogeneous one-module context has c-th coefficient divisible by

```text
d_c=c!(c-1)!.
```

An actual triangular detector attains this divisor. Therefore two nonempty tournaments agree in every context modulo M **iff** their coordinates agree modulo `M/gcd(M,d_c)`. Only the finite initial segment before `M|d_c` matters. Representative-free positive substitution transitions preserve these states. The exact universal contextual class counts at moduli 3,4,8,12 are respectively 5,3,14,11, with width-three witnesses for every class.

The narrower grammar of ordinal joins and unary triangle wraps misses H residues 7 and 11 modulo 12. Genuine multichild triangle branching produces the infinite escape families `11*25^k` and `55*25^k`. This is a separation between constructor grammars, not a proof of globally unbounded spectral width.

The stronger phase-two theorem uses K1, ordinal joins, and `B(T,U)=C3[T,U,K1]`. For every modulus M it realizes every residue compatible with odd H, by strong tournaments with arbitrarily large distinct exact H values. Thus the strong scalar profinite closure is the odd profinite integers. Local spectral width is three for every positive odd target other than one, whose positive-index local width is one.

The construction synchronizes CRT conditions using actual positive contexts. It does not infer simultaneous coverage from separate prime-power surjectivity. The mod-5 fixed diagonal, singular operator at 13, and restrictions on stripping final padding remain explicit safeguards.

**Local does not mean exact.** The inherited exact values `w(13)=w(19)=5` and `w(31)=6` are unchanged. In particular local width-three universality does not fill those exact width-three holes, or imply full spectrum completeness.

[Exact context proof](width/PROOF.md) · [All-modulus proof](width_phase2/PROOF.md) · [Independent positive-context audit](independent/WIDTH_UNIVERSALITY_AUDIT.md)

## 4. Actual constructors: saturated relaxations, thin bounded banks

For an actual marked source, the four-state vector is `v=(C,L,g,t)` and parent height is `h=g+t`. Actual states generate the full integer lattice, their closed conical hull is the nonnegative orthant, and their closed convex hull is exactly

```text
{v>=0 : g+t>=1}.
```

Consequently a linear operation that sends every actual source to an actual source has a nonnegative integral matrix. The converse is not asserted: lattice, parity, positivity, and convex-hull membership are not realization criteria.

The original marked-port cost survives every allowed right-local migration history. Above the original cost threshold, a bounded output must seal to an attained multiple of the original F. Thus H-neutral padding can destroy all nonmultiple lifts rather than preserve the usable resource. An exact cap-173 overflow quotient has 2,592 states, compared with 6,234 positive-saturation states.

Fix both height h and interior L. Every source has a finite actual base and an optional zero tail. Base parent order is at most

```text
h+1                         if L=0,
h+3+floor(log2 L)            if L>=1.
```

For a positive tail of length k,

```text
F=2^k q0+h-2L,  q=2^k q0,  g=0,      k>=1.
```

The k=0 base retains its original value separately. Hence every bounded `(h,L)` bank produces only `O(log X)` distinct outputs up to X. Unbounded resources are necessary for a universal constructor, not proved sufficient. The exact `h=L=1` language is `{1,3} union {3*2^k-1,5*2^k-1:k>=1}`.

[Source geometry and sealing](constructors/THEORY.md) · [Finite-base theorem](constructors/BOUNDED_INTERIOR.md) · [Independent review](independent/BOUNDED_INTERIOR_AUDIT.md)

## 5. Natural cycle carriers: a complete budget classification

A natural carrier has one vertex for each literal directed odd cycle, with adjacency recording intersection. Different cyclic orders on the same support remain distinct. For a connected natural carrier of packing number k, the sharp minimum cycle-vertex count is

```text
floor((3k-1)^2/4).
```

For a fixed number of components the exact minimum is the balanced sum of these costs. A connected m-cycle carrier has a strong source of order at most `1+floor(sqrt(4m+1))`, giving an exact finite recognition algorithm, not a polynomial-time claim.

This source bound makes the **29 connected carrier types through cycle-budget 15** complete over all source orders. Their disjoint-union compiler gives exactly **300 H values at total cycle cost at most 15**, including a minimum-cost source expression for each value. Closed-antineighborhoods of literal cycles are actual remainder carriers, so the link double count and forbidden disjointness degrees retain source realizability.

For H=31, minimum cycle cost is nine and is attained at source order seven. Minimum source order is six, but every order-six realization uses at least thirteen odd cycles. The two resource minima cannot coincide; neither should be identified with substitution width.

[Carrier proof and complete table](cycle_carriers/THEORY.md) · [Verification](cycle_carriers/VERIFICATION.md)

## 6. Coverage, edit barriers, and the remaining scalar equation

The common-strong-child overlap graph on strong n-tournament classes has radius at most `n-4` for `n>=4`. The exact common-core two-ear seam gives reciprocity, a lattice, and unbounded fixed-core divisibility. These prove connectivity and constraints on range hulls, not solidity of scalar intervals. The order-six physical unit-atlas obstruction has an explicit nonstrong-parent repair.

Complete order-seven replacement atlases expose an H147 deletion-deck orbit closed under zero/unit star moves; reaching the neighboring H145 level requires at least three arc reversals, with four minimum repairs. A separate literal leaf audit does not depend on the complete class census.

For `T=Q[C3,...,C3]` with k blocks, altering arcs incident to r selected vertices gives

```text
2*3^max(k-3r,0) | H(T')-H(T).
```

For q individual arc reversals the sharper exponent is `max(k-2q,0)`. These are all-order marked-boundary divisibility theorems, but do not rule out escaping by a sequence of zero-response moves.

Deleted-triangle cards have a unique homogeneous pair and canonical completion within the triangle-substitution family, with no strongness hypothesis needed for this structural result. Scalar H-restoring-cut uniqueness is a different question. It is verified through all strong quotients of order seven, but remains **OPEN** at arbitrary order.

The complete k=7 producer census is **FINITE-EXACT**: 353 quotients, 2,471 deletion-orbit rows, 40,484,864 profiles, and 2,591,031,296 represented cuts, with zero nonoriginal H-restoring cuts. Independent verification covers a complete selected k=6 literal-cut row, not the whole k=7 census. The selected Paley7 forward-DP coefficient reconstruction is additional producer cross-route evidence. Repeating the full producer program does not change this independence classification.

A formal positive capacity perturbation preserves the listed symmetries, divisibility, deletion floors, and layer means while creating additional zeros. It is not a tournament counterexample; it shows why those coarse invariants alone cannot prove scalar zero-rigidity.

[Coverage results](coverage/RESULT.md) · [H147 atlas](coverage/PHASE2.md) · [Edit theorem and finite zero census](coverage/SYMMETRY_EDIT.md) · [Card rigidity](coverage/TRIANGLE_CARD_RIGIDITY.md) · [Audit scope](coverage/AUDIT_STATUS.md)

## 7. Harmonic first-missing-speed chambers, including ties

For `C={1,...,q-1} union {q+b:b in A}`, with `A subset {1,...,q-1}`, let d be the first missing second-band index, or q if none is missing. The new chamber extends to `r0=1/(q+d)`. Writing W for the explicit primitive q-window atlas and `P_(q+d)` for the primitive `(q+d)`-grid, the complete closed safe sets are

```text
G_C(r)=W(r)                     for r0<r<=1/q,
G_C(r0)=W(r0) union P_(q+d).
```

The boundary retains isolated primitive-grid points, overlapping endpoints, and all tied owners. For appended harmonic banks `{qm,2qm,...,kqm}`, the theorem gives the exact maximum dichotomy, complete ceiling/maximizer census, corrected terminal mass, affine periods, and exact rational positive/negative certificates within the stated parameter domain.

The scalar first-tail special case and interval-core formulas are inherited and credited. Nonunit promotions and tied-owner data matter: `(q,d,E)` alone does not determine continuation below the chamber. No uniform small-denominator bank, arbitrary-speed theorem, or general below-floor classification follows.

[Harmonic theorem, formulas, and controls](harmonic_owners/README.md)

## Open frontier

The full positive-odd H-spectrum conjecture, global boundedness or unboundedness of exact spectral width, all-order strong-quotient scalar zero-rigidity, and LRC(14) remain open. The completed results sharpen these problems without replacing realizability by a relaxation or replacing proof by finite evidence.
