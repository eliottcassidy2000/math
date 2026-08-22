# Level-three Keller norm recursion: the order-seven boundary and the new prime

**Status: PROVED STRUCTURAL FIXED-MAP DIVISOR THEOREM + VERIFIED-EXACT
GLOBAL IDENTITY + FINITE-EXACT DEGREE-27 SPECIALIZATION; CONDITIONAL
ALL-LEVEL RECURSION.**  The structural proof and global reconstruction below
have an independent exact audit and are canonized as THM-3495.  The all-level
recursion remains conditional.

## Inheritance pass and live concepts

The closest proved mechanism is THM-2582 --
`odd-block-discriminant-tower-and-composite-jelonek-square-class`: an odd
norm-product carries one copy of the outer discriminant class.  The canonical
hostile is the old prediction that the level-two odd divisor was `LH`;
THM-2582 proves that the old `L` cancels and only `H` remains.  The corrected
near miss is to extrapolate that cancellation without measuring the next
norm.  The least-used relevant sidecar is the valuation of `N(H)` along the
old nonproperness divisor.

The live concept board for this pull was:

1. the degree-three norm as a pushforward of a prime divisor;
2. the Newton/Puiseux face at the generic point of `L=0`;
3. finiteness and etaleness away from the Jelonek boundary;
4. odd-block discriminant recursion; and
5. a single squarefree degree-27 specialization as a genericity certificate.

The Newton face proves polynomiality without printing a 66,146-term
numerator.  Finite-divisor geometry proves irreducibility without asking a
computer to factor that numerator.  The exact expansion then fixes the
constant, content, degrees, and transport hash.

Throughout, `F` is the fixed three-dimensional sporadic Keller map of
THM-2473 -- `sporadic-keller-branch-tower-depressed-trisection-anatomy`.
Target coordinates are `(a,b,c)`,

```text
L = 27a^2c^2 - 18abc + 16a + b^3c - b^2,
T = 4 - 3bc,
S = 27ac^2 - 9bc + 8,
E(w) = Lw^3 + Tw - 2c.                                  (1)
```

THM-2576 -- `composite-jelonek-image-divisor-and-two-component-
nonproperness-law` -- defines the primitive irreducible 361-term polynomial
`H`, with multidegree `(14,21,12)`, total degree `25`, and

```text
closure(F(V(L))) = V(H),              gcd(H,L)=1.         (2)
```

Write `N` for the cubic function-field norm induced by `F`.

## 1. The conditional all-level recursion

Let `Delta_r` be the discriminant square class of the degree-`3^r`
`x`-eliminant for `F^r`, in `K^*/K^{*2}` for `K=Q(a,b,c)`.  Suppose through
level `r+1` that the recursively defined inverse algebra is finite and
separable, each level-`r` block is separable, the three conjugate blocks are
pairwise coprime, and their product is the intended eliminant up to a nonzero
element of `K`.  THM-2582 gives

```text
[Delta_(r+1)] = [N(Delta_r)] [Delta_1]^(3^r)
                = [N(Delta_r)] [Delta_1].                 (3)
```

Therefore, conditionally at every licensed level,

```text
[Delta_r] = product_(j=0)^(r-1) [N^j(Delta_1)]
          = [(-1)^r product_(j=0)^(r-1) N^j(L)],          (4)
```

because `[Delta_1]=[-L]` and every norm degree is odd.  THM-2582 proves

```text
N(L)=H/(64L),                                             (5)
```

so `[Delta_2]=[H]`, while level three reduces to the one missing norm

```text
[Delta_3]=[-L N(H)].                                      (6)
```

Equation (4) is still conditional at arbitrary depth.  The rest of this note
discharges every hypothesis and divisor question in (6), only at level three.

## 2. Fixed-map level-three norm-divisor theorem

There is a primitive irreducible polynomial `J in Z[a,b,c]` such that

```text
N(H) = J/(2^35 L^7).                                      (7)
```

The scale is part of the statement: `J=2^35L^7N(H)` is primitive.  A
concurrent reflection that names a rationally rescaled equation for the same
hypersurface must be renormalized before comparing constants or square
classes.

It has

```text
number of terms                 66,146,
multidegree_(a,b,c)             (86,129,76),
total degree                    157,
gcd(J,LH)                       1,
coefficient-ledger SHA-256
  9aca78e67d33351b2f2fb4dbe8ab5bdff06373fdbd8ef9ec73d29b15bffedefe. (8)
```

Geometrically,

```text
closure(F(V(H))) = V(J),                                  (9)
```

and the restriction `V(H) -> V(J)` has generic degree one.  Thus `J` is the
third distinct Jelonek image prime, not merely a squarefree computational
factor.  Sections 3--5 separate the structural proof from the large exact
normalization.

## 3. Why the pole order is exactly seven

Work at the generic DVR of the irreducible prime `(L)`.  Put

```text
D = 18ac - 3b^2c + 2b,
K0 = 9ac-b,
M = 27ac^2-9bc+26,

Y0 = 81abc^2-72ac-15b^2c+16b,
A1 = 27abc^2+54ac-9b^2c+2b,
A2 = 27ab^2c^2+18abc-48a-9b^3c+10b^2.             (10)
```

The exact inverse chart, reduced in `K[w]/(E)`, is

```text
2S q_y = Y0 + 6Lw - 3K0 Lw^2,
8S q_z = Z0 + 6L A1 w - 9L A2 w^2,                 (11)
```

where direct expansion gives

```text
Y0 + 3K0 T                 = 2D,
-6A1 T - 18A2 c            = -24D,
Z0 + 9A2 T                 = -4ML.                  (12)
```

The fast structural companion rederives (11) from the rational inverse
coordinates by cross-multiplication modulo `E`; the displayed identities are
not inferred from numerical samples.

At the point

```text
(a,b,c)=(2/27,1,1),       L=0,
(c,T,S,D)=(1,1,1,1/3),                                  (13)
```

all four generic-DVR unit hypotheses hold simultaneously.  Thus `c,T,S,D`
are units at the generic point of `L`.  The Newton polygon of (1) has one
finite root with `v_L(w)=0` and two divergent roots with
`v_L(w)=-1/2`.

On a divergent root put `u=1/w`.  Equation (1) becomes

```text
L = -T u^2 + 2c u^3.                                    (14)
```

Substitute (14) into (11) and use (12):

```text
2S q_y = 2D - 6(T+K0 c)u + 12cu^2,
8S q_z = -24Du + (12A1c+4MT)u^2 - 8Mcu^3.          (15)
```

Consequently

```text
q_y = D/S + O(u),       q_z = -3(D/S)u + O(u^2),
3wq_z-2q_y = -11D/S + O(u).                         (16)
```

Now grade a monomial `x^i y^j z^k` of `H` by `i-k`.  Exact extraction from
the pinned 361-term polynomial gives

```text
max_H(i-k)=7,
in_(i-k)(H) = -63078912 x^7(3xz-2y)^3.              (17)
```

Equations (16)--(17) show that the leading coefficient is nonzero and each
divergent sheet has

```text
v_L(H(q(w)))=-7/2.                                      (18)
```

The finite root specializes at (13) to `w=2`.  The inverse point and exact
value are

```text
q=(2,5/6,-7/8),
H(q)=3393794313700733412883215882425216567
     /359414999291950792704 != 0.                       (19)
```

Hence the finite sheet has generic valuation zero.  Summing the three sheet
valuations proves, without a global expansion,

```text
v_L(N(H)) = -7/2-7/2+0 = -7.                            (20)
```

This is the mechanism behind the previously observed denominator `L^7`.

## 4. Polynomiality, one-prime support, and generic multiplicity one

Let `U=A^3\V(L)`.  THM-2473 proves that `S_F=V(L)`, that every fibre over
`U` has three points, and that the restriction over `U` is proper.  Since it
is also affine and everywhere etale (`det J_F=-2`),

```text
F_U:F^(-1)(U) -> U
```

is finite etale of degree three.  Therefore the norm of the regular function
`H` belongs to

```text
Q[a,b,c,L^(-1)].                                        (21)
```

Equation (20) now proves directly that

```text
G:=L^7N(H) is in Q[a,b,c] and gcd(G,L)=1.                (22)
```

This polynomiality argument depends only on the DVR calculation, not on
printing `G`.

The divisor `V(H)` is irreducible and is not contained in `V(L)`.  Its dense
open part in `F^(-1)(U)` is therefore irreducible of dimension two.  A finite
map preserves dimension, so its image is a closed irreducible surface in
`U`.  Write `V(P)` for its closure in `A^3`, with `P` an irreducible defining
polynomial.

For a finite map, the zero-divisor support of a norm is the image of the
zero-divisor support of the original function.  Thus on `U` the divisor of
`N(H)` is supported on `V(P)` and on no other prime.  After (22), unique
factorization gives

```text
G = u P^e,          u in Q^*,       1 <= e <= 3,          (23)
```

where `e` is the generic degree of `V(H) -> V(P)`; the upper bound comes from
the degree-three finite map.

The exact slice `(b,c)=(1,2)` from the earlier slice companion is

```text
N(H)=K_(1,2)/(2^21 L^7),                                (24)
```

where `K_(1,2)` is a nonconstant irreducible squarefree polynomial of degree
86.  Specializing (23) cannot turn an `e`-th power with `e>1` into (24): the
specialized `P` is neither zero nor constant because (24) is nonconstant.
Therefore

```text
e=1.                                                     (25)
```

This proves that `G` is irreducible up to a rational unit, proves (9), and
proves generic degree one.  The same slice has `gcd(K_(1,2),LH)=1`, so the
new prime is distinct from both old primes.  The global exact gcd computation
independently confirms this distinctness.

## 5. Exact global normalization without the degree-27 discriminant

The global companion uses FLINT integer multivariate arithmetic.  With
`Y=2S q_y` and `Z=8S q_z`, it forms only

```text
R=8^25 S^25 H(w,Y/(2S),Z/(8S))                           (26)
```

and reduces its `w`-degree 50 modulo the cubic (1).  The exact common
coefficient factor is

```text
gcd(coefficients after a common L^24 reduction)
   =2^40 L^21 S^24.                                      (27)
```

Thus, for integral polynomials `B_0,B_1,B_2`,

```text
H(q(w))=(B_0+B_1w+B_2w^2)/(2^35 S L^3).                 (28)
```

An explicit eight-term resultant, independently checked against the abstract
Sylvester resultant, gives

```text
Res_w(E,B_0+B_1w+B_2w^2)=2^70 L^4 S^3 J.                (29)
```

Because `E` has leading coefficient `L`, the norm of its quadratic argument
is the resultant divided by `L^2`.  Combining the cube of the denominator in
(28) with (29) yields exactly (7):

```text
N(H)=[2^70 L^4 S^3 J]/[L^2 (2^35 S L^3)^3]
    =J/(2^35 L^7).                                       (30)
```

The exact reconstruction proves the data in (8), finds `J` squarefree, and
recovers all three earlier one-parameter numerators coefficient-for-
coefficient.  In particular,

```text
J(a,1,2)=2^14 K_(1,2),
J(a,3,1)=K_(3,1),
J(a,1,3)=K_(1,3).                                       (31)
```

The first equality explains why that slice had denominator `2^21` rather
than the generic `2^35`.  Squarefreeness in the global computation is a
hostile control; irreducibility comes from Section 4, not from a squarefree
factorization routine.

### Normalization bridge to the concurrent residue sidecar

The concurrent reflection
`jc_level3_global_exact_pole_seven_and_prime_image_20260816.md` writes
`J_res:=L^7N(H)`.  In the primitive normalization of this note,

```text
J_res=J/2^35.                                            (31a)
```

Accordingly its square class `[-J_res]` is exactly `[-2J]`, not `[-J]` for
the primitive polynomial.  Its boundary-residue and explicit image-point
controls are compatible supplementary checks; equation (31a) must be applied
before importing any constant or square-class statement from that file.

## 6. The actual degree-27 class

A separate exact companion specializes the target to `(1,1,1)`, builds the
successive cubic inverse algebras of dimensions three and nine, and computes
the next norm-product polynomial by 28 exact `9 by 9` multiplication
determinants.  An unused 29th determinant checks the degree-27 interpolation
off-grid.  The primitive result has

```text
degree                         27,
gcd(P_3,P_3')                  1,
content                        1,
rational interpolation denominator
  32687327622020737269760000000000000,
ordered-coefficient SHA-256
  fa8ba9f1cb850116c347f6e31100d1902dea3ee1c11c2b6548f7280aa9f01d50. (32)
```

Since one specialization of the product of the three degree-nine blocks is
squarefree of full degree 27, generic separability of every block and generic
pairwise block coprimality follow.  This licenses THM-2582 for the actual
level-three eliminant.  Equations (6)--(7) now give

```text
[Delta_3]=[-L N(H)]
         =[-J/(2^35 L^6)]
         =[-2J]              in K^*/K^{*2}.              (33)
```

Thus the old `L` exponent is `1-7=-6`, `H` has valuation zero, and the sole
odd polynomial prime is the new irreducible `J`.  The constant class is
exactly `-2`.  No degree-27 discriminant was expanded.

## 7. The third nonproperness set

THM-2576's set-level composition law gives

```text
S_(F^3)=V(L) union F(V(L)) union F^2(V(L)).               (34)
```

The left side is closed.  Moreover, density of `F(V(L))` in `V(H)` and
continuity give
`closure(F(F(V(L))))=closure(F(V(H)))`.  Taking closures in the finite union
and using (2) and (9) therefore gives

```text
S_(F^3)=V(L) union V(H) union V(J)=V(LHJ).                (35)
```

Since `L,H,J` are pairwise distinct irreducibles, this set has exactly three
irreducible components.  This is a reduced set statement; it does not assign
scheme multiplicities to the nonproperness set.

## Scope and stopping boundary

- **PROVED, fixed map:** the order-seven pole, polynomiality, one-prime norm
  support, irreducibility, generic degree one, (9), (33), and (35).
- **VERIFIED-EXACT, independently replayed:** the global normalization (7),
  the coefficient data (8), the gcds and squarefreeness, and all three slice
  recoveries (31).
- **FINITE-EXACT with a proved generic consequence:** the one target in (32)
  and the resulting generic degree-27 separability/coprimality.
- **CONDITIONAL:** the all-level recursion (3)--(4) beyond the now-licensed
  third level.
- **OPEN:** whether later image primes remain distinct or generically
  degree-one, and whether the newest-factor cancellation persists at every
  level.
- **Out of scope:** this fixed map is a three-dimensional counterexample
  tower.  Nothing here bears on `JC(2)`, proves a statement for arbitrary
  Keller maps, or supplies a positive Jacobian-conjecture implication.

Reproduce with

```bash
python3 04-computation/keller_level_three_norm_divisor_structure_20260816.py
python3 -O 04-computation/keller_level_three_norm_divisor_structure_20260816.py
python3 04-computation/keller_level_three_global_norm_probe_20260816.py
python3 -O 04-computation/keller_level_three_global_norm_probe_20260816.py
python3 04-computation/keller_level_three_squarefree_tower_probe_20260816.py
python3 -O 04-computation/keller_level_three_squarefree_tower_probe_20260816.py
```

Ordinary, optimized, and stored outputs are byte-identical for all three
companions.  Script/output hashes are recorded with the matching artifacts in
`05-knowledge/results/INDEX.md`.
