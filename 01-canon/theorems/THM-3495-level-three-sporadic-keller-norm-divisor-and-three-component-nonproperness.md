---
id: THM-3495
title: "Level-three sporadic Keller norm divisor and three-component nonproperness"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF-AUDITED.  For the fixed
  sporadic Keller map of THM-2473, the third norm divisor is a new primitive
  absolutely irreducible polynomial J with N(H)=J/(2^35 L^7).  Consequently
  the third-iterate nonproperness set is V(LHJ), with exactly three
  irreducible components, and the degree-27 eliminant has square class
  [-2J].  This is a fixed-map theorem, not a Jacobian-conjecture result or an
  all-level newest-factor law.
source: codex/frontier-many/2026-08-16
audit: >
  The generic-DVR proof, norm support and image-multiplicity argument,
  resultant normalization, nonproperness composition law, and degree-27
  squarefree specialization were independently reconstructed and audited.
  A disjoint audit then computed the previously unused (b,c)=(1,1) norm
  slice, independently certified the normalization residual by a Rabin test,
  matched the level-three leading coefficient to N^2(L), and found a direct
  split 3-to-9 finite-field fibre whose nine cubic product is squarefree of
  degree 27.
  All five pinned companions replay byte-identically under ordinary and
  optimized Python; the global reconstruction independently reproduced the
  coefficient gcd, resultant content and degrees, J term count and degrees,
  old-factor gcds, and squarefreeness.
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2576-composite-jelonek-image-divisor-and-two-component-nonproperness-law
  - THM-2582-odd-block-discriminant-tower-and-composite-jelonek-square-class
related:
  - HYP-9033
  - MISTAKE-290
scripts:
  - 04-computation/keller_level_three_norm_divisor_structure_20260816.py
  - 04-computation/keller_level_three_global_norm_probe_20260816.py
  - 04-computation/keller_level_three_squarefree_tower_probe_20260816.py
  - 04-computation/jc_level3_global_norm_independent_audit_20260816.py
  - 04-computation/jc_level3_degree27_split_finite_field_audit_20260816.py
outputs:
  - 05-knowledge/results/keller_level_three_norm_divisor_structure_20260816.out
  - 05-knowledge/results/keller_level_three_global_norm_probe_20260816.out
  - 05-knowledge/results/keller_level_three_squarefree_tower_probe_20260816.out
  - 05-knowledge/results/jc_level3_global_norm_independent_audit_20260816.out
  - 05-knowledge/results/jc_level3_degree27_split_finite_field_audit_20260816.out
script_sha256:
  - f4255f6a6918458fb877523329061a66dfdef3a80b7d54cff74394b88c2f6628
  - 37bed904530acbddc17f3f612fcc0d4e8da85b8caddc1b3c86d49b80550e1559
  - fbb2d20388099377eb2498b1fb102f4e2a32785afa02a0ca00792f27c3d6bd3e
  - cb429497fbfbdf4cb538967bd472ab10051bd6fac33be10d14081b09e1543215
  - 9f8b6b4aad75b17eed8a6fd618cb661f3e61b85361b44ba676462148a2c5afc5
output_sha256:
  - 473cffeb80859a3bf91ab3c77fc6144089d2c6025735a81a7d782a8e774dbcf8
  - 1736ac542fd20d6782bc494b1199099e6d08f44830efd5d9c041764925d3db97
  - 4f387efdef50fe7611d51f394ab4dd416274955fe30dbe22368a8a66ba452b10
  - 62fef84710ef14c3a021b196bb4ba3589ad672a3c275a2cfb3c24a25836aea54
  - af2be3f27d524937836b93aacd5fa0700d5aa930c798dcf142b10206619fc520
coefficient_ledger_sha256: 9aca78e67d33351b2f2fb4dbe8ab5bdff06373fdbd8ef9ec73d29b15bffedefe
degree_27_ledger_sha256: fa8ba9f1cb850116c347f6e31100d1902dea3ee1c11c2b6548f7280aa9f01d50
hash_basis: raw LF bytes for files; ordered exact coefficient ledgers as stated
---

# THM-3495 -- the third norm divisor of the fixed sporadic Keller map

**PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF-AUDITED.**

Let `F:C^3->C^3` be the fixed sporadic Keller map of THM-2473, let

```text
L=27a^2c^2-18abc+16a+b^3c-b^2,
T=4-3bc,
S=27ac^2-9bc+8,
E(w)=Lw^3+Tw-2c,
```

and let `H` be THM-2576's primitive irreducible equation of
`closure(F(V(L)))`.  Write `N` for the cubic function-field norm induced by
`F`.

## 1. Theorem

There is a primitive absolutely irreducible `J in Z[a,b,c]` such that

```text
N(H)=J/(2^35 L^7).                                      (1)
```

Equivalently, the normalization in this theorem is exactly
`J=2^35 L^7 N(H)`.  Any reflection using a rationally rescaled equation for
the same hypersurface must be converted to this primitive normalization
before its constants or discriminant square class are compared with (1) and
(5).

The polynomial `J` has 66,146 terms, multidegree `(86,129,76)`, total degree
157, and

```text
gcd(J,LH)=1,
SHA256(lex coefficient ledger of J)
 =9aca78e67d33351b2f2fb4dbe8ab5bdff06373fdbd8ef9ec73d29b15bffedefe. (2)
```

Moreover,

```text
closure(F(V(H)))=V(J),                                  (3)
```

the dominant restriction `V(H)->V(J)` has generic degree one,

```text
S_(F^3)=V(LHJ),                                         (4)
```

with exactly three irreducible components, and the actual degree-27
`x`-eliminant of `F^3` is generically squarefree with discriminant square
class

```text
[Delta_3]=[-2J] in Q(a,b,c)^*/Q(a,b,c)^{*2}.             (5)
```

## 2. Generic boundary valuation

In `Q(a,b,c)[w]/(E)`, exact reduction of the inverse chart gives

```text
2S q_y=Y0+6Lw-3(9ac-b)Lw^2,
8S q_z=Z0+6LA1w-9LA2w^2,                                (6)
```

where, with

```text
D=18ac-3b^2c+2b,       M=27ac^2-9bc+26,
K0=9ac-b,
Y0=81abc^2-72ac-15b^2c+16b,
A1=27abc^2+54ac-9b^2c+2b,
A2=27ab^2c^2+18abc-48a-9b^3c+10b^2,
Z0=-2916a^3c^4+2916a^2bc^3-4536a^2c^2
   +621ab^3c^3-1026ab^2c^2+504abc+64a
   -207b^4c^2+454b^3c-256b^2,
```

one has the polynomial identities

```text
Y0+3K0T=2D,
-6A1T-18A2c=-24D,
Z0+9A2T=-4ML.                                           (7)
```

At the generic DVR of `(L)`, `c,T,S,D` are units: the simultaneous witness
`(a,b,c)=(2/27,1,1)` on `L=0` gives `(c,T,S,D)=(1,1,1,1/3)`.
The Newton polygon of `E` has one root of valuation zero and two roots of
valuation `-1/2`.  On either divergent root put `u=1/w`; then
`L=-Tu^2+2cu^3`, and (6)--(7) give

```text
q_y=D/S+O(u),        q_z=-3(D/S)u+O(u^2),
3wq_z-2q_y=-11D/S+O(u).                                 (8)
```

Exact extraction from the pinned polynomial `H` gives the Newton face

```text
max(i-k)=7,
in(H)=-63078912 x^7(3xz-2y)^3.                           (9)
```

Thus each divergent sheet contributes valuation `-7/2`.  The finite sheet
is generically a unit: at the same target its finite root is `w=2`, its
inverse point is `(2,5/6,-7/8)`, and

```text
H(2,5/6,-7/8)
=3393794313700733412883215882425216567
 /359414999291950792704 !=0.                            (10)
```

Consequently

```text
v_L(N(H))=-7.                                            (11)
```

## 3. Polynomiality and the unique image prime

Put `U=A^3\V(L)`.  THM-2473 proves that `F` has three points in every fibre
over `U`, is proper there, and is everywhere etale.  Hence
`F^{-1}(U)->U` is finite etale of degree three, so `N(H)` belongs to
`Q[a,b,c,L^{-1}]`.  Equation (11) implies

```text
G:=L^7N(H) in Q[a,b,c],       gcd(G,L)=1.                (12)
```

The exact point `p=(3,-1,0)` satisfies `H(p)=0`, while
`F(p)=(10,-46,33)` and `L(F(p))=-504`.  Hence `V(H)` is not contained in
`F^{-1}(V(L))`, so the open `V(H) intersect F^{-1}(U)` is dense,
irreducible, and two-dimensional.
Its finite image is one closed irreducible surface in `U`; write `V(P)` for
its closure, with `P` an irreducible defining polynomial.  The zero-divisor
support of a finite norm is precisely the finite image of the zero-divisor
support of its argument.  Hence unique factorization and (12) give

```text
G=uP^e,          u in Q^*,          1<=e<=3,             (13)
```

where `e` is the generic degree of the restriction `V(H)->V(P)`.

The exact slice `(b,c)=(1,2)` gives

```text
N(H)=K/(2^21L^7),                                       (14)
```

with `K` irreducible and squarefree of degree 86.  Specializing (13) shows
that `e>1` is impossible; the specialized `P` cannot be zero or constant
because (14) is nonconstant.  Thus `e=1`, `G` is irreducible up to a rational
unit, and `V(P)=V(G)`.  The same slice has `gcd(K,LH)=1`, proving distinctness
from both old primes.  This proves (3) and generic degree one before any
global numerator is printed.  Since the argument concerns the geometric
image of the irreducible complex hypersurface, `J` is absolutely
irreducible.

## 4. Exact normalization

Set `Y=2S q_y`, `Z=8S q_z` and reduce

```text
8^25S^25H(w,Y/(2S),Z/(8S))
```

modulo `E`.  Exact integer multivariate arithmetic finds common coefficient
factor `2^40L^21S^24` after putting the reduction over `L^24`, hence

```text
H(q(w))=(B0+B1w+B2w^2)/(2^35SL^3).                      (15)
```

The exact eight-term resultant is

```text
Res_w(E,B0+B1w+B2w^2)=2^70L^4S^3J,                      (16)
```

where `J` is primitive and has the data (2).  Dividing the resultant by the
`L^2` nonmonic norm factor and by the cube of the denominator in (15) gives
(1).  Exact gcd and squarefree checks pass.  Cross-specialization recovers

```text
J(a,1,2)=2^14K_(1,2),
J(a,3,1)=K_(3,1),
J(a,1,3)=K_(1,3)                                        (17)
```

coefficient-for-coefficient.  The irreducibility assertion comes from
Section 3, not from the squarefree computation.

## 5. Nonproperness and discriminant consequences

THM-2576's composition law gives

```text
S_(F^3)=V(L) union F(V(L)) union F^2(V(L)).              (18)
```

This set is closed.  Density of `F(V(L))` in `V(H)` and continuity give
`closure(F(F(V(L))))=closure(F(V(H)))`.  Taking closures and using (3)
yields (4).  Pairwise coprimality makes the three displayed hypersurfaces the
three irreducible components.

At target `(1,1,1)`, an exact triangular-algebra calculation produces a
primitive squarefree norm-product polynomial of full degree 27.  It uses 28
exact determinant samples plus an off-grid 29th determinant; its ordered
coefficient hash is

```text
fa8ba9f1cb850116c347f6e31100d1902dea3ee1c11c2b6548f7280aa9f01d50. (19)
```

One full-degree squarefree specialization proves generic block separability
and pairwise coprimality.  THM-2582 is therefore applicable to the actual
degree-27 eliminant and gives

```text
[Delta_3]=[-LN(H)]=[-J/(2^35L^6)]=[-2J],                 (20)
```

which is (5).  The old `L` has even exponent `-6`, and `H` is absent.

An independent good-reduction route reaches the same genericity conclusion
without the triangular algebra.  Over `F_101`, the target `(93,28,83)` has
three explicit first preimages and three explicit second preimages above each
of them.  All nine third cubic cores retain degree three; grouped by the
first preimage they give three squarefree pairwise-coprime degree-nine
blocks, and their direct product is squarefree of degree 27.  Since every
inverse-chart denominator used by this fibre is nonzero modulo 101, this is
a lawful good-reduction witness that the characteristic-zero generic
discriminant is not identically zero.

## 6. Exact companions and hostile controls

- `keller_level_three_norm_divisor_structure_20260816.py` cross-reduces the
  inverse chart, checks every identity in (6)--(9), and checks the finite
  point (10).  It is independent of the global `J` expansion.
- `keller_level_three_global_norm_probe_20260816.py` pins the raw `H` pickle
  and its coefficient ledger, asserts (15)--(16), freezes the `J` ledger,
  checks `gcd(J,LH)=1` and squarefreeness, and recovers all three old slices.
- `keller_level_three_squarefree_tower_probe_20260816.py` checks both inverse
  graphs, the degree-nine squarefree control, the full degree-27 squarefree
  polynomial, and an off-grid determinant.
- `jc_level3_global_norm_independent_audit_20260816.py` computes the fresh
  `(b,c)=(1,1)` norm by a closed regular-representation determinant.  Its
  reduced denominator is exactly `2^35L^7`; the primitive degree-86 slice is
  irreducible and coprime to `LH`.  It also proves the 527-term boundary
  residual irreducible by the Rabin criterion for `P(-1,lambda) mod 449`,
  checks both rational transverse residue values, and verifies
  `LC(P_3)=N^2(L)=J/(2^47L^6H)` at `(1,1,1)`.
- `jc_level3_degree27_split_finite_field_audit_20260816.py` supplies the
  independent `F_101` split-fibre witness above and checks the three block
  gcds plus the full degree-27 derivative gcd directly.

All five companions replay identically under ordinary and optimized Python.
The global reconstruction was independently reimplemented and reproduced the
coefficient gcd, resultant content and degrees, `J` degrees and term count,
both old-factor gcds, and squarefreeness.

## 7. Scope and next boundary

This theorem concerns one fixed polynomial map in dimension three and its
first three self-compositions.  It does not apply to arbitrary Keller maps,
does not bear on `JC(2)`, and does not prove the all-level newest-factor law.
At depth four one must again determine a new norm divisor, its boundary
valuation, image multiplicity, distinctness, and generic block separability.
