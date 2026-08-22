---
id: THM-3530
title: "Fixed Keller all-level image-prime and exact component tower"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the fixed
  sporadic Keller map, every raw cleared-norm rung is absolutely irreducible,
  is the reduced equation of the preceding rung's image closure, and maps to
  the next rung with generic degree one.  The nth self-iterate has exactly n
  reduced Jelonek components.  The mechanism is the finite-etale prime-power
  norm law plus coprime complete top-face exponents.
source: codex/packet-power-image-induction/2026-08-16
audit: >
  The independent audit verified the localized prime contraction, geometric
  integrality of the image, norm-valuation exponent, localized-unit removal,
  UFD power obstruction, all-level grade primitivity, image-closure induction,
  and reduced component indexing.  It supplied hostiles for ordinary versus
  absolute irreducibility, generic-degree norm powers, localized boundary
  units, and nonprimitive top grades.
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2570-jelonek-cusp-cylinder-normalization-and-conductor
  - THM-2576-composite-jelonek-image-divisor-and-two-component-nonproperness-law
  - THM-3528-fixed-keller-all-level-cleared-norm-polynomiality-and-finite-sheet-defect
  - THM-3529-fixed-keller-complete-packet-finite-sheet-unit
related:
  - THM-3531-fixed-keller-intrinsic-all-level-discriminant-square-class
  - THM-3504-level-four-sporadic-keller-image-prime-and-four-component-nonproperness
  - HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid
script: 04-computation/keller_packet_power_image_prime_induction_audit_20260816.py
output: 05-knowledge/results/keller_packet_power_image_prime_induction_audit_20260816.out
script_sha256: a9dca29c9a82ec97082adc87afc1b0695abbe3fcc3a5c5df31dabcf227f40b7a
output_sha256: 66f3c578d9a021df0ed07c8b260015ac7b6ab284dff320af24f34fe8d5794656
semantic_sha256: ab5aa0fe7f28607c766ea4e3d1d42a7b84b4e9843ce2e35df49ebefd06bace75
hash_basis: LF-normalized bytes
---

# THM-3530 -- primitive packet faces force an all-level prime image tower

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Retain the fixed sporadic map `F:C^3->C^3`, its irreducible boundary `L`, and
the raw cleared-norm tower

```text
P_0=L,
P_(n+1)=L^e_n N(P_n),
(e_0,m_0)=(1,0),
(e_(n+1),m_(n+1))=(7e_n-2m_n,3e_n-2m_n).              (1)
```

THM-3528 proves that every `P_n` is a nonzero polynomial with complete packet
`A(e_n,m_n)`.  Rational scalar normalizations do not affect any statement
below.

## 1. The theorem

For every `n>=0`:

```text
P_n is absolutely irreducible;                          (2)

closure(F(V(P_n)))=V(P_(n+1));                         (3)

V(P_n)->V(P_(n+1)) has generic degree one.             (4)
```

The polynomials `P_n` are pairwise nonassociate.  Consequently, for every
`r>=1`,

```text
S_(F^r)=V(P_0 P_1 ... P_(r-1))                         (5)
```

has exactly `r` irreducible components as a reduced algebraic set.  This
promotes the fixed-map
component-count prediction of HYP-9033 from depths one through four to every
depth.

## 2. The finite-etale norm of one prime is one prime power

Let `P` be an absolutely irreducible complete packet `A(e,m)`.  Put

```text
A_0=Q[a,b,c],
D=V(L),                 U=Spec(A_0[L^-1]),
X_U=F^-1(U).                                               (6)
```

THM-2473 proves that

```text
F_U:X_U->U                                                (7)
```

is finite etale of degree three.  THM-3529 proves that the
finite inverse divisor `F^-1(D)` is not a component of `V(P)`.  Hence

```text
E=V(P) intersect X_U                                     (8)
```

is a nonempty dense geometrically integral open subset of `V(P)` and a
closed prime divisor in `X_U`.  Because `(7)` is finite,
its image is a closed irreducible surface in `U`.  Write this equality and
its global closure as

```text
F_U(E)=V(R) intersect U,             closure(F_U(E))=V(R). (8a)
```

Choose `R` primitive with respect to `L`.  The closure is
Galois-stable because `F` and `P` are defined over `Q`, and it is geometrically
irreducible because it is the closure of the image of a geometrically
irreducible variety.  Write it as `V(R)`, with `R in Q[a,b,c]` absolutely
irreducible and not associated to `L`.

At the generic point of `V(R)`, the reduced prime divisor `E` has

```text
d=deg(E->V(R)) in {1,2,3}.                              (9)
```

The prime polynomial `P` has order one along `E`, and etaleness contributes
no ramification multiplier.  The divisor of the finite-algebra norm on `U`
is therefore the pushforward divisor

```text
div_U N(P)=d (V(R) intersect U).                        (10)
```

Equivalently, in the localization `A_0[L^-1]`,

```text
N(P)=u R^d,                                             (11)
```

where `u` is a unit.  Since `A_0` is factorial and `L` is irreducible, every
unit of `A_0[L^-1]` is `cL^k` for `c!=0` and `k in Z`.  Multiplying `(11)` by
the packet clearing power gives

```text
Q(P)=L^eN(P)=c L^(e+k) R^d.                             (12)
```

THM-3528 identifies `ord_L Q(P)` with the finite-sheet defect; THM-3529's
value zero makes `e+k=0`.  Thus

```text
Q(P)=cR^d,                 1<=d<=3.                    (13)
```

This is the exact geometric reduction.  It does not assume that a norm of an
irreducible polynomial is irreducible: the exponent `d` is precisely the
possible failure.

## 3. A primitive complete face forces `d=1`

By THM-3528, `Q(P)` has complete output packet

```text
A(e',m'),                (e',m')=(7e-2m,3e-2m).        (14)
```

Its complete maximum-`lambda` face is

```text
in_max-lambda Q(P)=c_0 x^e'(3xz-2y)^m',       c_0!=0. (15)
```

Initial forms multiply in the monomial-weight graded domain.  Equation
`(13)` therefore implies

```text
in_max-lambda Q(P)=c (in_max-lambda R)^d.              (16)
```

The two polynomials

```text
x,                       h=3xz-2y                     (17)
```

are coprime irreducibles in `C[x,y,z]`.  Unique factorization in `(15)`--
`(16)` forces

```text
d divides e',                 d divides m'.             (18)
```

For the actual raw orbit, the packet arithmetic gives

```text
gcd(e_n,m_n)=1 for every n>=0.                          (19)
```

Indeed, write `v_n=(e_n,m_n)^T=M^n(1,0)^T`.  If an odd prime `p`
divided both coordinates of `v_n`, then `det M=-8` would make `M` invertible
modulo `p`, forcing `(1,0)^T=M^(-n)v_n=0 mod p`, a contradiction.  Thus the
gcd is a power of two.  But

```text
e_(n+1)=7e_n-2m_n=e_n mod 2,             e_0=1,        (20)
```

so every `e_n` is odd and `(19)` follows.  The independent Cassini check

```text
e_n m_(n+1)-m_n e_(n+1)=3(-8)^n                         (21)
```

gives a second route together with `(e_n,m_n)=(1,3) mod 6` for `n>=1`.
When `Q(P_n)=P_(n+1)`, equations `(18)`--`(21)` leave only

```text
d=1.                                                    (22)
```

Thus the output is a scalar multiple of the absolutely irreducible image
equation `R`, and the restriction has generic degree one.  Since `E` is
dense in `V(P)`, continuity gives `F(V(P)) subseteq closure(F(E))`; the
reverse closure inclusion is immediate from `E subseteq V(P)`.  Hence
`closure(F(V(P)))=V(R)`, proving the inductive step `(2)`--`(4)`.

The base is `P_0=L`, whose absolute irreducibility follows from THM-2570's
geometrically integral cusp-cylinder normalization.  Induction now gives
`(2)`--`(4)` at every level.

## 4. Distinctness and the exact Jelonek component count

The packet cone has `0<=2m_n<=e_n`.  Therefore

```text
e_(n+1)-e_n=6e_n-2m_n >=5e_n>0.                        (23)
```

Associated polynomials have the same maximum-`lambda` weight.  Strict growth
in `(23)` makes all `P_n` pairwise nonassociate.  Since each is absolutely
irreducible, the hypersurfaces `V(P_n)` are pairwise distinct prime
divisors.

THM-2576 proves the exact set law

```text
S_(F^r)=union_(j=0)^(r-1) F^j(S_F),          S_F=V(P_0). (24)
```

Iterating `(3)` shows

```text
closure(F^j(V(P_0)))=V(P_j).                            (25)
```

Equation `(25)` follows by applying the continuity/dense-open argument in
Section 3 at each stage.  No equality between an unclosed image and its
closure is silently assumed.

The left side of `(24)` is closed.  Taking the closure of its finite union
and using `(25)` gives

```text
S_(F^r)=union_(j=0)^(r-1) V(P_j)
       =V(product_(j=0)^(r-1) P_j).                    (26)
```

Pairwise distinct absolute primes and `(23)` prove that `(26)` has exactly
`r` irreducible components, establishing `(5)`.

## 5. What closes, and what does not

The theorem closes four fixed-map questions at once:

1. every raw cleared norm is an absolute prime;
2. every newest prime is the closure of the preceding prime's image;
3. every newest image restriction has generic degree one; and
4. the `r`th self-iterate has exactly `r` Jelonek components.

It also shows that THM-3528's branch-transplant defect word is identically
zero on this raw orbit.  The transplant theorem remains a valid conditional
monoid statement for other packet inputs, but no raw tower return occurs.

The result does **not** prove:

- generic separability or full degree of every later coordinate eliminant;
- an all-level discriminant square-class recursion (subsequently proved
  intrinsically in THM-3531);
- exact positive discriminant multiplicities or singularity types;
- primitive integral normalizations, term counts, or global multidegrees;
- a classification of all maps in a fibre-degree grade;
- the same assertion for a tame conjugate without transporting its packet
  weights and finite divisor; or
- any statement about arbitrary Keller maps, `JC(2)`, `DC(2)`, LRC, or a
  general Jacobian-conjecture classification.

## 6. Equality boundary and hostile controls

The theorem has sharp failure modes.

1. **Primitive grade is essential.**  Since `L` has `A(1,0)`, the complete
   packets `L^2` and `L^3` have grades `(2,0)` and `(3,0)` and literal top
   faces `(16x)^2` and `(16x)^3`.  Packet completeness alone does not exclude
   image multiplicity.
2. **Input primeness is essential.**  For reducible `P`, different source
   components can have different image primes, so the norm need not be one
   prime power.
   Ordinary rational irreducibility is also insufficient:
   `x^2+y^2` is irreducible over `Q` but splits over `C`.
3. **The finite unit is essential.**  Without THM-3529, equation `(12)` may
   retain an old-`L` factor and the one-prime-power reduction fails.
4. **Etaleness is essential.**  Ramification could multiply the divisor
   order independently of the generic image degree.
   Even without ramification, the finite map `t=x^2` sends the prime `y=0`
   to `y=0` with generic degree two and norm `N(y)=y^2`; this is the sharp
   witness that the exponent in `(13)` really is the restriction degree.
5. **The complete face is essential.**  A selected leading monomial cannot
   support the UFD inference `(18)` if unrecorded equal-weight terms exist.

## Reproduction

```text
python -B 04-computation/keller_packet_power_image_prime_induction_audit_20260816.py
python -B -O 04-computation/keller_packet_power_image_prime_induction_audit_20260816.py
```

Normal and optimized transcripts match the stored output.  The companion
checks the packet orbit through depth 24, both primitivity routes, strict
grade growth, the two top-face primes, every square/cube divisibility gate,
the component indexing through depth 12, and the sharp `L^2/L^3` hostiles.
