---
id: THM-3530
title: "Fixed Keller all-level image-prime and exact component tower"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE UNDER INDEPENDENT HOSTILE AUDIT.
  Conditional on the reserved THM-3529 finite-sheet unit proof, the candidate
  uses finite-etale norm divisors and primitive complete top-face exponents
  to force generic image degree one at every raw rung.  No statement in this
  file is proved canon until both the dependency and this induction are
  independently promoted.
source: codex/packet-power-image-induction/2026-08-16
depends_on: []
related:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2576-composite-jelonek-image-divisor-and-two-component-nonproperness-law
  - THM-3504-level-four-sporadic-keller-image-prime-and-four-component-nonproperness
  - THM-3528-fixed-keller-all-level-cleared-norm-polynomiality-and-finite-sheet-defect
  - THM-3529-fixed-keller-complete-packet-finite-sheet-unit
script: 04-computation/keller_packet_power_image_prime_induction_audit_20260816.py
output: 05-knowledge/results/keller_packet_power_image_prime_induction_audit_20260816.out
script_sha256: b756f9707a5bce1885069a75e48f27c7a0163b321f9d606ce33c9e0e171fe67c
output_sha256: cae9264170d3fca662d7a2dcbc4380a26fff98be279e8f57fd600a2389150dc4
semantic_sha256: ab5aa0fe7f28607c766ea4e3d1d42a7b84b4e9843ce2e35df49ebefd06bace75
hash_basis: LF-normalized bytes
---

# THM-3530 -- primitive packet faces force an all-level prime image tower

**RESERVED / PROVISIONAL PROOF CANDIDATE UNDER INDEPENDENT HOSTILE
AUDIT.**  This file is not in the proved dependency graph.  Its proof uses
the still-reserved THM-3529 candidate and must not be cited as canon before
both audits are accepted.

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

## 1. Candidate theorem

For every `n>=0`, the candidate conclusions are:

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

has exactly `r` irreducible components.  This would promote the fixed-map
component-count prediction of HYP-9033 from depths one through four to every
depth.

## 2. The finite-etale norm of one prime is one prime power

Let `P` be an absolutely irreducible complete packet `A(e,m)`.  Put

```text
D=V(L),                 U=A^3\D,
X_U=F^-1(U).                                               (6)
```

THM-2473 proves that

```text
F_U:X_U->U                                                (7)
```

is finite etale of degree three.  Reserved THM-3529 would prove that the
finite inverse divisor `F^-1(D)` is not a component of `V(P)`.  Hence

```text
E=V(P) intersect X_U                                     (8)
```

is a nonempty dense irreducible divisor in `V(P)`.  Because `(7)` is finite,
its image is a closed irreducible surface in `U`.  Its closure is
Galois-stable because `F` and `P` are defined over `Q`, and it is geometrically
irreducible because it is the closure of the image of a geometrically
irreducible variety.  Write it as `V(R)`, with `R in Q[a,b,c]` absolutely
irreducible and not associated to `L`.

At the generic point of `V(R)`, the reduced prime divisor `E` has

```text
d=deg(E->V(R)) in {1,2,3}.                              (9)
```

Etaleness makes the order of `P` along `E` equal to one.  The divisor of the
finite-algebra norm on `U` is therefore the pushforward divisor

```text
div_U N(P)=d (V(R) intersect U).                        (10)
```

Equivalently, in the localization `A[L^-1]`,

```text
N(P)=u R^d,                                             (11)
```

where `u` is a unit.  Since `A` is factorial and `L` is irreducible, every
unit of `A[L^-1]` is `cL^k` for `c!=0` and `k in Z`.  Multiplying `(11)` by
the packet clearing power gives

```text
Q(P)=L^eN(P)=c L^(e+k) R^d.                             (12)
```

THM-3528 identifies `ord_L Q(P)` with the finite-sheet defect; THM-3529's
candidate value zero would make `e+k=0`.  Thus

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
gcd(e_n,m_n)=1 for every n>=1.                          (19)
```

Indeed the exact Cassini identity

```text
e_n m_(n+1)-m_n e_(n+1)=3(-8)^n                        (20)
```

restricts any common prime to `2` or `3`, while

```text
(e_n,m_n)=(1,3) mod 6,                 n>=1,           (21)
```

excludes both.  When `Q(P_n)=P_(n+1)`, equations `(18)`--`(21)` leave only

```text
d=1.                                                    (22)
```

Thus the output is a scalar multiple of the absolutely irreducible image
equation `R`, and the restriction has generic degree one.  This proves the
inductive step `(2)`--`(4)`.

The base is `P_0=L`, whose absolute irreducibility is part of the fixed
Jelonek theorem.  Induction now gives `(2)`--`(4)` at every level.

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

One way to keep `(25)` honest is to use constructibility: the dense image at
each stage contains a dense open of the irreducible target, and the next
generic-degree-one map remains dominant on that open.  No equality between
an unclosed image and its closure is silently assumed.

The left side of `(24)` is closed.  Taking the closure of its finite union
and using `(25)` gives

```text
S_(F^r)=union_(j=0)^(r-1) V(P_j)
       =V(product_(j=0)^(r-1) P_j).                    (26)
```

Pairwise distinct absolute primes and `(23)` prove that `(26)` has exactly
`r` irreducible components, establishing `(5)`.

## 5. What closes, and what does not

If promoted, the theorem would close four fixed-map questions at once:

1. every raw cleared norm is an absolute prime;
2. every newest prime is the closure of the preceding prime's image;
3. every newest image restriction has generic degree one; and
4. the `r`th self-iterate has exactly `r` Jelonek components.

It would also show that THM-3528's branch-transplant defect word is identically
zero on this raw orbit.  The transplant theorem remains a valid conditional
monoid statement for other packet inputs, but no raw tower return occurs.

The result still would **not** prove:

- generic separability or full degree of every later coordinate eliminant;
- an all-level discriminant square-class recursion;
- exact positive discriminant multiplicities or singularity types;
- primitive integral normalizations, term counts, or global multidegrees;
- a classification of all maps in a fibre-degree grade;
- the same assertion for a tame conjugate without transporting its packet
  weights and finite divisor; or
- any statement about arbitrary Keller maps, `JC(2)`, `DC(2)`, LRC, or a
  general Jacobian-conjecture classification.

## 6. Equality boundary and hostile controls

The candidate has sharp failure modes.

1. **Primitive grade is essential.**  Since `L` has `A(1,0)`, the complete
   packets `L^2` and `L^3` have grades `(2,0)` and `(3,0)` and literal top
   faces `(16x)^2` and `(16x)^3`.  Packet completeness alone does not exclude
   image multiplicity.
2. **Input primeness is essential.**  For reducible `P`, different source
   components can have different image primes, so the norm need not be one
   prime power.
3. **The finite unit is essential.**  Without THM-3529, equation `(12)` may
   retain an old-`L` factor and the one-prime-power reduction fails.
4. **Etaleness is essential.**  Ramification could multiply the divisor
   order independently of the generic image degree.
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
