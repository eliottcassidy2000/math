---
id: THM-3528
title: "Fixed Keller all-level cleared-norm polynomiality and finite-sheet defect identity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the fixed
  sporadic Keller inverse chart, a nonzero polynomial with a complete packet
  A(e,m) has L^e times its cubic norm in Q[a,b,c], and the result has packet
  A(7e-2m,3e-2m).  Its exact old-L multiplicity is the nonnegative order on
  the regular finite inverse sheet.  Consequently the raw cleared-norm tower
  from L is polynomial and carries the Pell-57 packet at every level.  This
  does not prove later L-coprimality, newest-image status, irreducibility,
  separability, distinct nonproperness components, an arbitrary-map law, or
  any general Jacobian-conjecture claim.
source: codex/all-level-cleared-norm/2026-08-16
audit: >
  A clean-room reciprocal-cubic audit independently verifies finite-etale
  norm regularity, the one-plus-ramified-pair Newton decomposition, complete
  max-lambda noncancellation, ramification weights, the nonmonic resultant
  hostile, rational scalar normalization, and 15,251 admissible packets.
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-3498-level-four-old-boundary-cancellation-and-degree81-discriminant-gate
  - THM-3506-fixed-keller-five-face-norm-transform-and-271-99-boundary
  - THM-3522-fixed-keller-five-face-renewal-propagation
related:
  - MISTAKE-415
  - THM-3521-fixed-R5-finite-sheet-unit-and-next-old-L-clearing
  - THM-3523-fixed-R6-finite-sheet-unit-and-next-old-L-clearing
  - THM-3527-fixed-R7-finite-sheet-unit-and-next-old-L-clearing
scripts:
  - 04-computation/keller_all_level_cleared_norm_packet_arithmetic_audit_20260816.py
  - 04-computation/keller_packet_monoid_branch_transplant_audit_20260816.py
  - 04-computation/keller_all_level_cleared_norm_polynomiality_independent_audit_20260816.py
outputs:
  - 05-knowledge/results/keller_all_level_cleared_norm_packet_arithmetic_audit_20260816.out
  - 05-knowledge/results/keller_packet_monoid_branch_transplant_audit_20260816.out
  - 05-knowledge/results/keller_all_level_cleared_norm_polynomiality_independent_audit_20260816.out
script_sha256:
  - 050d9ef31faa59c7ebb3b4dc0ca4df1774cbc4a10a5dd89c87358d1e73842fb6
  - 8256c7179c415e8588a0612608f7c253baf026d8a237cd8cbf01720758e8b5dc
  - 4e79d10d2a90cfdcbd6948d22b7385da6ea88b923f3483da96448af1ea1cdc77
output_sha256:
  - 2226933130c39b16b74d6805657a1a42acc87830c35282936da571bc62162a26
  - b934628ab80fefe5fdc6662d10f12bffa4f6f1fb25a3006b9254f3f4b7d204d6
  - e7247cb01a61558ee2af6ef662610e946364f6d9d32d44766841a501644ee30b
semantic_sha256:
  - a77811be1a53f2e0d0e0eeac3b4a4ecac358f79c8ed0b6ec5fa6f03f3bb0c826
  - f004cd7643933e81a2fbce73a3df6d72c9e4943702cca7ee17b32d6690be3874
  - 6fe70dcf5a0f1bd4f76ef8bc4986f79be1d74e4a0b6a71e40520dae67f4e456e
hash_basis: LF-normalized bytes
---

# THM-3528 -- all complete packets clear, but finite-sheet defects may remain

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Retain the fixed sporadic Keller map `F:C^3->C^3`, the irreducible target
polynomial

```text
L=27a^2c^2-18abc+16a+b^3c-b^2,
```

and the cubic function-field norm `N` of THM-2473 and THM-3495.  Put
`A=Q[a,b,c]`.  Target variables are renamed `(x,y,z)` after every norm, as in
THM-3506.

## 1. The theorem

Let `P in A` be nonzero and have the complete five-face packet `A(e,m)` of
THM-3506.  Thus its complete maximum-`lambda` face, for `lambda=i-k`, is

```text
in_max-lambda(P)=C x^e(3xz-2y)^m,       C!=0.          (1)
```

At the generic point of `(L)`, let `q_fin` denote the unique regular finite
inverse branch and define

```text
s_L(P)=v_L(P(q_fin)).                                  (2)
```

Then `s_L(P)` is a finite nonnegative integer and

```text
v_L(N(P))=-e+s_L(P),                                   (3)
Q(P):=L^e N(P) belongs to A,                           (4)
ord_L(Q(P))=s_L(P).                                    (5)
```

In particular, the finite-sheet unit condition `s_L(P)=0` is required for
`L`-coprimality, but it is **not** required for polynomiality.  Applying
THM-3522 to the now-unconditional polynomial (4) gives the complete packet

```text
A(e,m) --Q--> A(7e-2m,3e-2m).                         (6)
```

This is a closure theorem for complete packets under one fixed inverse chart,
not a theorem about arbitrary Keller maps.

## 2. The two divergent branches have exact order `-e/2`

Let `R=A_(L)` be the generic DVR of the old boundary.  THM-3498 proves that
the inverse cubic has one regular finite branch and two divergent geometric
branches.  On either divergent branch, with `u=1/w`, one has

```text
x=u^-1,
y=D/S+O(u),
z=-3(D/S)u+O(u^2),
3xz-2y=-11D/S+O(u),
v_L(u)=1/2,                                           (7)
```

where `D/S` is a unit.  A monomial `x^i y^j z^k` therefore has leading
`u`-order `k-i=-lambda`.  The complete face (1) evaluates to a nonzero unit
times `u^-e`; every monomial with smaller `lambda` has strictly larger
valuation.  Hence there is no equal-weight cancellation and

```text
v_L(P(q_div,+))=v_L(P(q_div,-))=-e/2.                 (8)
```

The exponent `m` contributes only the unit in (7).  This explains why the
old-boundary clearing exponent is the first packet coordinate rather than a
new scalar recurrence guessed from earlier rungs.

## 3. The finite branch supplies exactly the residual multiplicity

At the generic point of `(L)`, the coefficients `c,T,S,D` used in THM-3498
are units.  The linear residual root of the inverse cubic is simple, so
Hensel lifting places `q_fin` in the valuation ring of the completed base
field.  Because `P` is a polynomial in the source coordinates,

```text
v_L(P(q_fin))>=0.                                      (9)
```

The value is not identically zero.  Indeed, the generic source coordinate
field is a degree-three field extension of the target field; a nonzero
polynomial on the source is a nonzero field element, and remains nonzero in
each factor after completion.  Thus (2) is a finite integer.

After a splitting extension of the completion, the field norm is the product
of the three branch values.  Valuations add, with the two conjugate divergent
roots each counted once.  Equations (8)--(9) give

```text
v_L(N(P))=(-e/2)+(-e/2)+s_L(P)=-e+s_L(P),             (10)
```

which proves (3).  This calculation includes the local degree/ramification
weight: the pair of half-integral geometric orders is the degree-two local
factor and has integral total order `-e`.

There is also an exact divisor interpretation.  Let `D=V(L)` and let
`C_fin` be the closure in the source of the finite branch over the generic
point of `D`.  That branch has ramification and residue degree one, so

```text
s_L(P)=ord_(C_fin)(P)
      =length_R R/(P(q_fin)).                          (10a)
```

Consequently `s_L(P)>0` if and only if `C_fin` is a component of `V(P)`.
It counts generic primary thickness, not the number of returning reduced
components.

## 4. Global regularity turns the valuation bound into polynomiality

Set

```text
U=Spec(A[L^-1]).                                      (11)
```

THM-2473 identifies `V(L)` as the full nonproperness divisor of `F`; since
`F` is etale, the restriction `F^-1(U)->U` is finite etale of degree three.
The norm of the regular source element `P` is therefore the determinant of a
finite locally free multiplication map and belongs to

```text
N(P) in A[L^-1].                                      (12)
```

There are consequently no hidden `S`, discriminant, chart, or other
denominators.  Since `A` is a UFD and `L` is irreducible, (3) says precisely
that multiplying (12) by `L^e` leaves no denominator.  This proves (4), and
the same valuation identity gives (5).  The norm is nonzero because it is the
field norm of a nonzero field element.

Polynomiality was the only global hypothesis in THM-3522.  Applying that
theorem now proves (6), including all five complete output faces and their
nonzero overlap scalars.

## 5. The all-level raw packet tower

Define the unnormalized cleared-norm tower by

```text
P_0=L,                  (e_0,m_0)=(1,0),
P_(n+1)=L^e_n N(P_n),
(e_(n+1),m_(n+1))=(7e_n-2m_n,3e_n-2m_n).             (13)
```

THM-3506 proves that `L` has `A(1,0)`.  Equations (4) and (6), applied
inductively, prove for every `n>=0` that

```text
P_n is a nonzero polynomial with complete packet A(e_n,m_n).             (14)
```

Multiplying any rung by a nonzero rational scalar changes neither its packet
nor any `L`-valuation.  Suitable canonical scalar choices identify the first
eight rows with the already named polynomials:

| `n` | named scalar normalization | `(e_n,m_n)` |
|---:|---|---:|
| 0 | `L` | `(1,0)` |
| 1 | `H` | `(7,3)` |
| 2 | `J` | `(43,15)` |
| 3 | `G` | `(271,99)` |
| 4 | `R_5` | `(1699,615)` |
| 5 | `R_6` | `(10663,3867)` |
| 6 | `R_7` | `(66907,24255)` |
| 7 | `R_8` | `(419839,152211)` |

The next three actual raw packet rows are therefore

```text
(2634451,955095),
(16530967,5993163),
(103730443,37606575).                                  (15)
```

Equation (15) asserts polynomial complete packets.  It does not name their
irreducible factors or make them image equations.

## 6. Cone, Pell, and Cassini consequences

The seed lies in the stronger invariant cone

```text
e>0,        0<=2m<=e.                                  (16)
```

Writing `(e',m')=(7e-2m,3e-2m)`, one has

```text
m'>=2e>0,       e'-2m'=e+2m>=0.                       (17)
```

Thus (16) propagates.  The congruences `3|m` and `e=1 mod 3` also propagate,
so every packet in (13) is admissible.  Since the renewal matrix has trace
`5` and determinant `-8`, both coordinates obey

```text
u_(n+2)=5u_(n+1)+8u_n.                                (18)
```

The Pell-57 and Cassini identities formerly proved only for the abstract
matrix orbit are now realized by the raw polynomial tower for every `n`:

```text
3e_n^2-9e_nm_n+2m_n^2=3(-8)^n,                       (19)
e_nm_(n+1)-m_ne_(n+1)=3(-8)^n.                       (20)
```

Equivalently, with `X_n=6e_n-9m_n` and `Y_n=m_n`, multiplication by
`(5+sqrt(57))/2` gives

```text
X_n^2-57Y_n^2=36(-8)^n.                              (21)
```

The fibre-degree grading `3^n` and the packet-norm grading `(-8)^n` remain
different gradings sharing one iteration index.

## 7. The sharp boundary: polynomial packet versus new image prime

For the tower (13), put

```text
s_n=s_L(P_n).                                         (22)
```

Equation (5) becomes the exact divisor-return ledger

```text
ord_L(P_(n+1))=s_n.                                   (23)
```

The finite-sheet computations of THM-3498, THM-3506, THM-3521, THM-3523,
and THM-3527 prove `s_n=0` through input `R_7` and output `R_8`.  They remain
strictly stronger than this theorem: they prove the corresponding cleared
norm is coprime to `L`.  Starting with input `R_8`, the values `s_n` are open
and may be positive.  If that occurs, the next raw packet remains polynomial
but carries an old-`L` factor of exactly that multiplicity.

More precisely, if `Q(P)=L^sR_P` with `L` not dividing `R_P`, then

```text
div(Q(P))=sD+div(R_P),       (Q(P)):L^infinity=(R_P). (23a)
```

Thus `s>0` is a return of the finite source component to the old target
divisor.  It is not automatically a return of the whole hypersurface: that
requires the reduced source hypersurface to have no other components.  Nor
does `s=1` make `Q(P)` a reduced image equation; source reducedness, generic
image degree, component collisions, squarefreeness, and image-support
identification are separate gates.  The hostiles `P=t(t-1)`, `P=t^2` under
the identity map distinguish component return, whole-hypersurface return, and
multiplicity, while a generic degree-two finite map can send a reduced prime
to a squared norm.

Accordingly this theorem proves neither:

- `L`-coprimality or finite-sheet units after `R_8`;
- irreducibility, squarefreeness, primitivity over `Z`, or newest-image status
  for later raw rungs;
- separability, full degree, discriminant square classes, or distinct
  nonproperness components at later iterate levels;
- an all-level prime-factor or Jelonek-component classification;
- an analogous closure law for arbitrary Keller maps, a classification of
  Keller counterexamples, `JC(2)`, `DC(2)`, LRC, or the general Jacobian
  conjecture.

In particular, the raw all-level polynomial tower is not an all-level tower
of new image primes.

## 8. Packet monoid and the branch-transplant law

Complete packets form a multiplicative monoid modulo nonzero rational
scalars.  Indeed, every weighted initial form of a product is the product of
the two initial forms, so

```text
A(e,m) * A(f,n) = A(e+f,m+n).                         (24)
```

On this monoid define the fixed cleared-norm operator

```text
T(P)=L^e N(P)        when P has packet A(e,m).         (25)
```

The polynomiality theorem makes `T` an everywhere-defined endomorphism, and
norm multiplicativity plus grade additivity gives the exact laws

```text
T(PQ)=T(P)T(Q),       T(cP)=c^3T(P),                  (26)
grade(T(P))=M grade(P),
M=[[7,-2],[3,-2]].                                    (27)
```

Now suppose an old-boundary return occurs at a later raw rung:

```text
P_n=L^s R,       s=ord_L(P_n)>0.                      (28)
```

Initial forms are multiplicative, and the five initial forms of `L` are the
packet `A(1,0)`.  Dividing each complete face of `P_n` by the corresponding
face of `L^s` proves

```text
R has A(e_n-s,m_n),       0<s<=e_n-m_n.               (29)
```

The last bound is sharp at the packet level: the `z^(e-m)` exponent on the
minimum-`beta` face is the first one to reach zero.  Applying (26) repeatedly
to (28) gives the branch-transplant law

```text
P_(n+k)=P_k^s T^k(R)        for every k>=0.            (30)
```

Thus a returned `L^s` factor is carried along the canonical ancestry:

```text
L^s, H^s, J^s, G^s, R_5^s, ...                       (31)
```

up to the already recorded nonzero rational normalizations.  More precisely,
`s_j>0` puts `L^s` in the next rung `P_(j+1)`.  For every genuine later
output, `m_(j+1)>0`, so the quotient in (29) is nonconstant; that next rung
and all its descendants are reducible.  The seed factorization `P_0=L=L*1`
is the sharp hostile showing that positive multiplicity does not imply
reducibility when the cofactor is a unit.  The converse is false:
`s_j=0` proves only that `P_(j+1)` is `L`-coprime and does not prove its
irreducibility.

The defect word has a faithful formal encoding

```text
D(t)=sum_(n>=0) s_n t^n in N[[t]],                    (32)
```

and one ancestry shift is multiplication by `t`.  Its support is literally a
subset of the natural numbers.  Replacing (32) by a scalar harmonic subseries
loses that support: already on `{1,...,13}`, the `8192` subsets give only
`3712` reciprocal sums, `2944` of them collision values, with maximum
multiplicity three; in particular

```text
1/2=1/3+1/6=1/4+1/6+1/12.                            (33)
```

The formal series, not the harmonic scalar, is the correct ancestry sidecar.

## 9. Exact arithmetic and independent hostile audits

The first companion checks the matrix determinant, the eight named packet rows,
the next three raw rows, invariant cone, congruences, order-two recurrence,
Pell identity through `n=15`, Cassini identity through `n=14`, and the
symbolic valuation ledger

```text
(-e+s)+e=s>=0
```

for hostile values `s=0,1,2,e,e+1`.  Ordinary and optimized replays agree
line-for-line with the stored output.  These checks audit the arithmetic and
scope; Sections 2--4 are the geometric proof.  The second companion checks
all five packet-face exponent vectors under 256 orbit-row products, renewal
additivity, the sharp `L^s` quotient boundary through ten rows, three
synthetic transplant hostiles through five descendants, and the complete
`2^13` harmonic collision census.  Synthetic defects test (30); they are not
observed returns in the fixed Keller tower.  The independent implementation
starts from the reciprocal cubic, verifies the Newton slopes `-1/2,0`, both
local norm-weight conventions, `gcd(L,cTSD)=1`, the canonical finite branch,
and the nonmonic hostile

```text
Res_w(E,w)=2c,       N(w)=2c/L.                       (34)
```

It then checks all `15,251` admissible packets with `e<=300` and independently
returns `ACCEPT all-level raw polynomial packet induction at fixed-map scope`.

Reproduce with

```text
python -B 04-computation/keller_all_level_cleared_norm_packet_arithmetic_audit_20260816.py
python -B -O 04-computation/keller_all_level_cleared_norm_packet_arithmetic_audit_20260816.py
python -B 04-computation/keller_packet_monoid_branch_transplant_audit_20260816.py
python -B -O 04-computation/keller_packet_monoid_branch_transplant_audit_20260816.py
python -B 04-computation/keller_all_level_cleared_norm_polynomiality_independent_audit_20260816.py
python -B -O 04-computation/keller_all_level_cleared_norm_polynomiality_independent_audit_20260816.py
```

**QED.**
