---
id: THM-3186
title: "Full exterior continuant path convolution and cancellation wall"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Every arbitrary-length hidden-to-visible coefficient of the scalar
  factorial exterior transfer is an exact exit-time convolution of a Jacobi
  continuant with the invariant-line tail.  The formula includes the
  length-one zero boundary and recovers THM-3183 at length two.  At length
  three, two individually nonzero unit paths can cancel exactly, so the
  transfer graph and every local Smith profile still do not determine
  projected visibility.  In the factorial specialization the transverse
  continuant has a shifted order-two normalization with a first-order
  differential OGF, while the full selected visibility sequence satisfies an
  explicit cleared order-three polynomial-coefficient recurrence.  Hence both
  tails are P-recursive and evaluable in linear arithmetic time for fixed
  parameters, without a C-finite or PRS-depth conclusion.
audit: >
  The exact symbolic companion reconstructs the exterior transfer, compares
  direct matrix products with the convolution through seven steps, verifies
  an independent monomer-dimer expansion of every continuant, and checks the
  length-one and length-two boundaries.  It then verifies an exact factorial
  two-path cancellation whose transfers and path weights are all 11-adic
  units, and directly compares all three scalar and exterior Boolean support
  patterns with the same-Smith positive control.  It also checks seven terms
  of the normalized continuant, eight coefficients of its differential OGF,
  and five generic instances of the cleared visibility recurrence.  A second
  companion imports neither SymPy nor primary code: its custom sparse-
  polynomial ring and Fraction matrices independently repeat those recurrence
  checks, including eight indices across an actual vanishing exit coefficient.
  Normal/-O/stored replay and the independent immutable audits pass.
source: root/multiscale-newton-flag/2026-08-02
depends_on:
  - THM-3183-factorial-hecke-lattice-square-and-oriented-wedge-continuant
related:
  - THM-3176-six-step-prime-resonance-third-euclidean-newton-separation
  - THM-3185-iterated-factorial-frobenius-descent-and-witt-carry-reset-hierarchy
  - THM-3202-c3-repeated-join-moving-jet-formula-and-cfinite-obstruction
script: 04-computation/factorial_exterior_path_convolution_thm3186.py
output: 05-knowledge/results/factorial_exterior_path_convolution_thm3186.out
script_sha256: 5bacf9a0d7da1e467f9f22fc0d21d0c0b1968dab955736440e9010fb00e21eff
output_sha256: 66d600435bdb58b0fd0b3c82c49c7f9b2154b51e624a12e8b6bb3424ccac5776
independent_script: 04-computation/factorial_exterior_precursive_visibility_independent_audit_thm3186.py
independent_output: 05-knowledge/results/factorial_exterior_precursive_visibility_independent_audit_thm3186.out
independent_script_sha256: 05b7c8dc18693fa971c83af3bae6bd9a3bba2c18ca4184c309b210faafb58288
independent_output_sha256: e453bde9f739ecf2f3f4eb9470f7dc21365b51ba4b1f220f6c18fd2b917a0254
hash_basis: LF-normalized bytes
---

# THM-3186 -- full exterior continuant path convolution and cancellation wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3183 found the first two-step route from the hidden scalar wedge to the
visible moment wedge.  The whole finite tail admits a closed formula.  It is
not a product of local pivots: it is a signed sum over possible exit times,
and distinct nonzero paths can cancel.

## 1. Exterior transfer and notation

Retain THM-3183's scalar factorial state and write

```text
S_i=[a_i b_i c_i; 1 0 0; 0 0 d],
a_i=2(i+1)(2i+1)v,
b_i=i(i+1)Delta,
c_i=d-i-1,                    Delta=1-4dv.                 (1)
```

In the ordered wedge basis

```text
e_0=M_i wedge M_(i-1),
e_1=M_i wedge D_i,
e_2=M_(i-1) wedge D_i,                                      (2)
```

put

```text
A_i=Lambda^2 S_i
   =[u_i  -c_i    0    ]
    [0     alpha_i beta_i]
    [0     d       0    ],                                  (3)

u_i=-b_i,              alpha_i=a_i d,       beta_i=b_i d.   (4)
```

The argument below is formal over any commutative ring for matrices of the
shape `(3)`.  The visibility equivalences are stated over an integral domain
so that zero divisors do not create extra walls.

## 2. Transverse continuant

Fix `n>=1`.  Define

```text
t_0=0,                    t_1=beta_n,                       (5)
t_(j+1)=alpha_(n+j)t_j+d beta_(n+j)t_(j-1)       (j>=1).   (6)
```

Equivalently, let

```text
C_0=1,                    C_1=alpha_(n+1),                  (7)
C_r=alpha_(n+r)C_(r-1)+d beta_(n+r)C_(r-2)       (r>=2).   (8)
```

Then

```text
t_j=beta_n C_(j-1)                         (j>=1).          (9)
```

The polynomial `C_r` is the matching polynomial of a weighted path on the
vertices `1,...,r`: a monomer at `i` has weight `alpha_(n+i)`, and a dimer
joining `i-1` to `i` has weight `d beta_(n+i)`.  Thus

```text
C_r=sum_(monomer-dimer tilings T of [1,r]) weight(T).       (10)
```

This is not merely terminology.  Splitting the tilings according as the last
vertex is a monomer or belongs to the last dimer gives exactly `(8)`.

### Factorial normalization and differential generating law

The special coefficients `(4)` remove almost all of the apparent quadratic
growth in `(8)`.  Let `(n+2)^(overline r)` denote the rising factorial and
define

```text
G_(-1)=0,                 G_0=1,
G_r=2v(2n+2r+1)G_(r-1)+Delta G_(r-2)            (r>=1).   (N1)
```

Then the factorial continuant is exactly

```text
C_r=d^r (n+2)^(overline r) G_r.                           (N2)
```

Indeed, the monomer term in `(8)` contributes the missing factor `n+r+1`,
while the dimer term contributes the two missing factors
`(n+r)(n+r+1)`; after removing `d^r(n+2)^(overline r)`, they become the two
coefficients in `(N1)`.

The formal ordinary generating series

```text
G(z)=sum_(r>=0) G_r z^r
```

satisfies the first-order differential equation

```text
4v z^2 G'(z)
 +[2v(2n+3)z+Delta z^2-1]G(z)+1=0.                       (N3)
```

Coefficient extraction in `(N3)` is exactly `(N1)`.  Thus the normalized
continuant is P-recursive of order at most two and its formal OGF is D-finite.
The equation is a closed differential representation; it need not define a
convergent analytic OGF when the coefficients have factorial growth.

## 3. Full hidden-to-visible convolution

For `L>=1`, set

```text
A_[n,L]=A_(n+L-1)...A_(n+1)A_n                             (11)
```

and let `V_L` be the `e_0` coordinate of `A_[n,L]e_2`.  Then

```text
A_[n,L]e_2=(V_L,t_L,d t_(L-1))^t,                          (12)
```

where

```text
V_L=-sum_(j=1)^(L-1)
       c_(n+j)t_j prod_(h=j+1)^(L-1)u_(n+h)                (13)

   =-beta_n sum_(j=1)^(L-1)
       c_(n+j)C_(j-1) prod_(h=j+1)^(L-1)u_(n+h).           (14)
```

Empty sums are zero and empty products are one.  In particular,

```text
V_1=0,                                                     (15)
V_2=-c_(n+1)beta_n,                                       (16)
V_3=-beta_n[c_(n+1)u_(n+2)+c_(n+2)alpha_(n+1)].           (17)
```

Equation `(16)` is exactly THM-3183's pivot

```text
(n+2-d)n(n+1)Delta d.                                     (18)
```

### Proof

After one step, `(3)` sends `e_2` to `(0,beta_n,0)^t`, giving
`(12)` for `L=1`.  If `(12)` holds after `j` steps, one more multiplication
gives

```text
(V_(j+1),t_(j+1),d t_j)
 =(u_(n+j)V_j-c_(n+j)t_j,
   alpha_(n+j)t_j+beta_(n+j)d t_(j-1),
   d t_j).                                                  (19)
```

The last two coordinates are `(6)`.  Iterating the first coordinate from
`V_1=0` yields `(13)`, and `(9)` yields `(14)`.  This proves every asserted
identity.

There is also a directed-path proof.  The source must first take the edge
`e_2 -> e_1` of weight `beta_n`.  Before leaving the transverse pair it makes
a monomer-dimer walk counted by `C_(j-1)`.  It exits at the unique time `j`
through `e_1 -> e_0` of weight `-c_(n+j)`, then remains on the invariant line
through the weights `u_(n+j+1),...,u_(n+L-1)`.  Summing over the mutually
exclusive exit times gives `(14)`.

The cancellation wall is only a selected-chart wall on the local-unit
locus.  Indeed, put `q_j=(t_j,d t_(j-1))^t`.  Then

```text
q_(j+1)=[alpha_(n+j) beta_(n+j); d 0]q_j,
det=-d beta_(n+j).                                         (19a)
```

If `d beta_n prod_(j=1)^(L-1) beta_(n+j)` is nonzero, these `2x2`
transverse steps are invertible over the fraction field and
`q_1=(beta_n,0)^t` is nonzero.  Hence `q_L` cannot vanish.  In particular,
`V_L=0` never means death of the full exterior state on this locus: at least
one of the complementary `e_1,e_2` charts remains visible.  This is the
precise scalar-frame form of a Pluecker-chart transition, not yet an
identification with a coefficient-row PRS chart.

### Cleared P-recursive visibility law

Put `V_0=0`, the visible coordinate before any transfer.  The first two
coordinates in `(19)` give

```text
u_(n+j)V_j-V_(j+1)=c_(n+j)t_j,                 (P1)
t_(j+1)=alpha_(n+j)t_j+beta_(n+j)d t_(j-1).    (P2)
```

Eliminating `t_j,t_(j-1),t_(j+1)` and clearing the potentially vanishing
coefficients `c_i` yields, for every `j>=1`,

```text
c_(n+j)c_(n+j-1)[u_(n+j+1)V_(j+1)-V_(j+2)]
 =c_(n+j+1)c_(n+j-1)alpha_(n+j)
    [u_(n+j)V_j-V_(j+1)]
  +c_(n+j+1)c_(n+j)beta_(n+j)d
    [u_(n+j-1)V_(j-1)-V_j].                    (P3)
```

No division by a possibly vanishing `c_i` occurs, so `(P3)` remains valid at
the boundary indices where an exit edge disappears.  Under `(4)`, every
coefficient is a polynomial in `j`, and the coefficient of `V_(j+2)` is
`-c_(n+j)c_(n+j-1)`, a nonzero quadratic with leading coefficient `-1`.
Consequently `(V_L)_(L>=0)` is P-recursive of order at most three over the
parameter field.  For fixed rational `n,d,v`, either `(19)` or `(P3)` computes
the first `L` values in `O(L)` arithmetic operations, using the transfer at
isolated singular indices of `(P3)`.

This computational closure coexists with the exact selected-chart
cancellation below.  It neither asserts C-finiteness nor identifies the
coefficient-degree PRS transfer.

## 4. The complete cancellation wall

Over an integral domain, the source is visible after `L` steps exactly when

```text
beta_n E_L != 0,                                           (20)

E_L=sum_(j=1)^(L-1)
      c_(n+j)C_(j-1)prod_(h=j+1)^(L-1)u_(n+h).             (21)
```

The first alternative `beta_n=0` is a local entry wall.  The second,

```text
E_L=0,                                                     (22)
```

is a genuine signed path-cancellation wall.  For `L>=3` it need not be the
union of any local zero sets.

The smallest exact factorial hostile occurs already at length three.  Take

```text
n=1,                   d=5,                   v=4/105,
Delta=5/21.                                                 (23)
```

Then

```text
beta_1=50/21,
c_2=2,                u_3=-20/7,
c_3=1,                alpha_2=40/7.                         (24)
```

Both exit paths are nonzero, but their visible contributions are

```text
-beta_1 c_2 u_3= 2000/147,
-beta_1 c_3 alpha_2=-2000/147.                              (25)
```

Consequently

```text
V_2=-100/21 !=0,                  V_3=0.                    (26)
```

This is not caused by singular local dynamics.  Over `Z_11`, every number
displayed in `(24)--(25)` is a unit, and `S_1,S_2,S_3` (hence also their
exterior squares) are unimodular.  For the same `n=1,d=5` but `v=1`, the
three scalar and exterior transfers have the identical unimodular Smith
profiles and the same nonzero support graph, whereas

```text
V_3=115140                                                   (27)
```

is an `11`-adic unit.  Thus the local Smith profiles and transfer graph do
not determine projected visibility.

## 5. What the path formula changes

The generic transfer support graph has a persistent route from hidden to
visible after graph distance two.  Equations `(14)` and `(25)` show that this
graph is only the unsigned carrier.  The invariant needed for projected
Newton data is the oriented exit-time amplitude `E_L`, including its
arithmetic cancellations.  In holotopy language, paths may persist while
their signed section dies on a cancellation hypersurface; no internal
topology of the support graph detects that event.

The formula is therefore a complete finite-tail sidecar for the bare scalar
exterior transfer, but it is not yet the coefficient-degree PRS transfer.
THM-3183 proves that the offset-six PRS walls `H,J,K` live in additional
Schur-complement and connection data.  Identifying those data with a gauged
version of `(14)`, and proving the empirical first-separation depth
`floor(s/2)`, remain open.

## 6. Exact evidence and scope

Run

```text
python 04-computation/factorial_exterior_path_convolution_thm3186.py
python -O 04-computation/factorial_exterior_path_convolution_thm3186.py
python 04-computation/factorial_exterior_precursive_visibility_independent_audit_thm3186.py
python -O 04-computation/factorial_exterior_precursive_visibility_independent_audit_thm3186.py
```

and compare LF-normalized bytes with the declared output.  The companion uses
exact symbolic and rational arithmetic only.  It reconstructs `(3)` from all
`2x2` minors of `S_i`, compares direct products with `(12)--(14)` through
seven steps, and independently enumerates the monomer-dimer tilings in
`(10)`.  It also verifies `(N1)--(N3)`, five free-coefficient instances of
`(P3)`, `(15)--(18)`, and every rational and `11`-adic claim
in `(23)--(27)`.  There is no floating point, random sampling, imported
executable, or assertion-sensitive test.  The independent companion repeats
the recurrence arguments in a custom sparse-polynomial ring, checks that the
leading `V_(j+2)` coefficient has degree two and leading coefficient `-1`, and
uses `Fraction` matrices to replay eight indices across the actual boundary
`n=1,d=5,c_(n+3)=0`, including both hostile controls.

The theorem proves no general PRS-depth law, fixed-offset closure, `NC(2)`,
`GMC(2)`, or `LRC(14)` consequence.

**QED.**
