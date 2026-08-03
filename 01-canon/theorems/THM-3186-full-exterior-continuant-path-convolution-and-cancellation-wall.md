---
id: THM-3186
title: "Full exterior continuant path convolution and cancellation wall"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT IMMUTABLE AUDIT PENDING.
  Every arbitrary-length hidden-to-visible coefficient of the scalar
  factorial exterior transfer is an exact exit-time convolution of a Jacobi
  continuant with the invariant-line tail.  The formula includes the
  length-one zero boundary and recovers THM-3183 at length two.  At length
  three, two individually nonzero unit paths can cancel exactly, so the
  transfer graph and every local Smith profile still do not determine
  projected visibility.
audit: >
  The exact symbolic companion reconstructs the exterior transfer, compares
  direct matrix products with the convolution through seven steps, verifies
  an independent monomer-dimer expansion of every continuant, and checks the
  length-one and length-two boundaries.  It then verifies an exact factorial
  two-path cancellation whose transfers and path weights are all 11-adic
  units.  Normal/-O/stored replay and independent immutable audit are pending.
source: root/multiscale-newton-flag/2026-08-02
depends_on:
  - THM-3183-factorial-hecke-lattice-square-and-oriented-wedge-continuant
related:
  - THM-3176-six-step-prime-resonance-third-euclidean-newton-separation
  - THM-3185-iterated-factorial-frobenius-descent-and-witt-carry-reset-hierarchy
script: 04-computation/factorial_exterior_path_convolution_thm3186.py
output: 05-knowledge/results/factorial_exterior_path_convolution_thm3186.out
script_sha256: 5ee9c3924113be5e6c9769d2c2bf67dbffae3d63d4917101f86278b1c3e8813c
output_sha256: 8966f0979ed97b4acfc0105e9781f7cdacbe654c9f8ae6729c50166aa9559754
hash_basis: LF-normalized bytes
---

# THM-3186 -- full exterior continuant path convolution and cancellation wall

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT IMMUTABLE AUDIT PENDING.**

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
```

and compare LF-normalized bytes with the declared output.  The companion uses
exact symbolic and rational arithmetic only.  It reconstructs `(3)` from all
`2x2` minors of `S_i`, compares direct products with `(12)--(14)` through
seven steps, and independently enumerates the monomer-dimer tilings in
`(10)`.  It also verifies `(15)--(18)` and every rational and `11`-adic claim
in `(23)--(27)`.  There is no floating point, random sampling, imported
executable, or assertion-sensitive test.

The theorem proves no general PRS-depth law, fixed-offset closure, `NC(2)`,
`GMC(2)`, or `LRC(14)` consequence.

**QED.**
