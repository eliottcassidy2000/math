---
id: THM-3458
title: "Rule 30 right-edge 2-adic odometer and the moving-observer boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Right-edge packing
  conjugates the Rule 30 single-seed evolution to a compatible tower of
  unitriangular Boolean permutations.  Its seed-orbit closure is a 2-adic
  odometer and every fixed edge-offset observable is periodic, while the
  prize center is given by a moving diagonal observer rather than by a fixed
  edge-coordinate readout; an exact width-six collision proves real quotient loss.
  Periodic-ring and de Bruijn compilers are exact but settle none
  of the three Rule 30 prizes.
source: root-rule30-260813291-20260815
audit: >
  independent recurrence, quotient-bijection, period-lift, inverse-limit,
  moving-observer, ring-prefix, de Bruijn, circulant-Jacobian, prize-scope,
  hash, and ordinary/optimized/stored byte-replay audit; documentation clean
depends_on: []
related:
  - THM-2005-support-dirichlet-automatic-tournament-atlas
  - THM-2050-period14-top-germs-do-not-determine-global-loneliness
  - THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar
  - THM-3364-cyclotomic-boolean-clocks-berggren-t4-xor-and-crt-phase
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
  - THM-3459-rule30-ternary-intersection-factorial-truth-lift-and-keller-boundaries
script: 04-computation/rule30_right_edge_odometer_thm3458.py
output: 05-knowledge/results/rule30_right_edge_odometer_thm3458.out
script_sha256: 8b9a6d029419079f5507d3b153fc43af760e846ddbc93d199392ff4da81640ec
output_sha256: ae7f2f6dd2c21bb3e4fc6d9d2779080477563d74b169d1a2b7123f062f0ac3d3
hash_basis: raw bytes
---

# THM-3458 -- Rule 30 right-edge 2-adic odometer and the moving-observer boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem solves an exact edge subsystem of the distinguished Rule 30
seed orbit.  Its principal boundary is equally important: the center column
is presented as a moving coordinate rather than one of the fixed edge-
coordinate observables.

## 1. Packed right-edge recurrence

Let `a_t(j)` be the standard Rule 30 evolution from the single seed

```text
a_0(0)=1,       a_0(j)=0 for j!=0,                     (1)
```

using

```text
f(l,c,r)=l xor c xor r xor c r.                         (2)
```

Read row `t` inward from its right edge and pack it as

```text
R_t=sum_(k=0)^(2t) a_t(t-k) 2^k.                       (3)
```

Then

```text
R_0=1,
R_(t+1)=Phi(R_t),
Phi(x)=x xor ((2x) or (4x)).                            (4)
```

In particular, if `x_k` is bit `k` and exterior negative bits are zero,

```text
(Phi x)_k=x_k xor (x_(k-1) or x_(k-2)).                 (5)
```

The Rule 30 center bit is the diagonal readout

```text
c_t=a_t(0)=bit_t(R_t).                                  (6)
```

Equation (4) follows by applying (2) at every site and matching the three
parents after the right-edge shift.  The first packed rows are

```text
1, 7, 25, 111, 401, 1783, 6409, 28479, 102849,... .    (7)
```

The highest set bit of `R_t` is exactly `2t`, since the extreme-left boundary
bit propagates uniquely and remains one.

## 2. Finite triangular permutations and the period-lift law

For `w>=1`, let `Phi_w` be (4) modulo `2^w`.  Equation (5) is unitriangular:
bit `k` of the output is bit `k` of the input plus a function of lower bits.
It follows that `Phi_w` is a permutation of `Z/2^w Z`, with inverse recovered
from low to high bits.

Let `P_w` be the least positive seed return time

```text
Phi_w^(P_w)(1)=1.                                       (8)
```

Then

```text
P_1=1,
P_w | P_(w+1),
P_(w+1)/P_w in {1,2},                                  (9)
```

so every `P_w` is a power of two.  More exactly, put

```text
epsilon_w=bit_w(R_(P_w)).                               (10)
```

Then

```text
P_(w+1)=2^(epsilon_w) P_w.                              (11)
```

Proof: after `P_w` steps the lower `w` bits return.  On the two lifts of the
finite seed cycle, the new high bit is translated by the fixed bit
`epsilon_w` after one circuit.  It either returns immediately or toggles and
returns after the second circuit, giving (9)--(11).

The periods are unbounded.  Indeed, if `2P_w<w`, the top bit
`2P_w` of `R_(P_w)` lies inside the width-`w` quotient and is one, whereas the
same bit of the seed is zero.  Hence

```text
ceil(w/2) <= P_w <= 2^(w-1).                            (12)
```

The upper bound is the size of the seed fibre modulo two, and the lower bound
is enough to make the period tower cofinal among powers of two.

## 3. The seed-orbit closure is a 2-adic odometer

The formulas (5) define a compatible homeomorphism `Phi` of `Z_2`.  Let

```text
X=closure{Phi^t(1):t>=0} in Z_2.                        (13)
```

At level `w`, let `C_w` be the embedded seed cycle of length `P_w`.  Then

```text
X = inverse_limit C_w  ~=  inverse_limit Z/P_w Z.       (14)
```

Because the `P_w` are unbounded powers of two, this inverse limit is `Z_2`.
The compatible maps

```text
t mod P_w |-> Phi_w^t(1)                               (15)
```

give a homeomorphism

```text
iota:Z_2 -> X,          iota(t+1)=Phi(iota(t)).         (16)
```

Thus `Phi|X` is conjugate to the 2-adic adding machine.  This is an orbit-
closure statement, not a claim that `Phi` is ergodic on all of `Z_2`.

## 4. Fixed observers are periodic; the center observer moves

For fixed `k`, define the edge-offset column

```text
b_k(t)=bit_k(R_t).                                      (17)
```

It factors through `Phi_(k+1)`, so it is purely periodic with period dividing
`P_(k+1)`.  If one period is `v_0,...,v_(p-1)`, then

```text
sum_(t>=0) b_k(t)z^t
  =(sum_(r=0)^(p-1) v_r z^r)/(1-z^p).                  (18)
```

For example,

```text
b_5 = (00010101)^infinity,
sum_(t>=0)b_5(t)z^t=(z^3+z^5+z^7)/(1-z^8),             (19)
```

whose period density is `3/8`, not `1/2`.

The prize sequence (6) is not any fixed `b_k`: it reads the diagonal entry
`b_t(t)` of the periodic-column array.  Its definition therefore requires the
external observer address `t`, which a fixed truncation does not itself carry.
There is also an exact failure of naive state-only recovery at width six:

```text
R_7 = R_15 mod 2^6 = 63,
c_7=0,                     c_15=1.                     (20)
```

Hence identical width-six odometer state can give different center bits when
the observation depth moves.  This single collision is not extrapolated to
all widths: proving that every fixed quotient fails even coincidentally would
already bear on prize nonperiodicity.  What is proved generally is the type
separation between fixed continuous observables and the moving readout.
Fixed-observer periodicity, rational OGFs, and Haar properties of the odometer
do not decide eventual periodicity or density of the diagonal sequence.

This is the same general warning as a phase/owner sidecar in the LRC work:
closing the state dynamics does not close a target predicate whose observer
address grows with time.

## 5. Periodic-ring finite-prefix compilers

On a ring of length `N`, Rule 30 is a map of a finite set.  Every center orbit
is eventually periodic.  If its preperiod is `mu`, period is `lambda`, and
the center bits are `d_t`, its OGF is exactly

```text
sum_(t>=0)d_t z^t
 =sum_(t<mu)d_t z^t
  +z^mu (sum_(r<lambda)d_(mu+r)z^r)/(1-z^lambda).       (21)
```

Lift the ring seed periodically to the line.  Its nearest extra seeds are at
distance `N`, so the ring center agrees with the isolated-seed center for all

```text
0<=t<N.                                                 (22)
```

Consequently every finite prefix of the prize sequence has a genuine
eventually periodic Rule 30 ring extension.  No finite prefix, fitted
recurrence, or fitted density can by itself settle either of the first two
prizes.

## 6. De Bruijn inverse-count sequence compiler

Use the ordered states `00,01,10,11` and let

```text
A_0 = [[1,0,0,0],
       [0,0,0,0],
       [0,1,0,0],
       [0,0,1,1]],

A_1 = [[0,1,0,0],
       [0,0,1,1],
       [1,0,0,0],
       [0,0,0,0]].                                     (23)
```

An edge `(a,b)->(b,c)` is counted in `A_e` exactly when
`f(a,b,c)=e`.  Therefore a cyclic output word
`y=y_0...y_(N-1)` has exactly

```text
tr(A_(y_0) A_(y_1) ... A_(y_(N-1)))                   (24)
```

periodic predecessors.

For a fixed block `w`, put `M_w=product_(e in w) A_e`.  The number of cyclic
predecessors of `w^k` is

```text
q_k=tr(M_w^k).                                          (25)
```

By Cayley--Hamilton, `(q_k)` is C-finite of order at most four and has the
usual eigenvalue/characteristic-polynomial closed form.  Arbitrary one-step
cyclic inverse counts require only `O(N)` sparse `4 x 4` matrix
multiplications.  This is an exact spatial inverse-count compiler, not an
algorithm for the time-`n` center bit.

Two boundary cases are

```text
#preimages(0^N)=2,
#preimages(1^N)=3 if 3|N, and 0 otherwise.              (26)
```

The first predecessors are `0^N,1^N`; the second are the three rotations of
`(100)^(N/3)`.  Thus open-block surjectivity and periodic gluing are different
typed problems.

Equivalently,

```text
char(A_0)(lambda)=lambda^2(lambda-1)^2,
char(A_1)(lambda)=lambda(lambda-1)(lambda^2+lambda+1),   (26a)
```

so taking traces of powers proves both formulas in (26) for every `N`.

## 7. Characteristic-zero Jacobian boundary

Lift the periodic Rule 30 ANF to characteristic zero:

```text
F_i=x_(i-1)+x_i+x_(i+1)+x_i x_(i+1),                  (27)
```

with cyclic indices.  At the zero point, its Jacobian determinant is zero
when `3|N` and has absolute value three otherwise.  At the all-minus-one
point its determinant has absolute value one.  Hence the determinant is
nonconstant for every `N>=3`; the periodic ANF lift is not Keller.

Indeed, at zero the Jacobian is the circulant `I+S+S^(-1)`, so

```text
det J(0)=product_(zeta^N=1)(1+zeta+zeta^(-1)).          (27a)
```

This vanishes exactly when a primitive cube root is an `N`th root and has
absolute value three otherwise.  At the all-minus-one point the Jacobian is
the cyclic shift `S^(-1)`, whose determinant has absolute value one.

By contrast, the open-edge recurrence (5), read as an ANF polynomial map, is
unit triangular and hence a tame Keller automorphism.  Closing the boundary
turns a triangular dependency order into a cycle and destroys that property.
The faithful real multilinear lift and the characteristic-two Frobenius
representative have the further incompatible behaviours recorded in
THM-3459.  Boundary condition, characteristic, and polynomial representative
are all load-bearing sidecars.  No planar Jacobian-conjecture reduction
follows.

## 8. Relation to the Rule 30 prizes

The three advertised questions concern the center column of the isolated
single-seed evolution: eventual periodicity, limiting density one half, and a
specified lower bound on computation.  This theorem supplies:

```text
PROVED: fixed-edge odometer and fixed-offset rational sequences;
PROVED: moving-observer typing and the width-six obstruction (20);
PROVED: finite-ring and spatial inverse-count compilers;
OPEN:   all three center-column prize questions.                        (28)
```

The packed recurrence computes the first `n` rows by `n` iterations on
growing integers.  It is an exact upper-bound representation only.  It proves
no lower bound in any machine or bit-cost model.  In particular, no inference
is made from the announcement's informal big-O phrasing.

## 9. Preservation and loss ledger

```text
source:     isolated-seed Rule 30 space-time diagram
target:     right-edge integers; finite triangular cycles; odometer closure;
            periodic rings; de Bruijn inverse counts
map:        (3)--(6), quotient mod 2^w, (21), (23)--(25)
preserves:  every packed row; fixed lower bits; finite seed-cycle phase;
            ring prefixes before wrap; exact one-step predecessor counts
destroys:   center observer address under fixed quotient; infinite boundary
            under a ring; forward time history under inverse-count matrices
sidecar:    moving depth k=t; isolated boundary; forward seed phase;
            characteristic and lift representative
consequence: odometer/rational/C-finite compilers and exact no-gos, not a
             Rule 30 prize, LRC(14), FC(3), or JC(2) result
```

## 10. Exact audit

The companion compares direct and packed rows through time `512`, inverts
every state through width twelve, computes seed periods through width thirty,
checks the lift bits (10), freezes fixed-offset words through offset twelve,
verifies (20), audits single-seed rings `2<=N<=12`, checks the de Bruijn
formula on every binary output word through length eight, and evaluates the
periodic Jacobian controls through `N=14`.

```bash
python3 04-computation/rule30_right_edge_odometer_thm3458.py
python3 -O 04-computation/rule30_right_edge_odometer_thm3458.py
```

Both runs reproduce the stored output byte for byte.  The finite period prefix
is evidence for the implementation only; the recurrence (11), rather than a
fit to that prefix, proves the universal period-tower statement.
