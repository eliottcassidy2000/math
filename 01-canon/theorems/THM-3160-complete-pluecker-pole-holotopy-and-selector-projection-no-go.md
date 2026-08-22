---
id: THM-3160
title: "Complete Pluecker pole holotopy and same-degree selector projection no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Common virtual-pole subtraction acts functorially on the full exterior
  square of the two endpoint profiles.  The usual selector currents are only
  same-degree Pluecker rows and are not an invariant projection: two
  one-letter endpoint pairs have the identical zero selector profile in every
  degree, while one common pole subtraction gives different degree-two
  currents.  Cross-degree endpoint minors restore an exact linear depth map,
  but neither Hasse positivity nor pole-local response-compatible stochastic
  transport.
source: root/multiscale-newton-flag/low-child-flag-extension-2026-08-02
audit: >
  Two independent hostile audits rederived the row-image pole-matrix
  convention, exterior-square recurrence, commuting prefix squares,
  four-term top-coordinate formula, all-degree one-letter hostile,
  bounded-lag obstruction, bifiltration maps, and the terminal-Markov and
  LRC-polarization scope boundaries.  Fresh normal and optimized executions
  byte-match the stored transcript and both LF-normalized hashes.
depends_on:
  - THM-3137-finite-stochastic-pole-selector-polytope-and-portability-wall
  - THM-3154-sharp-small-part-profile-jet-and-central-facet-reconstruction
related:
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2420-affine-shell-cross-reference-composition-and-complete-zero-reference-hostile
  - THM-3144-mixed-depth-selector-persistence-death-barcode
  - THM-3149-depth-three-selector-persistence-and-cross-support-wall
  - THM-3155-sharp-depth-four-selector-resurrection-through-degree-eleven
script: 04-computation/gmc_complete_pluecker_pole_holotopy_thm3160.py
output: 05-knowledge/results/gmc_complete_pluecker_pole_holotopy_thm3160.out
script_sha256: b132c2746581af89b30d961038d08f7908696f1ae400a1a83c47de181a0898f0
output_sha256: cc0e1b1578040dd997fa433eb5814c7e3bd8941269b75303dc0637de33081260
hash_basis: LF-normalized bytes
---

# THM-3160 -- complete Pluecker pole holotopy and same-degree selector projection no-go

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The selector staircase has two superficially conflicting features.  The
complete profile jet of THM-3154 observes every partition coordinate and
every Hasse facet, yet adding one more physical pole can resurrect a common
selector after the preceding depth cap has died.  The missing distinction is
between *coordinates* and *transport*.  Pole subtraction acts linearly on
each endpoint, hence on their complete exterior square.  The selector keeps
only one same-degree row of that exterior square.  That row is not closed
under the pole action.

This gives a precise discrete holotopy model.  Prefix paths are flat on the
full Pluecker tensor because pole subtractions commute.  Projecting to the
usual selector current destroys the cross-degree minors needed to move along
an edge.  The result below is universal and algebraic; it does not claim that
the full sidecar is positive or supplies a pole-local, response-compatible
stochastic transition.

## 1. Endpoint space and the pole action

Fix a cutoff `D` and let `Lambda_<=D` be the symmetric functions of degrees at
most `D`, including the constant `1`.  For a scalar `M`, define

```text
T_M f[X]=f[X-M].                                             (1)
```

This is a degree-lowering triangular linear operator.  If `A,B` are two
linear endpoint functionals on `Lambda_<=D`, put

```text
A^M(f)=A(T_M f),                    B^M(f)=B(T_M f).          (2)
```

For the GMC selector, `A` is the signed-bank endpoint `Phi` and `B` is
evaluation on the distinguished residual alphabet `Q`; subtracting a common
physical prefix applies exactly `(2)` to both.

THM-3154 gives the same operator in complete profile coordinates.  With
length marker `u` and small-part markers `w_1,...,w_k`, its one-pole kernel is

```text
K_M^[k](t,u,w)
 =(1-Mt)/[1+(u-1)Mt
            +u(1-Mt) sum_(r=1)^k (w_r-1)(Mt)^r].            (3)
```

Multiplication by `(3)`, truncated through degree `D`, is a unipotent
triangular matrix for `(1)`.  For a prefix multiset `sigma`, the matrix is the
product of the matrices `K_M^[k]`, `M in sigma`.  The factors commute, so the
result depends only on the multiset and not on a chosen deletion order.

## 2. The complete Pluecker tensor is a flat depth object

Define the alternating endpoint tensor

```text
Omega_(A,B)(f,g)=A(f)B(g)-A(g)B(f).                         (4)
```

Equation `(2)` gives the exact pole recurrence

```text
Omega_(A^M,B^M)(f,g)=Omega_(A,B)(T_M f,T_M g).              (5)
```

Fix a coefficient basis `(e_i)` and use the row-image convention

```text
T_M e_i=sum_j (K_M)_(ij)e_j.
```

The endpoint coordinate columns `a_i=A(e_i)` and `b_i=B(e_i)` therefore
transform as `a^M=K_Ma` and `b^M=K_Mb`.  If
`Omega=ab^T-ba^T`, then

```text
Omega^M=K_M Omega K_M^T.                                   (6)
```

Under the opposite function-column convention all displayed matrices are
transposed; the invariant statement is `(5)`.  Thus the full truncated
exterior square is an exact linear sidecar for every prefix edge.
Commutativity of the `T_M` makes `(6)` path-independent around
every square of the prefix-multiset graph.  This is the promised discrete
pole holotopy: the full object has a flat connection.

## 3. The selector is a non-invariant same-degree projection

For a partition `lambda` of `N`, the selector current is

```text
G_N(lambda)
 =A(h_N)B(m_lambda)-A(m_lambda)B(h_N)
 =Omega_(A,B)(h_N,m_lambda).                                (7)
```

The faithful profile code of THM-3154 reconstructs every coordinate in
`(7)`.  It does not retain the minors in `(4)` whose two arguments have
different degrees.

Already the top partition exposes the loss.  Since

```text
T_M h_N=h_N-Mh_(N-1),             T_M m_(N)=m_(N)-M^N,      (8)
```

equations `(5)--(7)` give the four-term recurrence

```text
G_N^M((N))
 =G_N((N))
  -M^N Omega(h_N,1)
  -M Omega(h_(N-1),m_(N))
  +M^(N+1) Omega(h_(N-1),1).                               (9)
```

The last three terms are not same-degree selector coordinates.  More
generally, expanding `T_M m_lambda` in `(5)` uses cross-degree minors of the
two complete endpoint profiles.  Therefore the complete profile jet is a
complete *observer* at a fixed prefix but not an autonomous state variable
for the depth evolution.

## 4. Exact one-letter no-go

The failure is not merely a dimension count.  Work over the rationals and
let `A_x,B_y` be evaluation on the one-letter alphabets `[x]` and `[y]`.
For every `N>=1` and every partition `lambda` of `N`,

```text
G_N^(x,y)(lambda)=0.                                       (10)
```

Indeed, if `lambda=(N)`, then `h_N[x]=m_(N)[x]=x^N` and the
determinant vanishes.  If `lambda` has more than one part, both monomial
evaluations vanish.  Thus *the entire same-degree selector profile in every
degree* is zero for every pair `(x,y)`.

After one common subtraction, however,

```text
h_2[x-M]=x^2-Mx,                    m_(2)[x-M]=x^2-M^2,

G_2^M((2))=M(x-M)(y-M)(x-y).                              (11)
```

The same hostile is not confined to degree two.  For every `N>=2`,

```text
G_N^M((N))
 =x^(N-1)(x-M)(y^N-M^N)
  -(x^N-M^N)y^(N-1)(y-M).                                 (11a)
```

For `(M,x,y)=(1,2,3)`, the values in degrees two through eight are

```text
(-2,-22,-170,-1150,-7322,-45262,-275690),                 (11b)
```

while `(1,2,2)` gives zero in every degree.  Both parents have the identical
input `(10)`, but their children have different selector currents, already at
the GMC starting degree five.  Consequently there is no universal map --
linear, positive, or even set-theoretic -- from the collection of all
same-degree selector currents to the selector current after one pole
subtraction.

The sign also rules out positivity preservation.  The top singleton upset
`{(N)}` has negative mass for each degree in `(11b)`, so the zero parent is
Hasse-positive while every displayed child is not a nonnegative
fine-to-coarse Hasse boundary.

The missing datum is visible in `(9)`: for the first fixture,

```text
Omega(h_1,1)=x-y=-1,                                       (12)
```

whereas it is zero for the second.  Hence at least one cross-degree endpoint
minor is genuinely necessary.  Adjoining all such minors is sufficient by
`(5)--(6)`; no claim of dimension-minimality is made.

## 5. No bounded-lag sidecar works uniformly

The full tensor is finite at a fixed degree cutoff, but no depth-uniform
bounded degree lag suffices universally.  Fix `L>=0`, choose
`N>max(L,1)`, let `A` extract the coefficient of `m_(N)` in the monomial
basis, and let `B` extract the constant coefficient.  Compare this pair with
the zero functional and the same `B`.

Every same-degree selector current is zero for both pairs.  Moreover every
minor `Omega(f,g)` with `f,g` homogeneous and

```text
|deg(f)-deg(g)|<=L
```

is zero for both: the only possible nonzero pairing has degrees `N` and zero.
But

```text
Omega(h_N,1)=1,                    G_N^M((N))=-M^N
```

for every `M!=0`, while the zero control remains zero.  Hence no fixed lag
cap makes the selector projection a universal depth state.  The
degree-`N`/degree-zero minor is genuinely required in this fixture.

## 6. The honest two-parameter selector bifiltration

Fix a finite physical pole multiset and let `S_<=d` be its legal prefixes of
depth at most `d`.  For a lower degree `N_0`, define

```text
C_(D,d)={lambda in Delta(S_<=d):
         G_N(lambda) is Hasse-positive for N_0<=N<=D}.      (13)
```

There are two exact, elementary maps:

```text
C_(D+1,d) subset C_(D,d),                                  (14)

i_d:C_(D,d) -> C_(D,d+1),
    lambda |-> lambda extended by zero on depth d+1.        (15)
```

The maps commute.  Thus horizon is an obstruction direction and depth is a
resource direction; the nonempty cells form a monotone region in the
`(D,d)` grid.  THM-3155 supplies the proved sharp corner

```text
C_(11,3)=empty,       C_(11,4)!=empty,       C_(12,4)=empty. (16)
```

This is a genuine feasibility bifiltration, not a one-dimensional persistence
module generated by pole deletion.  Equations `(9)--(11)` prove why the
tempting stronger interpretation fails: zero extension of a *law* is lawful,
but there is no closed depth transition on its projected selector currents.

There is no obstruction to realizing an arbitrary finite terminal law by an
abstract state-dependent Markov deletion chain: sample the terminal subset,
then reveal a uniformly random order and use posterior transitions.  That
observation is formal and does not transport responses.  The missing object
here is specifically a pole-local, value-dependent, response-compatible
transition whose current evolves by the claimed depth rule.

## 7. Exact verification

Run

```text
python 04-computation/gmc_complete_pluecker_pole_holotopy_thm3160.py
python -O 04-computation/gmc_complete_pluecker_pole_holotopy_thm3160.py
```

and compare byte-for-byte with

```text
05-knowledge/results/gmc_complete_pluecker_pole_holotopy_thm3160.out.
```

The companion uses integer arithmetic only.  It checks the degree-two pole
matrix and its commuting prefix squares, the exterior-square recurrence on a
finite endpoint bank, the four-term top-coordinate identity, all one-letter
parent zero currents through degree eight, the factorization `(11)` on an
integer fixture grid, and the two exact hostile children.

## 8. Scope and surviving sidecar

The theorem is universal endpoint algebra.  Its hostile pairs are valid
one-letter alphabet evaluations, but they are not asserted to be states of
the fixed support-`(1,3)`, bank-`I2` selector polytope.  Therefore the no-go
excludes a universal current-only depth map; it does not rule out an
additional accidental recurrence on one particular finite bank.

The full cross-degree Pluecker tensor is an exact linear and path-independent
repair.  Its pole matrices contain signed virtual-subtraction coefficients,
and `(11)` shows that its projection need not preserve the Hasse cone.  It
does not provide a probability law, pole-local response-compatible stopping
transport, an original-response decomposition, arbitrary-radial NC2, or the
Gaussian Moment Conjecture.

The mechanism is the alternating-polarization analogue of the LRC boundaries
in THM-2380 and THM-2420.  There, separate spectra or self-Grams lose a charged
cross-word or same-shell phase, and a Hermitian cross-correlation restores it.
Here, same-degree determinant rows lose cross-degree endpoint phase, and the
exterior tensor restores it.  This is a related algebraic pattern only: no
LRC observable or endpoint service is constructed from the GMC tensor.

QED.
