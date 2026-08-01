---
id: THM-3044
title: "Pointed correspondence determinant, Hall dual, and cycle boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  A
  supplied physical bijection between selected-head and genuinely later root
  torsors closes each THM-2545 joint table into an endomorphism.  The first
  formal log-determinant coefficient is exactly semantic arrival, and its
  minimum over prescribed margins/support has an exact min-cost transport
  dual.  Hall deficiency is the positivity gate but can be a strict
  quantitative lower bound.  Higher ghosts detect balanced circulation, not
  arrival.  THM-2549 closes conditionally when its semantic later root is
  physically identified with the carry-corrected ancestry root; no such
  pointing is constructed here, so no LRC row is excluded.
source: kind-pasteur-2026-08-01-pointed-arrival-holotopy
audit: >
  Three independent immutable-file audits ACCEPTED the typed pointed closure,
  independent-gauge covariance, formal log-trace identity, exact min-cost
  transportation dual, Hall zero gate and strict quantitative hostiles,
  balanced-cycle/acyclic-support boundary, unpointed orbit no-go, and the
  conditional THM-2549 subpacket closure.  They independently replayed normal,
  optimized, and stored output, matched both LF hashes, and passed the
  documentation checker.  A stronger independent census covered 110,722 gauge
  cells, 18,636 log identities, and 54,234 small transport instances.  The
  audits caught and repaired an overbroad
  shortest-cycle quantifier and two executable-evidence wording defects before
  promotion.
depends_on:
  - THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile
  - THM-2549-future-pullback-target-neutrality-and-cemetery-hall-boundary
related:
  - THM-2829-q11-semantic-reselection-and-fine-ancestry-phase-obstruction
  - THM-2835-q11-semantic-word-horn-and-bockstein-blind-support-no-go
  - THM-2870-prime-power-convolution-versus-physical-diagonal-intertwiner-obstruction
  - THM-2889-dicyclic-reverse-action-joint-carrier-and-skew-lift-separation
  - THM-2894-unmarked-residual-semilattice-order-and-group-clutch-no-go
  - THM-3040-formal-corner-resultant-width-quotient-and-all-order-bernoulli-law
script: 04-computation/lrc14_pointed_correspondence_determinant_hall_dual_thm3044.py
output: 05-knowledge/results/lrc14_pointed_correspondence_determinant_hall_dual_thm3044.out
script_sha256: 659fd435c949e4dc2756225c27f390c85afd970a674f1b8906cafd53d34908e6
output_sha256: 6a4efabb85472c492692ea2da971a30959ba807ab1b8f6f6bf22acb7a674ba21
hash_basis: LF-normalized bytes
---

# THM-3044 -- pointed correspondence determinant, Hall dual, and cycle boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2545 identifies the exact finite object behind semantic arrival: for each
terminal word, a nonnegative joint table from a selected empty-head root to a
genuinely later target-active root.  Its diagonal is the desired mass.  The
two root sets are nevertheless different torsors.  A square matrix written by
giving both torsors the same numerical labels is not yet an endomorphism.

The missing operation is a physical **pointing** between the torsors.  Once it
is supplied, the correspondence closes to an endomorphism, its first
log-determinant ghost is exactly arrival, and every higher ghost has a precise
cycle meaning.  This is the useful part of the proposed ``holotopy'' transfer:
close a directed correspondence before taking a loop invariant.

## 1. Closing the typed correspondence

Fix one THM-2545 word stratum `sigma`.  Let `H_sigma` be its selected-head
root torsor, `B_sigma` its genuinely later root torsor, and let `partial` be
the cemetery symbol.  Both root torsors have the same finite cardinality, but
are not identified.  Write

```text
C^sigma_(h,b) >= 0,
h in H_sigma,        b in B_sigma union {partial},       (1)
```

for a feasible joint table, with margins `p^sigma,q^sigma`, total mass
`M_sigma`, and allowed support graph `E_sigma` as in THM-2545.

Assume in addition a supplied physical root pointing

```text
iota_sigma:H_sigma -> B_sigma                            (2)
```

which is a bijection.  Let `C_R^sigma` be the root-to-root block of `(1)` and
view it as the real linear map `R^B_sigma -> R^H_sigma`.  Let

```text
J_iota:R^H_sigma -> R^B_sigma,
J_iota e_h=e_(iota(h)).                                  (3)
```

The pointed closure is the actual endomorphism

```text
A_sigma=C_R^sigma J_iota:R^H_sigma -> R^H_sigma,
(A_sigma)_(h,h')=C^sigma_(h,iota(h')).                   (4)
```

Define the word-product determinant germ

```text
D_iota(t)=prod_sigma det(I+t A_sigma).                   (5)
```

It has constant term one.  Therefore its formal logarithm is defined, and

```text
[t] log D_iota(t)
 =sum_sigma tr(A_sigma)
 =sum_(sigma,h) C^sigma_(h,iota(h)).                     (6)
```

The right side is exactly the semantic arrival mass relative to the supplied
pointing.  In coordinates in which `(2)` is the identity, `(6)` is
THM-2545's diagonal `H`.

### Independent-gauge covariance

Relabel `H_sigma` and `B_sigma` independently by permutation matrices `P,Q`.
Then the typed objects transform as

```text
C_R -> P C_R Q^(-1),       J_iota -> Q J_iota P^(-1),
A_sigma -> P A_sigma P^(-1).                              (7)
```

Thus `(5)` and every trace in `(6)` are gauge invariant.  The pointing is not
an arbitrary coordinate choice: it is the sidecar which reduces the
independent left/right gauge to a simultaneous gauge.  Merely writing the
same symbols `0,...,12` on both torsors does not supply `(2)`.

## 2. Exact minimum-arrival dual

For `(h,b) in E_sigma`, put

```text
c_iota(h,b)=1 if b=iota(h), and 0 otherwise.             (8)
```

Cemetery edges have cost zero.  Among all nonnegative tables with the fixed
margins and support, let

```text
m_sigma=min sum_((h,b) in E_sigma)
                  c_iota(h,b) C_(h,b).                   (9)
```

Then the exact dual formula is

```text
m_sigma
 =max { sum_h alpha_h p_h + sum_b beta_b q_b :
         alpha_h+beta_b <= c_iota(h,b)
         for every (h,b) in E_sigma }.                  (10)
```

Consequently the minimum forced first ghost over all word strata is

```text
min [t]log D_iota(t)=sum_sigma m_sigma.                  (11)
```

### Proof

Use the finite transportation network with source-to-head capacities `p_h`,
allowed head-to-later edges `(h,b)`, and later-to-sink capacities `q_b`.
Give a pointed edge cost one and every other edge cost zero.  A full flow is
exactly a feasible table and its cost is `(9)`.  The standard finite min-cost
flow dual has one unrestricted potential on each head and later vertex; after
eliminating the source and sink potentials it is precisely `(10)`.  Feasibility
makes the primal polytope nonempty and compact, so both optima are attained and
finite.  Rational margins have rational optimal tables; integer margins have
integer optima by network total unimodularity.  Word strata are disjoint, so
their minima add, proving `(11)`.  QED.

This dual is useful because it permits several Hall cuts to combine in one
certificate.  That combination is invisible if one records only the largest
single cut deficiency.

## 3. Hall is the exact positivity gate, not always the exact amount

Delete the pointed edges from `E_sigma`, obtaining `E^0_sigma`, and let
`N^0_sigma(S)` be the right neighbourhood of `S subset H_sigma`.  THM-2545,
with its diagonal transported through `iota`, gives

```text
m_sigma=0
 iff p_sigma(S) <= q_sigma(N^0_sigma(S)) for every S.    (12)
```

Hence a strict Hall deficiency is exactly equivalent to a forced positive
first ghost.  The quantitative number

```text
Delta_sigma=max_S
 (p_sigma(S)-q_sigma(N^0_sigma(S)))_+                    (13)
```

always satisfies

```text
m_sigma >= Delta_sigma.                                 (14)
```

Indeed the dual potentials `alpha=1_S`, `beta=-1_(N^0(S))` are feasible in
`(10)` and have objective equal to that deficiency.  But `(14)` can be
strict.

### Sharp restricted-support hostile

Take two roots, no cemetery, identity pointing, margins

```text
p=(1,1),       q=(1,1),
E={(0,0),(0,1),(1,1)}.                                  (15)
```

There is one feasible table:

```text
C=[[1,0],[0,1]].                                         (16)
```

Its arrival is `2`, while the largest Hall deficiency is only `1` (from
`S={1}` or `S={0,1}`).  The dual choice

```text
alpha=(0,1),       beta=(1,0)                            (17)
```

has value `2`, proving exact optimality.  Thus the tempting statement
``Hall deficiency equals the minimum ghost'' is false on restricted graphs;
the strongest survivor is the exact dual `(10)` plus the positivity gate
`(12)`.

The strictness is not an artifact of a root lacking an off-diagonal option.
On `Z/5Z`, take unit margins and

```text
E^0={(0,1),(0,2),(1,4),(2,3),(3,0),(4,0)},
E=E^0 union {(h,h):h in Z/5Z}.                            (17a)
```

Every row and column is incident to an off-diagonal edge.  All Hall cuts have
deficiency at most one, but the exact minimum is two.  The permutation table
with row images `(1,4,2,3,0)` gives cost two, while the dual potentials

```text
alpha=(-2,-1,-1,0,0),       beta=(0,2,2,1,1)             (17b)
```

are feasible and also have value two.  This is a genuine multilevel/cascade
dual obstruction missed by every single cut.

When all root and cemetery edges are allowed, THM-2545's special formula is
recovered:

```text
m_sigma=max_h
 (p_h+q_(iota(h))-M_sigma)_+.                            (18)
```

Here a single overload really is exact.  The support restriction in `(15)` is
what makes several dual constraints accumulate.

## 4. Higher ghosts are cycle detectors

For every word define the nonnegative cycle moments

```text
gamma_(sigma,r)=tr(A_sigma^r).                           (19)
```

The complete formal expansion is

```text
log det(I+tA_sigma)
 =sum_(r>=1) (-1)^(r+1) gamma_(sigma,r)t^r/r.            (20)
```

Thus `gamma_1` is arrival.  For `r>=2`, `gamma_r>0` exactly when the pointed
support digraph contains a positive closed walk of length `r`.

Suppose a stratum has no cemetery mass and its pointed margins are balanced:

```text
p_h=q_(iota(h)) for every h.                             (21)
```

Then `(4)` is a finite nonnegative circulation.  If it is nonzero, its support
contains a directed cycle.  Therefore either

```text
gamma_1>0,
```

or there is some `2<=r<=|H_sigma|` with `gamma_r>0`.  This follows by starting
on a positive edge and repeatedly following a positive outgoing edge;
conservation prevents a dead end and finiteness forces a repeated vertex.  A
loop contributes to `gamma_1`; otherwise the first extracted cycle has length
at least two, and its positive edge product occurs in `tr(A^r)`.

More sharply, work in the zero-first-ghost branch `gamma_1=0`.  Nonnegativity
then says that there are no positive loops.  Let `ell` be the length of the
shortest directed cycle.  Then `gamma_r=0` for `1<=r<ell` and

```text
[t^ell]log det(I+tA)
 =(-1)^(ell+1)gamma_ell/ell != 0.                         (21a)
```

Conversely, if the allowed pointed off-diagonal support is acyclic, balance
forces every off-diagonal entry to vanish: a nonzero balanced off-diagonal
flow would contain a cycle.  In that case the first ghost equals the entire
root-to-root mass.  This acyclic-support criterion is a genuine sufficient
arrival theorem, but the current LRC compatibility graph is not known to have
that property.

Balance is load-bearing.  A single directed path has nonzero mass but is
nilpotent, so every `gamma_r` vanishes.  More importantly, a loopless balanced
cycle can have `gamma_1=0` and a nonzero higher ghost.  Higher ghosts certify
circulation, not semantic arrival.

## 5. The no-pointing hostile is an orbit identity

Let `R` have at least two roots and take mass `u>0`.  Relative to an arbitrary
temporary common numbering, compare

```text
C_aligned=uI,          C_pi=uP_pi                         (22)
```

for a fixed-point-free permutation `pi`.  These tables have the same uniform
margins and the same two Gram matrices.  More strongly, they lie in the same
orbit under independent left/right root gauges.  Yet relative to the temporary
identity pointing their arrivals are

```text
u|R|       and       0.                                  (23)
```

The two-root identity/swap pair is minimal.  With three roots and a 3-cycle,

```text
det(I+tuI)=(1+tu)^3,
det(I+tuP_(123))=1+(tu)^3.                               (24)
```

Equation `(24)` does not contradict the orbit statement: computing
`det(I+tC)` has silently used one common basis and therefore smuggled in a
pointing.  Under independent gauges `C` is a correspondence, not an
endomorphism, and `I+tC` is ill typed.  The lawful invariant is `(5)`, in which
`J_iota` transforms with `C`.

This gives the exact stopping rule for the unpointed determinant programme:
all marginals, singular data, and lawful independent-left/right-gauge
invariants may survive while semantic arrival changes.  An ordinary
square-matrix class function already presupposes the missing common basis.
THM-2829 and THM-2835 supply paths and one-way
cospans, not the return map `(3)`; THM-2870 rules out the tempting invertible
same-mask convolution/physical-diagonal intertwiner; THM-2889 and THM-2894
show why conjugate or unmarked endpoint data cannot manufacture a pointing.

## 6. The exact THM-2549 conditional closure

THM-2549 constructs, for a unit role `w`, the carry-corrected old-action
ancestry root

```text
b_w(u,x)=w^(-1)(d_L({wx})-floor(wu)) mod 13
```

and proves pointwise that

```text
b_w(u,x)=h.                                               (25)
```

If a genuinely later semantically target-active role is certified to have
root `b_w` on the selected packet, or a proved physical intertwiner identifies
its root with `b_w`, then `(25)` supplies the pointing and the table is
pointed-diagonal.  On exactly the subpacket where this certification holds,
consequently

```text
[t]log D_iota(t)=sum_sigma M_sigma,                       (26)
```

the whole positive packet mass.  Every Hall cut collapses for the strongest
possible reason: each atom itself is a pointed arrival.

If certification holds only on a positive subpacket, apply `(26)` after
restricting `C,D_iota`, and `M_sigma` to that subpacket; on the unreduced
packet this gives only the corresponding positive lower bound.

The premise of this paragraph is still open in the live LRC proof.  A future
base digit, a numerical root label, a Bockstein path, or a one-way cospan is
not that premise.  Nothing in this theorem constructs the genuinely later
semantic root map, a physical `J_iota`, a target/owner phase, or a scalar-cover
row exclusion.  Under the conditional premise the pointed support is diagonal,
so its all-head Hall cut has deficiency `M_sigma`; this is stronger than mere
positivity.  The LRC ledger is unchanged.

## 7. Why the GMC quotient suggested this operation

THM-3040 succeeds because its width-dependent resultant factors through one
fixed nonzero lower resultant.  Taking a width quotient cancels that fixed
piece; its nonvanishing nevertheless supplies the basepoint at which the
formal logarithm is defined, and summing the relative character leaves only
one integration constant.

The LRC correspondence has the complementary defect.  It already has many
relative paths, transition tables, and loop-like spectral statistics, but its
head and later-root endpoints live in independently gauged torsors.  The map
`J_iota` is the exact analogue of the missing basepoint trivialization: it
closes a path into an endomorphism before the determinant is taken.  This is a
transfer of operation, not a reduction from GMC to LRC.  It explains both the
positive theorem `(6)` and the negative theorem `(22)--(24)`.

## 8. Exact companion and scope

The dependency-free companion uses only integer and rational arithmetic.  It
checks:

- `432` independent-gauge covariance cells and `36` formal log/trace cells;
- every support graph and every margin pair of total mass `1..4` for two heads
  versus two roots plus cemetery: `2667` feasible transportation instances,
  all matching the min-cost dual and Hall zero gate;
- `110` small instances where the largest single Hall deficiency is strictly
  below the exact minimum, including `(15)--(17)`, plus the separate five-root
  robust hostile `(17a)--(17b)`;
- all `1214` nonzero balanced `3x3` tables with entries `0,1,2`; each of the
  `44` loopless tables has a length-two or length-three cycle ghost;
- aligned/permutation double-gauge hostiles, the unbalanced nilpotent-path
  boundary, and the THM-2549 diagonal control.

Reproduce with

```text
python 04-computation/lrc14_pointed_correspondence_determinant_hall_dual_thm3044.py
python -O 04-computation/lrc14_pointed_correspondence_determinant_hall_dual_thm3044.py
```

Both executions equal the stored nine-line transcript byte-for-byte.  The
result is a finite correspondence/transport theorem plus a conditional typing
bridge.  It does not prove that the required physical pointing exists.

**QED.**
