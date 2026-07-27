---
id: THM-2545
title: "Word-stratified Hall arrival criterion and owner-word transportation hostile"
status: >
  PROVED + VERIFIED-EXACT.  The theorem is an abstract finite weighted
  transportation result and an exact hostile at the current THM-2461/2537
  semantic interface.  It proves that the recorded owner-word, selected-head,
  and character data, even supplemented by fixed later target-root marginals,
  do not force the positive target hit in THM-2537 equation (56).  It does not
  assert that the hostile is a physical scalar-cover row, construct a later
  target-active field, exclude a row, or prove LRC(14).
source: codex-2026-07-27-word-stratified-arrival-hall
depends_on:
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
related:
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
script: 04-computation/lrc14_word_stratified_hall_arrival_thm2545.py
output: 05-knowledge/results/lrc14_word_stratified_hall_arrival_thm2545.out
script_sha256: 80799fea9575e055682f2df60d3b3ffc9e1724703fe261b0c379ebf87f4fa1ab
output_sha256: c0adf0ef75e458486bb1750f0bd45f209cc66c939b715f012ce3c1da3cf330e4
hash_basis: working-tree bytes (LF)
---

# THM-2545 -- word-stratified Hall arrival criterion and owner-word transportation hostile

THM-2537 reduces its semantic remainder to the positive hitting problem

```text
integral g(z) Psi_tau(e(z))
         A_tar(z,t_tau(e(z))) dz > 0.                         (1)
```

There is an exact finite object behind (1): a word-stratified joint table
between the selected empty-head root and the root of a genuinely later
target-active role.  Separate root marginals do not determine its diagonal.
The exact condition which makes the diagonal unavoidable is a weighted Hall
deficiency after deleting the diagonal from the allowed ancestry graph.

In the maximally permissive case, the condition collapses to one root
overload:

```text
p_sigma(t)+q_sigma(t)>M_sigma.                                (2)
```

The current canon proves neither (2) nor a restrictive support graph.  In
fact a two-by-two aligned/swap pair, independently in every terminal-word
stratum, preserves all one-point data currently connected to (1) and changes
the hit from full mass to zero.  This remains true after tensoring with the
new 42-cut ancestry pairing: its nonzero coefficient contraction lives on
typed-different sheets and is not yet an integral of a positive product on
one actual ancestry.

## 1. The typed joint object

Let `R` be a finite root set and add a cemetery symbol `partial`.  Work on the
selected nonconstant packet as a finite measure space `(Omega,nu)`, where in
the THM-2537 application

```text
d nu(z)=g(z) Psi_tau(e(z)) dz.                                (3)
```

The weight is nonnegative.  Let

```text
sigma:Omega -> Sigma={ {a},{c},{a,c} }
```

be the retained terminal owner word,

```text
h:Omega -> R
```

the selected empty-head root, and

```text
b:Omega -> R union {partial}                                 (4)
```

the root carrying a categorically selected **genuinely later** target-active
role.  Set `b=partial` when that later role is absent.  Thus (4) is not the
empty endpoint of the same predecessor wall.  Constructing (4) from a live
scalar cover is still open; the theorem characterizes what would be needed
once it is constructed.

For every word `sigma`, define

```text
C^sigma_(t,s)
 =nu({z:sigma(z)=sigma, h(z)=t, b(z)=s}),

p^sigma_t=sum_s C^sigma_(t,s),
q^sigma_s=sum_t C^sigma_(t,s),
M_sigma=sum_t p^sigma_t=sum_s q^sigma_s.                     (5)
```

The desired word-stratified target hit is

```text
H_sigma=sum_(t in R) C^sigma_(t,t),
H=sum_sigma H_sigma.                                         (6)
```

If `A_tar(z,t)=1_(b(z)=t)`, then (6) is exactly (1).  Cemetery mass is
genuinely target-inactive and never contributes to (6).

Let

```text
E_sigma subset R x (R union {partial})                       (7)
```

be the compatibility graph allowed by whatever ancestry sidecar is retained,
and assume the actual joint table is supported in `E_sigma`.  Put

```text
D={(t,t):t in R},
E^0_sigma=E_sigma minus D.                                   (8)
```

The deletion in (8), not a Fourier census, is the decisive operation.

## 2. Word-stratified weighted Hall theorem

Fix the margins in (5) and the graphs (7).  Assume at least one joint table
with those margins is supported in every `E_sigma`.  For `S subset R`, write

```text
N^0_sigma(S)
 ={s:(t,s) in E^0_sigma for some t in S}.                    (9)
```

Then the following are equivalent:

1. there is an admissible joint table with `H=0`;
2. for every word `sigma` and every `S subset R`,

```text
p^sigma(S) <= q^sigma(N^0_sigma(S)).                         (10)
```

Consequently the retained data force a positive target hit in every
admissible realization if and only if there are a word `sigma` and a root
set `S` for which

```text
p^sigma(S) > q^sigma(N^0_sigma(S)).                          (11)
```

Moreover every admissible table obeys the quantitative lower bound

```text
H >= sum_sigma Delta_sigma,

Delta_sigma
 =max_(S subset R)
   (p^sigma(S)-q^sigma(N^0_sigma(S)))_+.                     (12)
```

Equation (11) is the exact minimal **margin plus ancestry-support** condition.
For a known joint table, the even more primitive condition is simply
`C^sigma_(t,t)>0` for some `(sigma,t)`: that diagonal cell is precisely the
piece of information lost when one retains only one-point marginals.

### Proof

Fix `sigma`.  Build a bipartite flow network with source capacities
`p^sigma_t`, infinite (or `M_sigma`) capacities on the edges in
`E^0_sigma`, and sink capacities `q^sigma_s`.  A flow of value `M_sigma`
is exactly a zero-diagonal coupling with the prescribed margins.  The
max-flow/min-cut theorem says that such a flow exists exactly when (10)
holds for every left vertex set `S`.  The word strata are disjoint, so a
global zero-hit table exists exactly when this is possible separately in
every stratum.  This proves the equivalence and (11).

For any admissible table and any `S`, split the mass leaving `S` according
to whether it lands in `N^0_sigma(S)`.  The first part is at most
`q^sigma(N^0_sigma(S))`.  Every supported edge landing outside that
neighbourhood must be diagonal, since every off-diagonal supported edge was
put into `E^0_sigma`.  Therefore

```text
H_sigma >= p^sigma(S)-q^sigma(N^0_sigma(S)).                 (13)
```

Maximize (13) over `S` and sum over the disjoint words to obtain (12).
All statements hold for arbitrary nonnegative real weights.  Rational
weights may equivalently be proved after clearing denominators.  Conversely,
every finite rational table is realized by a finite atomic measured space,
so this is an exact semantic model and not merely matrix notation.  QED.

## 3. Complete-off-diagonal specialization and normalization

Suppose `|R|>=2` and every pair in

```text
R x (R union {partial})                                      (14)
```

is compatible.  After deleting the diagonal, the neighbourhood of a
singleton `{t}` is every right vertex except `t`; the neighbourhood of every
set with at least two roots is the entire right side.  Hence (10) reduces
exactly to

```text
p^sigma_t+q^sigma_t <= M_sigma       for every t in R.       (15)
```

More sharply, among all joint tables with the margins in (5),

```text
min H_sigma
 =delta_sigma
 =max_(t in R)(p^sigma_t+q^sigma_t-M_sigma)_+.               (16)
```

At most one root can have a positive overload, since two strict overloads
would use more than the total `2M_sigma` of the two margins.  The lower bound
in (16) is inclusion-exclusion.  If `t_*` has overload `delta_sigma`, place
that mass at `(t_*,t_*)`, subtract it from both corresponding margins, and
apply (15) to the remaining mass.  This constructs an off-diagonal coupling
of everything else and proves equality.

For `M_sigma>0`, put

```text
P^sigma_t=p^sigma_t/M_sigma,
Q^sigma_t=q^sigma_t/M_sigma.                                 (17)
```

The exact normalized minimum is

```text
(min H_sigma)/M_sigma
 =max_(t in R)(P^sigma_t+Q^sigma_t-1)_+.                     (18)
```

For `M_sigma=0`, both sides are defined to be zero before normalization.
Cemetery mass is included in `M_sigma`; it makes (15) harder, not easier, to
violate.  Thus merely proving that a target-active role occurs somewhere is
strictly weaker than proving that it must meet the selected head.

## 4. What THM-2461 and THM-2537 actually retain

THM-2461 proves the prescribed-clock atom-to-word coupling

```text
W_(k_j)(tau,sigma)
 =mu(P_tau intersection E_j intersection T^(-k_j)Q_(j,sigma)).
                                                                    (19)
```

On the complete local bank, (19) has one source-atom row and entries
`rho_(j,sigma)`.  It proves a largest word mass at least `rho_j/3`, but it
does not make the word a function of the source atom.  Its strongest survivor
explicitly asks for a common-base semantic/root intertwiner.

THM-2537 then constructs the selected source and empty-head packets.  The
source packet retains the terminal word and late owner, while its selected
head is disjoint from the old carrier and does **not** inherit that word.  Its
equation (56), reproduced as (1), asks whether the head meets the unique
target-active first-failure role; four of the five roles are target-neutral.

Thus the live canon currently supplies neither:

- the later categorical root map `b` in (4);
- the joint cells `C^sigma_(t,s)` in (5);
- a root compatibility graph together with a Hall-deficient subset; nor
- a word-stratified overload (2).

Even granting the first item, granting word inheritance at the later stage,
and granting target activity everywhere still does not force (1).

## 5. Minimal two-by-two owner-word hostile in `F_13`

Let every word have any prescribed positive mass `m_sigma`.  Use only roots
`0,1` inside `F_13`, split every word stratum into two equal atoms, and compare

```text
                 later root                         later root
               0          1                       0          1

aligned   0  m_sigma/2    0            swap  0    0       m_sigma/2
head      1     0      m_sigma/2        head  1 m_sigma/2     0 .     (20)
```

Both tables have exactly the same data

```text
p^sigma=q^sigma=(m_sigma/2,m_sigma/2,0,...,0),
M_sigma=m_sigma,                                             (21)
```

the same word mass, the same unique THM-2461 source atom, and the same late
owner label.  Put a genuine later target-active role on every atom, so both
models have zero cemetery mass.  Keep its time coordinate, target covector,
threshold layer, and weight `g Psi_tau(e)` identical.  Only the pairing of
its root with the already selected empty head changes.

The complete upstream selected-wall data can be held fixed as well.  For
`tau=1`, use the singleton masks `delta_12` and `delta_0`; their selected
occupied-to-empty walls have heads `0` and `1`.  The pointwise Cayley identity,
selected source/head masses, owner-word labels, old carrier, and any attached
deepest-triangle sidecars are therefore identical in the two models.  The
common two-root marginal has all twelve nontrivial `F_13` Fourier colours,
because

```text
1+zeta^alpha !=0             (alpha=1,...,12)                (22)
```

for a root of unity of odd order thirteen.

Nevertheless

```text
H_aligned=sum_sigma m_sigma,
H_swap=0.                                                     (23)
```

The swap saturates (15) at equality for both roots, so Hall predicts exactly
this freedom.  One root cannot support such a hostile when the target-active
mass is present: its only root-to-root edge is diagonal.  Hence two roots are
minimal.  A three-by-three cyclic permutation is a redundant larger version.

This is an exact hostile to an inference from the **recorded interface data**,
not a claim that either table is already realized by a live scalar-cover row.
It identifies the missing coordinate without inventing a physical example.

## 6. The 42-cut ancestry pairing is a positive control, not the diagonal

The exact artifact

```text
04-computation/lrc14_cut_bundle_ancestry_pairing_opus_20260727.py
05-knowledge/results/lrc14_cut_bundle_ancestry_pairing_opus_20260727.out
                                                                    (24)
```

is `FINITE-EXACT`.  It computes `432/432` nonzero lawful placements of

```text
P_(tau,a)(beta)
 =sum_(alpha in F_13^*) Chat(-alpha)
                         Psi_(tau,a)(alpha,beta),             (25)
```

retains target charge and a positive one-entry owner loop, and verifies both
the root-kappa and cut-torsor collapses are zero.  This is major positive
structure.  Its own type audit, however, says that the temporal-collision
root `u` and static quotient-stalk row `h` are different sheets.  The licensed
map is the character-label alignment `k=alpha`; the septimal cut torsor is not
transported to the stalk.

Accordingly, (25) is a nonzero algebraic contraction in
`Q(zeta_91)`.  It is not yet an identity of the form

```text
P=integral_(one actual ancestry) X(z)Y(z) dnu(z),             (26)
```

much less an integral of a nonnegative product.  Fourier labels may be paired
coefficientwise while their geometric variables remain on different sheets.
The vanishing of both invariant torsor collapses makes this distinction
visible: forgetting either sidecar destroys the pairing.

The hostile (20) can be tensored, word by word, with the entire finite object
in (24).  Keep its cut coefficient, stalk root, owner atom, graft delay,
target covector, and both torsor phases fixed, and attach the later semantic
root by the aligned or swapped table.  The word totals, selected-head and
granted later-target-root marginals, all `432` nonzero coefficient slots,
owner loop, target charge, and both zero torsor collapses remain unchanged;
(23) still changes from full to zero.  The exact companion verifies this
independent-tensor control.

Thus the 42-cut result does not fail.  It lands exactly one categorical level
before (1): coefficient ancestry has been aligned, while event ancestry has
not yet been diagonally supported.

## 7. The extra ancestry condition which defeats every swap

There are three equivalent levels at which a successful sidecar can be stated.

1. **Actual joint support.**  Construct one word/root cell
   `C^sigma_(t,t)>0` on a common physical base.  This directly proves (1).
2. **Support graph plus margins.**  Construct `E_sigma,p^sigma,q^sigma` and
   prove the strict Hall deficiency (11).  This is the exact weakest condition
   visible to those data.
3. **Unrestricted compatibility.**  Prove the root overload (2).  By (16),
   its excess is the exact forced contribution to (1).

A semantic two-cell or base-dependent affine intertwiner of the kind isolated
by THM-2461 and THM-2542 is useful precisely when it supplies level 1 or
restricts `E_sigma` enough to supply level 2.  A character-label equality,
nonzero coefficient product, owner marginal, target covector, or complete root
spectrum alone supplies none of them.

The highest-leverage live test is therefore concrete: on a positive
THM-2537 threshold/word layer, construct the genuinely later selector (4),
compute its root-resolved compatibility graph with the selected head, and
look first for a Hall-deficient subset.  This allows a multi-root ancestry
restriction to win even when no singleton overload exists.

## 8. Exact dependency-free referee

Run

```bash
python3 04-computation/lrc14_word_stratified_hall_arrival_thm2545.py
python3 -O 04-computation/lrc14_word_stratified_hall_arrival_thm2545.py
```

Both executions reproduce

```text
05-knowledge/results/lrc14_word_stratified_hall_arrival_thm2545.out.
```

The referee independently compares weighted Hall inequalities with integer
max flow on every support graph in `189440` three-by-three cases, `1440`
two-by-two cases, and `8896` two-by-three cases with a cemetery column.  It
checks (16)--(18) in `80018` complete-off-diagonal margin pairs and the
cemetery extension in another `12614` pairs, realizes all three owner-word
hostiles, verifies their twelve nontrivial root colours and two-root
minimality, and checks the `42`-cut/`432`-slot typed-sidecar tensor control
together with both exact character collapses.
