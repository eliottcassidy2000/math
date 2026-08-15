---
id: THM-3391
title: "Weighted common-source cyclic support capacity"
status: >
  PROVED analytic + FINITE-EXACT + INDEPENDENTLY VERIFIED-EXACT.  A
  nonnegative weighted finite source may map to several cyclic quotients.
  Exact pullback-band capacities union-bound on the common source, the exact
  Boolean hypergraph resolves equality, and the correlated cover locus
  retains the endpoint-grid sidecar.  Complete located cosets give a
  multi-high terminal.  The analytic statements in Sections 1--7 are proved;
  the application statuses recorded in Section 5 do not change, and Section
  8 is a prospective selector composition only.  No refined-ledger decrement, physical entry,
  rung, arbitrary-k statement, or LRC(14) conclusion is claimed.
source: root/lrc14-weighted-common-source/2026-08-14
depends_on:
  - THM-3351-projected-k3-z216-sixteen-family-translated-two-high-closure
  - THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph
related:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-3385-odd-fibre-doubling-projection-and-half-even-complement-clocks
  - THM-3388-three-sheet-phase-triangle-cover-clutter
  - THM-3389-four-sheet-typed-cover-clutter
  - THM-3390-complete-z216-projected-wall-by-hybrid-torsion-support-terminal
script: 04-computation/lrc14_weighted_common_source_cyclic_support_capacity_thm3391.py
output: 05-knowledge/results/lrc14_weighted_common_source_cyclic_support_capacity_thm3391.out
script_sha256: 22c2ea187e3d43ca55dd61611a0f6d8a70cf7b1111b1f01cb7338bc1aef7e195
output_sha256: 9cc8b652ae3552f970fae1b8f46f3b6c1d4316a5170d2f9a37eb5e59495e3062
audit_script: 04-computation/lrc14_weighted_common_source_cyclic_support_capacity_independent_audit_thm3391.py
audit_output: 05-knowledge/results/lrc14_weighted_common_source_cyclic_support_capacity_independent_audit_thm3391.out
audit_script_sha256: c95c8463a21890b9ad55e89a3d4ebcbbc69fc3635e5ff8317fc7646686b8d9a7
audit_output_sha256: 6e968c9d2b50f66c2a0b28f7c986923eae95aa90b2f9a225691f419ec6da4b96
semantic_sha256: 8b8ba2d55374ce28698b685bde04f8e16e8cf4dafeb5436b419e0b9fcd642ef7
audit_semantic_sha256: 78e1da67f5662a7f756e221ed18590486600e361cc327bf33187ffb8f4391cca
hash_basis: LF-normalized bytes
---

# THM-3391 -- weighted common-source cyclic support capacity

**PROVED analytic + FINITE-EXACT + INDEPENDENTLY VERIFIED-EXACT.**

## 1. The common-source theorem

Let `C` be a finite ground set with a nonnegative weight `w(c)` and

```text
W=sum_(c in C) w(c)>0.                                  (1)
```

Zero-weight elements may be discarded.  A literal multiset is represented
either by distinct occurrences of weight one or by aggregating coincident
occurrences into an integer weight.

For each blocker `i`, retain a cyclic quotient and a map

```text
phi_i:C -> X_(n_i)=Z/n_iZ.                              (2)
```

At the LRC radius define the pullback band

```text
E_i(theta)={c in C: ||theta+phi_i(c)/n_i||<1/14},       (3)
lambda_i=max_theta sum_(c in E_i(theta)) w(c).           (4)
```

Then for arbitrary, independently chosen phases,

```text
sum_i lambda_i<W                                         (5)
```

leaves one positive-weight source occurrence outside every blocker.  Indeed,

```text
w(union_i E_i(theta_i))
 <=sum_i w(E_i(theta_i))
 <=sum_i lambda_i<W.                                    (6)
```

The survivor in `(6)` is one **common source element** for all blockers.  It
is not a separately optimized survivor in each quotient.

This formulation is stronger and better typed than requiring one common
cyclic modulus.  Arbitrary denominators `n_i` are lawful when the actual
source set `C` and all maps `phi_i` survive.  An `lcm(n_i)` lift is a
sufficient reconstruction only after `C` has been discarded, and then it
must be built from the same literal source occurrences with multiplicity.
Independently optimized quotient supports do not identify a common survivor.

## 2. Exact cyclic-window formula

Specialize to `C=X_n`, let the quotient map be multiplication by an integer
`a`, and retain any nonnegative weight `w` on `X_n`.  Push the weight forward:

```text
b_a(r)=sum_(x: a*x=r mod n) w(x).                        (7)
```

Thus collisions of a nonunit multiplier and pre-existing multiplicities are
both retained.  Put

```text
q_n=ceil(n/7)=floor((n-1)/7)+1.                          (8)
```

Then the exact translation-uniform capacity is

```text
lambda_w(a)
 =max_(t in X_n) sum_(j=0)^(q_n-1) b_a(t+j).            (9)
```

To prove `(9)`, multiply the strict band by `n`.  Its trace is an open cyclic
interval of length `n/7`.  Integer residues in that interval have strict
span

```text
7*(last-first)<n,                                       (10)
```

so they lie in `q_n` consecutive cyclic positions.  Conversely `q_n`
consecutive positions have span

```text
q_n-1=floor((n-1)/7)<n/7,                               (11)
```

and therefore fit inside a suitable open translate.  Nonnegative weights
make the maximum full-window weight exact.

Equivalently, circularly sort the multiset of image residues `a*x mod n`,
double it by appending a copy shifted by `n`, and use a two-pointer window
whose endpoints satisfy `(10)`.  Repeated image residues are repeated in the
sorted list.  This costs `O(|C| log |C|)`; a pushforward bitmap and sliding
sum costs `O(n+|C|)`.

For a distinct located support `S subset X_n`, a high label whose exact
denominator is `n` has a unit numerator.  The exact universal one-high test is

```text
max_(a in (Z/nZ)^*) lambda_(1_S)(a)<|S|.                (12)
```

The maximum over units is load-bearing.  In `X_8`, `S={0,3}` survives the
unit `a=1`, but `a=3` maps it to the adjacent pair `{0,1}`, which one band
contains.

## 3. Boolean exactness and the equality boundary

For blocker `i`, let

```text
H_i={E_i(theta):theta in T}                              (13)
```

be its pullback-band hypergraph on the positive-weight part of `C`.  Under
independent phases, universal survival is equivalent to the absence of a
choice

```text
E_i in H_i with union_i E_i=C.                          (14)
```

Thus `(5)` is a sufficient scalar certificate, while `(14)` is the exact
Boolean certificate.

There is a particularly cheap equality test.  If

```text
sum_i lambda_i=W,                                       (15)
```

then a cover exists if and only if maximum-weight edges can be chosen whose
positive-weight parts partition `C`.  Any cover forces equality throughout

```text
W=w(union_i E_i)<=sum_i w(E_i)<=sum_i lambda_i=W,       (16)
```

so every chosen edge is maximum and every positive-weight source element is
hit exactly once.  The converse is immediate.

For example, in `X_8` the support `{0,1,4,5}` has one-blocker capacity two.
Two unit bands attain the equality sum four and cover the adjacent pairs
`{0,1}` and `{4,5}`.  Strict inequality in `(5)` cannot be weakened.

## 4. Correlated phases and the cover locus

Suppose the phases are not independent but arise from a common base parameter
`y`:

```text
theta_i=theta_i(y).                                     (17)
```

Define the exact projected full-cover locus

```text
B_C={y: C subset union_i E_i(theta_i(y))}.              (18)
```

If every member of `C` is already safe from the fixed low clauses, the
projected safe parameter set is exactly

```text
P=T minus B_C.                                          (19)
```

The independent-phase product hypergraph contains the correlated family
`(17)`.  Hence `(5)` proves `B_C=empty`, but failure or equality of `(5)` says
nothing by itself.  An exact rational event sweep of the correlated edge
tuple may still prove either `B_C=empty` or a sufficiently small measure.

For the THM-2941 projected terminal, a completion would require safe mass
below `36/91`.  Therefore the more general closure target is

```text
mu(B_C)<55/91,                 equivalently mu(P)>36/91. (20)
```

This is the measure-level survivor after Boolean noncover fails.

THM-3387 uses the same object on the full sheet fibre: its `B_q(U)` is `(18)`
for the correlated one-parameter sheet-edge tuple.  Put
`A_C=union_c D_c` and `Gamma_D=D^(-1)Z/Z`.  Pointwise safe-image equality is
equivalent to `B_q(U) subset A_C`.  For an aligned open-cell target, the
corrected criterion is only `B_q(U)\A_C subset Gamma_D`; endpoint-only covers
can survive on the deleted grid (MISTAKE-382), although the literal `<=14`
atlas has none.  The projected terminal instead has no remaining core mask on
`C` and uses pointwise emptiness or the measure threshold `(20)`; finite
endpoint residuals do not affect the latter but do refute an exact `P=T`
claim.

## 5. Exact specialization to the proved theorems

### 5.1 THM-3385: uniform full fibre

Take `C=X_q`, `w=1`, and `phi_u(k)=u*k mod q`.  If

```text
g=gcd(u,q),                     m=q/g,                  (21)
```

then the pushforward bitmap has `m` equally spaced nonzero positions, each
of weight `g`.  Equation `(9)` becomes

```text
lambda(u)=g ceil(m/7),                                  (22)
```

exactly THM-3385's sheet capacity.  Its sum bound is robust against arbitrary
independent phases, even though the physical phases are correlated through
the common base point.

### 5.2 THM-3351: weighted actual-cell source

THM-3351 Section 4 explicitly takes `C` to be the **multiset of actual
complete cells** fixed-safe for the first and sole low drift.  Multiplicity
is retained.  For denominator `d`, its quantity

```text
B_C(d)=max_(u in (Z/dZ)^*) lambda_C(u)                  (23)
```

is exactly `(9)` after pushing the actual cells through

```text
c |-> u*c mod d.                                        (24)
```

Its equation `(7)` is

```text
number of cells safe from every high at y
 >=|C|-sum_i B_C(d_i).                                  (25)
```

The maxima in `(23)` make `(25)` uniform in every allowed unit numerator,
height, and common local coordinate.  Its quantifier order is

```text
for every tuple of high labels, for every y,
there exists one actual cell c safe from all clauses.   (26)
```

The cell may depend on `y`; it is nevertheless one common source cell for
all blockers at that `y`.

The exact twenty-one-case boundary is:

```text
 3  coarse common-modulus capacity certificates,
15  strict weighted-window certificates, weakest surplus 1723,
 3  denominator-two equality cases closed by safe measure 3/7. (27)
```

Thus the weighted common-source theorem is already a proved projected
mechanism, not merely a proposed extension.

### 5.3 Prospective THM-3390 one-high terminal

The current projected complete-wall verifier remains under exact replay and
THM-3390 remains `RESERVED / UNPROVED`.  For one high denominator `d`, it
forms the **deduplicated** support

```text
S={c mod d:c is a fixed-clean actual cell}.             (28)
```

Deduplication is necessary before comparing support cardinality with the
unweighted ambient capacity `ceil(d/7)`: one residue may have arbitrarily
many actual-cell representatives.  The exact weighted predicate
`lambda_C(a)<|C|` is equivalent to `lambda_(1_S)(a)<|S|` for one blocker,
but their numerical margins differ.

The intended terminal's two certificates are corollaries of `(12)`:

1. `|S|>ceil(d/7)` bounds every unit window by the ambient capacity.
2. A located pair of image order `2..7` cannot lie in one strict band under
   any unit numerator.

Exact evaluation of `(12)` strictly dominates both.  It can settle scattered
supports at ambient equality and locate a maximizing hostile unit/window.
No change to the running verifier is needed to state this analytic
strengthening; its existing certificates already imply `(12)` in every case
they close.

### 5.4 THM-3388: the correlated q=3 equality stratum

On the full source `C=X_3`, each unit-speed blocker has capacity one.  Three
blockers therefore reach the scalar equality `sum lambda_i=W=3`.  Section 3
says that an independent-phase cover must partition the three sheets into
three singleton maximum edges, but it does not decide whether those three
owners occur at one physical base phase.

The provisional THM-3388 candidate retains precisely the coordinate erased by
that product-hypergraph quotient.  After assigning speeds `u,v,w` to the
three sheets, it records signed affine gaps

```text
p in P(u,v), q in P(v,w), r in P(w,u)
```

and asks for the circulation law

```text
w*p+u*q+v*r=0.                                         (28a)
```

The hostile `(1,4,7)` has every pair shadow but no closing circulation;
`(1,4,5)` closes.  Thus THM-3391 supplies the equality partition, while the
q=3 program supplies the correlated phase sidecar.  THM-3388 remains related,
not a dependency, until its independent audit promotes it.

### 5.5 THM-3389: the next typed equality stratum

For `X_4`, unit blockers have singleton capacity, whereas a multiplier with
gcd two collapses to an antipodal two-sheet image and has capacity two.  Thus
the nonzero-multiplier scalar equality has only the capacity partitions

```text
1+1+1+1,                 2+1+1,                 2+2.      (28b)
```

The Boolean equality criterion sharpens `(28b)`: four singleton edges must
use all four sheets; a pair edge and two singleton edges must use its
complementary sheets; two pair edges must choose the two complementary
antipodal fibres.  Correlated physical phases may still delete some of these
product-hypergraph covers.
This is the exact transverse type signature reserved by THM-3389's four-sheet
clutter stub.  The abstract `X_4` hypergraph additionally permits the zero
multiplier, whose one maximum edge is all four sheets, but that trivial
capacity-four species is not transverse and is not assigned to THM-3389.
No correlated q=4 classification is claimed here.

### 5.6 Prime-sheet equality beyond triangles

For a prime `q<=7`, every nonzero multiplier is a unit and has capacity one.
Consequently `m<q` blockers leave a sheet by the strict sum theorem, while
`m=q` reaches equality and can cover only by assigning one maximum singleton
edge to every sheet.  The independent-phase equality skeleton is therefore a
`q`-uniform clutter on blocker speeds.  At `q=2` it is the THM-3387 graph; at
`q=3` it is the phase-triangle program of THM-3388.  At `q=5,7`, the next
honest object is a correlated phase-cycle clutter, not a graph of pair tests.
THM-3391 proves this scalar/Boolean reduction but does not classify those
higher correlated cycles.

## 6. Located cosets and a multi-high terminal

Suppose a distinct support contains a complete coset `T` of a cyclic subgroup
of order `r`.  Restricting multiplier `a` to `T` gives the exact uniform
capacity

```text
kappa_T(a)
 =gcd(a,r) ceil((r/gcd(a,r))/7).                        (29)
```

For a native exact denominator `d`, every allowed numerator is a unit modulo
`d`, hence also modulo every `r|d`.  Equation `(29)` simplifies to

```text
kappa_T(a)=ceil(r/7).                                   (30)
```

Consequently `m` arbitrary high blockers with that native denominator leave
a point of `T` whenever

```text
m ceil(r/7)<r.                                          (31)
```

For `2<=r<=7`, `(31)` is exactly `m<r`.  An order-three full coset defeats
two highs; an order-four coset defeats three; and an order-seven coset defeats
six.  This can bypass a duplicate-two-high scalar gap in a future projected
wall.

A pair alone does not give `(31)`: different blockers may kill its two
endpoints.  A full located coset or the exact common-source hypergraph is
required.  For mixed native denominators, retain the common source maps and
use `(4)--(6)`; do not silently treat the lifted multipliers as units.

## 7. The denominator-two equality boundary

In each of THM-3351's three `d_1=d_2=2` equality cases, `(25)` is zero.  Since

```text
B_C(2)=max(weight of even cells, weight of odd cells),  (32)
```

equality forces the two parity fibres of the actual-cell multiset to have
equal weight.

Write a concrete high drift as

```text
z_i=(L/2)p_i,                     p_i=2h_i+1 odd.       (33)
```

On parity sheet `k`, its local phase is

```text
k/2+p_i*y/2 = p_i*(y+k)/2 mod 1.                       (34)
```

For fixed heights this is exactly THM-3387's `q=2` sheet graph.  The two
highs cover both parity sheets somewhere if and only if

```text
p_1+p_2>7 gcd(p_1,p_2).                                (35)
```

Thus graph nonedges give `B_C=empty`.  Edges can make `B_C` nonempty, but
they need not make it large.  The denominator-only THM-3351 state had erased
the odd numerators and had to quantify over both genera: `(p_1,p_2)=(1,3)`
is a nonedge, while `(1,9)` is an edge whose exact cover-locus measure is
`4/63`.

THM-3351's measure proof is the uniform survivor.  For any one fixed cell,
each odd high is dangerous on at most `2/7` of the local-coordinate circle.
If all cells are covered, that fixed cell is covered, so

```text
mu(B_C)<=4/7=52/91<55/91,
mu(P)>=3/7=39/91>36/91.                                (36)
```

This closes both graph edges and nonedges and is strictly stronger than
requiring `B_C=empty`.

## 8. Packetwise selector quantifiers

Let `P` be a packet and `J(P)` a set of lawful selector addresses.  Suppose
address `j` retains a common weighted source `C_(P,j)` and all blocker maps.
For numerators fixed by the packet, the lawful certificate is

```text
for every P, there exists j in J(P) such that
sum_i lambda_(C_(P,j),phi_i)(a_i)<W_(P,j).              (37)
```

If numerator directions remain universally quantified, replace each term by
the maximum over its allowed numerator set.  The conclusion has order

```text
for every P, there exists j, for every phase tuple,
there exists c in C_(P,j) safe from every blocker.      (38)
```

One may not use a different selector `j` or a different source survivor for
each blocker.  The final survivor may depend on the common phase tuple.

This is presently only a meta-connection to the coherent-shift multicell
rescue.  That result verifies exact literal interval unions at individual
body-safe cells, but it has not constructed packet-indexed cyclic supports
and common quotient maps of the form required by `(37)`.  A literal lift is
needed before theorem-level composition.

## 9. Sharp hostiles and quotient losses

1. **Ambient equality goes both ways.**  In `X_14`, `{2,3}` is contained in
   one band, whereas the antipodal `{0,7}` is not.
2. **The maximum over units is necessary.**  `X_8,{0,3}` survives unit one
   but is contained after multiplication by unit three.
3. **Nonunits collapse torsion.**  In `X_4`, multiplier two sends `{0,2}` to
   one image point.
4. **Multiplicity cannot feed an unweighted bound.**  One hundred source
   cells at residue zero still have band load one hundred, not at most
   `ceil(14/7)=2`.
5. **Sum equality can cover.**  Two bands partition `{0,1,4,5}` in `X_8`.
6. **Correlation can rescue equality.**  The `q=2` pairs `(1,3)` and `(1,9)`
   have the same independent capacity sum, but only the latter has nonempty
   correlated cover locus.
7. **Separate quotient survivors are invalid.**  Two blockers may leave
   different source elements while their union covers all of `C`.
8. **Strict openness is load-bearing.**  An order-seven pair is safe at
   exact separation `1/7`; closed bands would change the boundary.
9. **A removed grid is a real sidecar.**  THM-3387's corrected witness has
   `B_2(U)` outside the core only at `{3/14,11/14}` on the deleted grid: the
   open-cell completion holds while pointwise safe-image equality fails.

## 10. Exact and independent verification

The canonical companion checks:

- `22,230` exact lambda instances: every subset of size at most three in
  `X_n`, every multiplier including zero and nonunits, for `2<=n<=14`;
- sort/double/two-pointer equality with the pushforward bitmap/window formula;
- `819` full-fibre gcd capacities through modulus forty;
- `4,574` located image-order pair controls;
- all `400` ordered odd-speed pairs through `39` against the `q=2` gcd graph
  and exact rational cover-locus measure;
- all `54` order-`2..7` coset/blocker threshold controls; and
- weighted multiplicity, unit-max, nonunit, equality-partition, and lost-
  source-identity hostiles.

Reproduce with

```text
python3 04-computation/lrc14_weighted_common_source_cyclic_support_capacity_thm3391.py
python3 -O 04-computation/lrc14_weighted_common_source_cyclic_support_capacity_thm3391.py
python3 04-computation/lrc14_weighted_common_source_cyclic_support_capacity_independent_audit_thm3391.py
python3 -O 04-computation/lrc14_weighted_common_source_cyclic_support_capacity_independent_audit_thm3391.py
```

Each ordinary/optimized pair byte-matches its stored output.  The independent
implementation does not import the primary probe: it sweeps exact phase
arrangements, uses a common rational grid for the q=2 cover locus, and
reproduces the deleted-grid endpoint residual `{3/14,11/14}` with zero
residual open cells.  LF-normalized hashes are frozen in the frontmatter.
Its exact audit counts are

```text
strict windows  18,677       equality cases   5,724
mixed maps         216       unit supports   32,751
cosets            6,238      q=2 pairs          400
```

## 11. Scope and further transfers

The analytic statements in Sections 1--7 are the proved theorem; the cited
application programs retain their own statuses.  Sections 9--10 freeze its
hostile boundaries and two exact audit paths.  Section 8 is a lawful selector
statement but remains prospective: the coherent-shift multicell rescue has
not yet been lifted to one packet-indexed common cyclic source.

The meta-pattern **Keep the common source under cyclic blocker quotients** is
now routed in `00-navigation/META-PATTERNS.md`.  Its trigger is several
cyclic blockers acting on one weighted source; its counterindications are
lost source identity, erased numerator gcd type, changed endpoint convention,
or a deleted target grid.  The next concrete extensions are the q=3
circulation sidecar, the q=4 typed equality clutter, and exact
`lambda_S(a)` replacement of the two coarse one-high terminals in THM-3390.

This theorem changes no refined ledger.  It proves no physical entry, rung,
arbitrary projected sector, arbitrary-residue sector, or LRC(14).
