---
id: THM-2850
title: "Paid distinguished-coloop round budget and slope-grammar rank defect"
status: >
  PROVED + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the paid
  labelled seed/edge-row forcing semantics, the coordinate change
  (a,b)->(a+b,a-b) identifies forcing with distinguished-coloop peeling in
  A=[B|D_rho B].  The round headroom is cycle surplus minus
  slope-compatible bicirculations; q nonempty rounds force score at least
  1+1/q and two rounds force at least 3/2.  With independently variable
  row slopes, generic full paid-row rank is equivalent to partition into
  two B-independent sets (two forests after adjoining ground); under
  label-sharing constraints this condition is necessary but need not be
  sufficient.  The theorem does not certify the private Arithmetic Kakeya
  verifier or compile mode-three quotient witnesses into the intuitive
  same-H recursion.
source: root/ak-vandermonde-round-rank-2026-07-28
depends_on: []
related: []
script: 04-computation/ak_distinguished_coloop_round_rank_audit.py
output: 05-knowledge/results/ak_distinguished_coloop_round_rank_audit.out
script_sha256: b9212c6f27d839b591ab675149f3c2168ceb51d10cd9d828a670eeb26ad3b90c
output_sha256: ef4ba726489818f19d85d1b74bf64b27893025524d8175590fadbc10f215daae
secondary_script: 04-computation/ak_projective_matroid_audit.py
secondary_output: 05-knowledge/results/ak_projective_matroid_audit.out
secondary_script_sha256: f5159f97717a09dcd49fccd2a28e39a1c535c46309f23e0dbf483eed93fb7f05
secondary_output_sha256: 1616fc1fd9d8101429ac1a00389289df915861477f0f5f40061bf6e2ac1f0605
independent_script: 04-computation/ak_distinguished_coloop_round_rank_independent_audit.py
independent_output: 05-knowledge/results/ak_distinguished_coloop_round_rank_independent_audit.out
independent_script_sha256: 410525b4e69d2785f35d029daa97de2ca5f6b817ddd74e58e596630470fadbc5
independent_output_sha256: eed5875acd62998e467052916012a24738c0af9a8197fb37b2fdd6941b2934e1
hash_basis: LF-normalized bytes
---

# THM-2850 -- paid distinguished-coloop round budget and slope-grammar rank defect

**PROVED + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This is a theorem about the paid labelled-graph semantics described below.
It is not a theorem about the private Arithmetic Kakeya verifier, whose code
is unavailable, nor a constructibility proof for the mode-three records.

## Statement

Let a forcing certificate have `u=n-t>0` initially live vertices and `g`
charged generator slots, counted with multiplicity.  All ranks, kernels, and
spans below are over `Q`.  Assume:

1. every nonzero generator is a labelled seed or a labelled physical edge,
   and the active nonzero generator-row multiset has cardinality
   `e_1<=g`;
2. every nonzero label `(a,b)` has `a+b != 0`;
3. forcing succeeds in `q>=1` simultaneous nonempty rounds of sizes
   `f_1,...,f_q>0`, with `sum_j f_j=u`.

This includes the strict/equal-suffix game, the per-suffix physical-edge
game, and any quotient multigraph obtained by identifying vertices, provided
that every pushed-forward row comes from a charged slot.  Loops may become
zero and disappear, parallel rows retain their charged multiplicity, and
always `e_j<=e_1<=g`.

For the strict AK accounting, every nonzero `f_i(prefix)` is charged
`h=d_(i+1)...d_k` in `m` and supplies exactly `h` equal-suffix physical edge
rows; every seed is charged once and supplies one row.  Mode three charges
each physical edge occurrence, while quotienting can only identify or zero
rows.  Hence `e_1<=g` and `r_1<=e_1<=g`; tied or duplicate rows can only
lower rank.  The suffix-unconstrained loose rule is excluded because one
cost-`h` entry can supply a complete bipartite row bank of span dimension as
large as `2h-1`.

At the start of round `j`, let:

- `u_j` be the number of live vertices;
- `A_j` be the generator matrix restricted to the live columns, in coordinates
  `(s,d)=(a+b,a-b)`;
- `r_j=rank(A_j)`, `H_j=r_j-u_j`, and `nu_j=2u_j-r_j`;
- `p_j` be the additional rank lost by deleting the `s`-columns of the fired
  vertices after their `d`-columns have been deleted.

Then

```text
rank(B_j)=u_j,
0 <= p_j <= f_j,
r_{j+1}=r_j-f_j-p_j,
H_{j+1}=H_j-p_j,
nu_{j+1}=nu_j-f_j+p_j,
f_j <= H_j <= H_1 <= g-u.
```

Here `B_j` is the signed incidence block in
`A_j=[B_j | D_rho B_j]`.  Consequently,

```text
r_1 = u + sum_j p_j,
H_1 = sum_j p_j,
max_j f_j <= H_1,
g >= u + max_j f_j,
g/u >= 1 + 1/q.
```

There is also an exact cycle-space interpretation.  Let `e_j` be the number
of generator rows that remain nonzero after restriction, let
`C_j=ker(B_j^T)`, let `c_j=e_j-u_j=dim(C_j)`, and let `z_j` be the number
of bridges (rows absent from the support of every circulation).  Then

```text
sigma_j := e_j-r_j
         = dim { y in C_j : D_rho y in C_j },
H_j = c_j-sigma_j,
sigma_j >= max(0, 2c_j-(e_j-z_j)),
f_j <= H_j <= min(c_j, u_j-z_j).
```

Thus `H_j` is the ordinary cycle surplus after subtracting the
**slope-compatible bicirculation space**.  A large bridge forest or a
projectively invariant circulation is not usable forcing headroom.

The initial full-row-rank question separates topology from
**independently variable row slopes** exactly.  For a fixed incidence matrix
`B` with row set `E`, the polynomial matrix

```text
A(rho) = [B | diag(rho_e) B]
```

has full row rank for some (equivalently, Zariski-generic) independently
chosen rational row slopes
if and only if `E` can be partitioned as `E=E_0 disjoint_union E_1` so that
the rows `B[E_0]` and `B[E_1]` are each independent.  In the augmented graph
where seeds are edges to ground, these are two forests.  In particular every
row subset `S` must satisfy

```text
|S| <= 2 rank(B[S]).
```

If this local condition fails, no choice of labels can recover the lost
rank.  If a two-forest partition exists, only a proper algebraic cancellation
locus in the full row-slope space is bad.  A certificate grammar may constrain
many rows to share one slope.  On that parameter subspace the determinant can
vanish identically—for example a triangle with one common slope has
`[B|rho B]` of rank two—so the two-forest condition remains necessary but is
not sufficient without checking the restricted determinant.

More generally, in the unrestricted row-slope space,

```text
max_rho rank A(rho)
  = min_(S subset E) (|E\S|+2 rank(B[S]))
  = |E|-delta(B),
delta(B):=max_(S subset E)(|S|-2 rank(B[S])).
```

Thus the unrestricted generic defect is always exactly `delta`; full row
rank occurs exactly when `delta=0`, equivalently when the whole row set
partitions into two `B`-independent sets.

Precisely, let `L` be a rational affine-linear subspace of row-slope space
with Zariski-dense rational points; coordinate-equality constraints are the
main label-sharing case.  There is a full-row-rank point in `L` if and only if
some row-full minor polynomial restricts nontrivially to `L`.  When it exists,
full rank holds on a nonempty Zariski-open dense subset of `L`.

More quantitatively, define the topology-only defect

```text
delta_j = max_(S subset E_j) (|S|-2 rank(B_j[S])).
```

Then `sigma_j>=delta_j`: indeed
`rank(A_j)<=|E_j\\S|+2rank(B_j[S])`.  For a legal slope space `L_j`, put

```text
sigma_L_j^gen=e_j-max_(rho in L_j) rank A_j(rho),
gamma_L_j=sigma_L_j^gen-delta_j,
kappa_j(rho)=sigma_j(rho)-sigma_L_j^gen.
```

Then all three defects are nonnegative and

```text
sigma_j(rho)=delta_j+gamma_L_j+kappa_j(rho),
H_j=c_j-delta_j-gamma_L_j-kappa_j(rho).
```

Here `delta` is forced by topology/density, `gamma` is forced by the
label-sharing grammar, and `kappa` is special-coefficient cancellation
inside that grammar.  Thus the older undifferentiated residual
`epsilon_j=sigma_j-delta_j` splits exactly as
`gamma_L_j+kappa_j(rho)`.  In the unrestricted row-slope space `gamma=0`
for every `B`; the two-forest condition is the sharper condition
`delta=0`.

Equivalently, without naming the legal slope space,

```text
H_j = c_j-delta_j-epsilon_j,
epsilon_j=gamma_L_j+kappa_j(rho).
```

For exactly two rounds, writing `f_1+f_2=u`,

```text
p_2=f_2,
r_1=u+f_2+p_1,
g/u >= 3/2.
```

If additionally `g/u <= 5/3`, both round sizes obey the necessary balance
window

```text
2u-g <= f_1,f_2 <= g-u.
```

More sharply, if `lambda=g-r_1` is unused paid-cost/rank slack, both endpoints
sharpen:

```text
2u-g+lambda <= f_1,f_2 <= g-u-lambda.
```

Thus the smallest exact `5/3` search shapes with zero slack have only the
following two-round windows:

```text
(g,u)=(5,3): each f_j is in [1,2];
(g,u)=(10,6): each f_j is in [2,4];
(g,u)=(15,9): each f_j is in [3,6];
(g,u)=(20,12): each f_j is in [4,8].
```

## Proof

Rule (2) depends only on the rational row span: a rational witness clears
denominators to an integral combination, while every integral combination is
rational.  After the vertexwise change `(a,b)->(s,d)`—invertible over `Q`
with determinant `-2`—scale each physical row by its nonzero `s=a+b`.  A
seed is an edge to ground, and a labelled edge with projective slope
`rho=d/s` has row

```text
(ordinary signed incidence | rho * ordinary signed incidence).
```

Hence `A_j=[B_j | D_rho B_j]`.

Rule (2) fires `v` exactly when the row space contains a nonzero singleton
on `d_v`.  By annihilator duality, this is equivalent to every right-kernel
vector having `d_v`-coordinate zero, which is equivalent to the `d_v` column
not lying in the span of the other live columns.  Thus firing is exactly the
condition that `d_v` is a coloop of the represented column matroid.

The cost hypothesis matches the strict grammar row by row.  A nonzero
`f_i(prefix)` has suffix multiplicity `h=d_(i+1)...d_k`; its contribution to
`m` is `h` and it creates exactly those `h` physical rows.  Each seed
contributes one to `r` and one row.  In the per-suffix/mode-three grammar,
each edge occurrence is separately charged.  Pushing rows through a quotient
can turn them into zero rows or duplicates but cannot create a new uncharged
row.  Therefore the number of active nonzero rows satisfies `e_j<=e_1<=g`,
and in particular `r_1<=g`.

Every live connected component must be pinned by a seed or an edge to an
already deleted vertex.  Otherwise the right-kernel vector which is zero in
all `s` coordinates and constant `1` in the component's `d` coordinates
annihilates every row.  It gives a column dependence using every `d_v` in
that component, so none of those distinguished columns is a coloop.  The
component can never begin to fire, contradicting success.  A signed incidence
matrix with every component pinned has full column rank, proving
`rank(B_j)=u_j`.

Each fired `d_v` is a column-matroid coloop of `A_j`.  Deleting the `f_j`
fired `d`-columns therefore lowers rank by exactly `f_j`.  The remaining
matrix still contains the full-rank `B_j`, so
`r_j-f_j >= u_j`, i.e. `f_j<=H_j`.  Deleting the matching `f_j` `s`-columns
loses some additional `p_j` with `0<=p_j<=f_j`, giving the three transition
identities.

For the cycle refinement, the left kernel of `A_j` consists exactly of
those circulations `y in C_j` for which the slope-weighted vector
`D_rho y` is another circulation.  This proves the formulas for `sigma_j`
and `H_j`.  Every circulation vanishes on each bridge, so both `C_j` and
`D_rho C_j` lie in an ambient coordinate space of dimension `e_j-z_j`.
The map

```text
C_j -> (Q^(e_j-z_j))/C_j,    y |-> D_rho y mod C_j
```

has kernel the bicirculation space and target dimension
`(e_j-z_j)-c_j`.  Rank-nullity gives
`sigma_j >= 2c_j-(e_j-z_j)`, and hence
`H_j<=min(c_j,u_j-z_j)`.

For the topology/slope separation, suppose first that the rows of `A(rho)`
are independent.  Some square row-full minor is nonzero.  Expand its
determinant according to which rows use columns from the second block.
Every nonzero term partitions the rows into two sets on which the
corresponding `B` minors are nonzero, giving the two independent sets.
Conversely, choose nonsingular `B` minors for a two-independent-set
partition.  In the corresponding square minor of `A(rho)`, that partition
contributes the squarefree monomial `product_(e in E_1) rho_e` with nonzero
coefficient.  Distinct row partitions give distinct monomials, so this
determinant polynomial is not identically zero.  Rational points avoid its
zero set, and every rational `rho=p/q` lifts to the legal integral label
`(q+p,q-p)`.

The same determinant expansion applied to every row subset identifies the
generic row matroid with the union of two copies of the row matroid of `B`.
The matroid-union rank formula gives

```text
max_rho rank A(rho)
 =min_(S subset E)(|E\S|+2 rank(B[S]))
 =|E|-max_S(|S|-2 rank(B[S]))
 =|E|-delta.
```

This proves `sigma^gen=delta` in the unrestricted slope space, including
the rank-deficient case `delta>0`.

For a legal rational affine-linear slope space `L`, restrict every row-full
minor polynomial to `L`.  Full row rank at one point is equivalent to one
restricted minor being nonzero there, hence to that restricted polynomial
not being identically zero.  Its nonvanishing locus is Zariski open and
nonempty; density of `L(Q)` supplies rational points.  Conversely, if every
restricted minor vanishes identically, no point of `L` has full row rank.
The maximum rank on `L` is attained on a nonempty Zariski-open set, which
proves the `sigma_L^gen`, `gamma_L`, and `kappa` decomposition above.
Coordinate ties can identify distinct forest monomials after restriction,
so their coefficients may cancel identically even though the unrestricted
minor is nonzero.

At termination both rank and live count are zero.  Summing the transitions
gives `r_1=sum_j(f_j+p_j)=u+sum_j p_j`; also `r_1<=g`.  The global bounds
follow.  In the last of two rounds all `2f_2` live columns are independent:
the `f_2` distinguished columns are coloops and the surviving incidence block
has rank `f_2`.  Thus `r_2=2f_2` and `p_2=f_2`.  The remaining two-round
claims follow from `f_1,f_2<=H_1<=g-u` and `f_1+f_2=u`.
If `lambda=g-r_1`, then `H_1=r_1-u=g-u-lambda`.  Hence
`f_i<=g-u-lambda`, while the complementary-round identity gives
`f_i=u-f_(3-i)>=2u-g+lambda`, proving the sharpened window.

## Exact witness replay

Reproduce with:

```bash
python3 04-computation/ak_distinguished_coloop_round_rank_audit.py
python3 -O 04-computation/ak_distinguished_coloop_round_rank_audit.py
```

Normal and optimized output are byte-identical.  Exact profiles:

```text
strict 13/7:     r0=13, H0=6, f=(2,5),   p=(1,5)
per-suffix 7/4:  r0=14, H0=6, f=(1,4,3), p=(0,3,3)
quotient 12/7:   r0=12, H0=5, f=(1,3,3), p=(0,2,3)
quotient 9/5:    r0=9,  H0=4, f=(2,3),   p=(1,3)
```

The script's zero-based names `r0,H0` are the theorem's initial
`r_1,H_1`.

All four have zero paid-row/rank slack, zero initial bicirculation dimension,
and zero bridges in every round.  Later rounds acquire bicirculations:
their active `sigma_j` profiles are respectively `(0,1)`, `(0,1,1)`,
`(0,1,1)`, and `(0,2)`.  Exact dense-subset witnesses show
`sigma_j=delta_j` in every round, so both `gamma_L_j=0` and
`kappa_j(rho)=0` throughout: none of these records loses rank from its
label-sharing grammar or from accidental coefficient cancellation beyond
the topology-forced two-arboricity defect.  The minimal
`(subset size)/(B-rank)` cores after the initial round are:

```text
strict 13/7:     7/3
per-suffix 7/4:  7/3, 7/3
quotient 12/7:   5/2, 7/3
quotient 9/5:    6/2   (defect 2)
```

The literal suffix-unconstrained rule is outside the theorem: one charged
`f_i` entry supplies a complete-bipartite bank of generator rows while `m`
charges only a matching-sized layer.  This is exactly the hypothesis that
fails in the score-to-one exploit family.

The quotient extension is algebraically sound but deliberately limited:
it audits the pushed-forward labelled multigraph after identifications.  It
does not establish recursive same-`H` legality, transport prefix/suffix label
dependence, or decide whether the benchmark counts duplicate/loop slots in
the same way.  Those remain the format-compilation obstruction.
