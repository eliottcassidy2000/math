# Selected-root activation is a boundary port birth, not a rearm

**Status: PROVED integration of accepted incoming formulas and our marked-module
identities; conditional completion criteria remain conditional. No new census.**

Incoming source: commit `d0383cbc2486b985007db238677aab6c827fdb02`,
`05-knowledge/results/seed_recursion_births_20261111/PROOFS.md`, R10–R11,
and its `selected_root_bridge.json` / `residual_grids.json`. These were read
via `git show origin/main:path`; origin/main had advanced to
`46b06e8a864c048cc62bc640447b9dfc23b79ba6` at inspection. R11 explicitly
credits Session260's constructor ladders and independently places them at two
previously classified deletion roots. Its new useful contribution here is
that activation/deletion-root interpretation, not discovery of those ladders.

## 1. A literal marked isomorphism, not only matching numbers

In B=(24,1,3,6,14), vertices2→1 have identical exterior neighborhoods.
Vertices4→3 are another such consecutive twin pair. Coarsen either pair into
one transitive module. This permits redistributing how many clones are called
“the first original vertex” versus “the second original vertex” without
changing the resulting marked tournament.

The two four-vertex quotients, listing the three exterior vertices followed
by the module port, may be written

    Q_I=(6,8,10,1),   Q_II=(8,1,3,6),  port=3.

Inflate the port as

    P_u → x(selected) → G_g → y(selected) → S_v,
    u,g,v≥0.

The exterior vertex orders are (0,3,4) for Q_I and (0,1,2) for Q_II.
At u=g=v=0 this recovers the original marked B with selected sets {1,2}
and {3,4}, respectively, up to the displayed module ordering.

The exact unified formulas are

    H=12·2^(u+g+v)−1,
    D=0^g+2^g−1−6·2^(u+v)(3^g−2^(g+1)+1),
    q=U=C·2^u3^v(2^(g+1)+3^(g+1))−e·0^u,             (1)

where (C,e)=(6,1) for Q_I and (10,3) for Q_II, and0^0=1.

No new coefficient computation is needed. For u≥1, use coverage/PROOF.md's
v=1 or3 family with a=u−1,b=v,c=g. For u=0, use the all-position theorem
in constructors/RESEARCH.md §3 with j=g,k=v. The two pieces are precisely
the same coarsened module described above. The formulas agree on their
actual common graph interpretations. In particular D=1 iff g=0 or1.

R11 is exactly the subfamily u=g=0,v=m−1. It gives

    q_I=30·3^(m−1)−1=10·3^m−1,
    q_II=50·3^(m−1)−3.

For m≥2 these are our former v=2 and4 boundary presentations with
a=c=0,b=m−2. This correspondence is a marked graph isomorphism: the joint
transitive block has its first two vertices selected and all remaining
vertices unselected. The m=1 boundary supplies q29 and47.

## 2. Proved corollary: two roots, commuting directions, and one legal gap birth

Put δ=e·0^u. The previous L/R laws remain

    L:(q,δ)→(2q+2δ,0),       u→u+1,
    R:(q,δ)→(3q+2δ,δ),       v→v+1.

At g=0 there is also a legal one-time birth G between the selected pair:

    G:(q,δ,g=0)→((13q+8δ)/5,δ,g=1).                    (2)

The division is integral because q+δ=5C·2^u3^v before G. The map is
physical insertion of one unselected clone between x and y. It commutes
with L and R as an ordered marked graph operation. Explicitly, the two
orders of L/R give (6q+6δ,0); the two orders of L/G give
((26q+26δ)/5,0,g=1); and the two orders of R/G give
((39q+34δ)/5,δ,g=1). The latter two comparisons begin at g=0.
Thus this is both a marked-graph and an exact response-state commutation
statement, with the one-use G guard retained. A second gap birth is illegal for unit seeds: g=2 gives

    D=3−12·2^(u+v)<0.

Thus our previous eight starting presentations are just first-generation
L/R descendants of four boundary states, and those four are themselves
obtained from two roots by optionally applying G:

    q29 --G--> q77,       q47 --G--> q127.

For example q29 has R-child89 and L-child60; q47 has R-child147 and
L-child100. These are the initial states in our former four-vertex-indexed
families. All states in this unified carrier admit a literal deletion path
to one of the two roots, decreasing u+g+v. This is a complete descent **inside
this specified carrier**, not a descent of arbitrary unit seeds.

R11's first R birth creates a fresh unselected suffix vertex whose further
inflation is stable. The corresponding selected-root activation is real and
escapes the old-vertex/unselected-port hypothesis of incoming R2. The gap
birth G illustrates a different outcome: its new unselected vertex cannot
be inflated again while preserving D=1. “New vertex” does not automatically
mean “new stable port.”

## 3. Why this still supplies neither same-q rearming nor general deflation

In (1), δ>0 occurs only at u=0. Then q is odd. At u≥1, δ=0 and q is even.
Thus the entire coarsened carrier satisfies the sharp separation

    δ>0 ⇒ q odd,        δ=0 ⇒ q even.

A same-q rearm from δ0 toδ1 after L, needed by our conditional binary-digit
lemma, is therefore impossible **inside this carrier**. R11 preserves its
positive boundary defect while creating an R-port; it does not regenerate
that defect after it has been erased by L. Coarsening the twin block changes
the presentation, not this parity separation or the actual response q.

There is useful but limited internal halving. If u≥2, reversing L produces
q/2 and remains at δ0. At u=1, reversing L instead produces

    (q/2−e, δ=e),

not q/2. For instance the last inverse step is60→29, not60→30, and
100→47, not100→50. At g=1 the corresponding steps are156→77 and260→127.
This locates the boundary correction exactly. It explains why arbitrary
internal powers of two do not imply an odd-head deflation theorem.
These are carrier restrictions, not all-tournament exclusions of30,50,78,
or130.

The odd slopes of the enlarged coarsened carrier form exactly the following
eight sequences, indexed by v≥0:

    15·3^v, 39·3^v, 25·3^v, 65·3^v,
    30·3^v−1, 78·3^v−1, 50·3^v−3, 130·3^v−3.

This absorbs the redundancies in our former twelve-sequence presentation
and adds its four missing boundary starts. It remains O(log X) slopes≤X.
For the four odd boundary sequences the whole s≥1 ray is already supplied
by the tail identity. That local head completion does not extend to all
odd slopes or to arbitrary δ0 states.

## 4. The third residual grid is a direct head replacement

Incoming R10 proves a new fixed-port seed family

    U(a)=35·2^a+32,

and, after a=r+5, covers H=1+(35·2^r+1)2^s for s≥6. Its separate
Hamilton kernel covers s=1 for every r≥1. The uniform levels s=2,3,4,5
remain open in that packet. This is correctly distinct from the already
complete d=3 and d=21 grids.

Neither operation is our missing same-q rearm. The seed-port update is
U→2U−32, a different physical interface, and the s=1 kernel constructs a
new tournament with the desired H. It does not assert a unit seed with
multiplier35·2^r+1 or deflate the multiplier32(35·2^r+1).

It nevertheless suggests a weaker finishing requirement than universal
seed deflation. Suppose the finite-bank/rearm criterion of coverage/PROOF.md
has produced every even multiplier q≥2B. Then every count1+4p with p≥B
already exists. To finish all sufficiently large odd counts, it is enough
to give **direct Hamilton constructors** for1+2p for every odd p≥B;
these need not arise by deflating a unit seed. Together the two statements
cover every odd N≥4B+1, leaving a finite head check. R10 exemplifies this
second, replacement-style strategy on the single family p=35·2^r+1.
It does not establish the universal odd-p hypothesis.

## 5. Preserved conclusion

The incoming work strengthens the physical interpretation: selected births
can activate fresh stable ports at genuine deletion roots, and separate
Hamilton kernels can complete a head without any seed-halving operation.
After merging the marked twin presentations, the exact rearming obstruction
is even clearer. No same-q rearm, universal odd deflation, finite complete
binary seed bank, or Chan–Pak proof is supplied by these incoming results.
