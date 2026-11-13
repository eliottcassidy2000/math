# Uniform-module root transport: exact target, strict convexity, and arithmetic limits

**Status: PROVED here, with FINITE-EXACT independent controls.** This is a theorem
about an actual marked tournament and a specified family of neutral arc flips.
It is not a classification of the full Hamiltonian-count fibre or the native
Hamilton-path complement metagraph. Chan–Pak remains open.

The inherited operation is the marked transitive-module preparation in
[Session 298](../../session298_exact_target_geometry_20261113/README.md), extending
the H17 construction. The new statements below apply to **every fixed uniform
exterior**, classify the entire root-position response by positive run kernels,
and prove a metric law, nonnegative forward differences of every order, and
strict convexity. The two-width formulas are special cases.
The kernels were first derived in this session's reusable-repairs lane; their
general-context proof and consequences are recorded here.

## 1. Actual objects and the preserved target

Let Q be any finite tournament with a distinguished vertex v. Replace v by
a transitive module C=TT_m: every outside vertex sees every vertex of C in
the same direction, as prescribed by Q. Write E=Q-v. Fix a named root o in C.
Keep the relative order of the other m-1 named vertices fixed and insert o
after l of them, where r=m-1-l remain after it. Call this marked tournament T_l.

Use the exact state

```
h = H(T_l),             s_l = Hamilton paths of T_l starting at o,
q = H(T_l-o),           F_l = ordered spanning pairs (P,Q'):
                              P nonempty, ending at o; Q' possibly empty.
a_l = h-s_l,            d_l = F_l-q,
D_l = d_l-a_l = F_l+s_l-h-q.
```

H(empty)=1. The letter Q' denotes a path, not the quotient tournament Q.
**D is the controller imbalance; it is not the ready-state defect q-s.**

All T_l have the same unrooted isomorphism type. Deleting the named root gives
the same labelled tournament at every l. Consequently **the full target h
and the full deletion count q are exactly constant**. Moving l to l+1 changes
only the arc between o and its next neighbour in the total order. These two
vertices have identical incidences to every other vertex at that step.
Thus each preparation is a literal one-arc operation with no uncounted paths.

If Q is strong and has at least three vertices, all these substitutions are
strong. Outside vertices communicate by paths lifted from Q; paths from one
module vertex to another can leave through a directed cycle of Q through v
and return. Strongness is not needed for the response theorem itself.

## 2. Positive run-kernel representation

Define, for j>=1 and l,r>=0,

```
S_j(l,r) = sum_{i=0}^{j-1} (-1)^(j-1-i) binom(j-1,i) i^l (i+1)^r,
E_j(l,r) = S_j(r,l),
```

using 0^0=1. These are nonnegative integers. S_j counts covers of TT_(l+r+1)
by j ordered, nonempty increasing paths with the distinguished first path
starting at o. Place o in bin 0. Every predecessor of o must be assigned to
one of the other j-1 bins; every successor can use any bin. Inclusion–exclusion
requires each other bin to be nonempty and yields S_j. Reversing the total
order gives E_j, for the distinguished path ending at o.

There are nonnegative integers C_j and B_j, depending only on the exterior
Q and its marked module position, such that for **every** m,l,r,

```
s_l = sum_j C_j S_j(l,r),          F_l = sum_j B_j E_j(l,r).        (1)
```

Only j<=|E|+1 is needed in the first sum and j<=|E|+2 in the second.
To prove (1), take each actual Hamilton path or two-path cover and contract
each maximal consecutive C-run to a module symbol. Record the order of all
outside vertices, the path divider for F, and which module run contains o.
For s that run is first; for F it is last in the first path. Number that run
0 and number the other runs in their skeleton order. The outside incidences
are independent of the choice of vertices assigned to the runs precisely
because C is a module. Each valid skeleton contributes exactly S_j or E_j.
The contraction and expansion are inverse, so no realization is omitted or
double-counted. Empty global second paths are retained. Adjacent module
symbols across the path divider are allowed; within a single path they would
form one maximal run and are excluded.

Two coefficients have uniform positive lower bounds:

```
B_1 >= H(E) >= 1,                  B_2 >= H(Q) >= 1.             (2)
```

For B_1 put the whole module in P and a Hamilton path of E in Q'. For B_2
put the distinguished module run alone in P, and put a second module run
into each Hamilton path of the quotient Q to form Q'. The elementary
insertion argument gives a Hamilton path in every finite tournament:
insert a new vertex before the first path vertex it beats, or append it if
there is none. The required arc from its predecessor follows from minimality.
The empty exterior has one empty Hamilton path by convention.

## 3. All forward differences are nonnegative; exact contraction and equality

The adjacent-twin bijections were already proved in
[O4, return6325 observer](../../return6325_observer_20261103/PROOFS.md): for
adjacent homogeneous vertices u->v and their marked contraction C,

```
F_v(T)-F_u(T) = F_root(C),       s_u(T)-s_v(T) = s_root(C).
```

This is inherited, not a new twin identity. Swapping u and v injects the
covers ending at u into those ending at v; the omitted covers use u->v and
contract to covers ending at the contracted root. The start-path identity
is the dual injection. The new all-order consequences below combine this
physical contraction with the fixed-size H-neutral root-phase orbit.

Write T(m,l) for the marked graph at phase l in a module of size m, and
Delta for the fixed-size step l->l+1. Iterating the inherited identities gives

```
Delta^t F(m,l) = F(m-t,l),
Delta^t s(m,l) = (-1)^t s(m-t,l),
Delta^t D(m,l) = F(m-t,l)+(-1)^t s(m-t,l),                 (3)
```

for t>=1 and l+t<=m-1. Contraction removes a module successor; exterior
incidences and the root's number l of predecessors are unchanged. The
constants h,q disappear after the first difference. Equivalently, directly
from the positive kernels in Section2,

```
Delta E_j(l,r) = E_j(l,r-1),       Delta S_j(l,r) = -S_j(l,r-1),
Delta^t D(l,r) = sum_j B_j E_j(l,r-t)
                 + (-1)^t sum_j C_j S_j(l,r-t).
```

For every nonempty rooted tournament U,

```
F_root(U) >= H(U-root) >= s_root(U).
```

The first inequality uses the covers whose first path is the singleton root;
the second deletes the first vertex of a path starting there. If the root
has an in-neighbour u, the additional covers P=(u,root), Q'=any path of the
remaining tournament prove F>H(U-root). If the root is a source, every cover
ending there has singleton first path and every deletion path can be prefixed
by the root, so F=H(U-root)=s. Consequently:

- **Every forward difference of D of positive order is nonnegative.**
- Every even-order difference is strictly positive, since F+s>=1.
- An odd-order difference vanishes exactly when the root of T(m-t,l) is a
  source: equivalently l=0 and the original module has no outside in-neighbour.
- In particular, if the quotient Q is strong with at least three vertices,
  **all** existing positive-order differences of D are strictly positive.

Thus D is nondecreasing and strictly discretely convex when m>=3. If an
outside in-neighbour exists, D is strictly increasing. If none exists, its
first two entries are equal (when m>=2), and it strictly increases thereafter.
Its minimum is at the first position, its maximum at the last; its sublevel
sets are initial intervals. A value D* is attained at at most two positions,
with a double occurrence only in the first two positions in the source case.
In a strong context a value is attained at most once. Binary search finds
all positions with D=D* using O(log m) state evaluations; this is an oracle
count, not a claim that each arbitrary-context state evaluation is cheap.

The boundary m=1 has no differences; m=2 has no second difference. At the
contracted singleton F=s=q=1, so odd differences vanish and even differences
are 2. No strong-context assertion includes that source case. Uniform exterior
incidences remain essential.

### Exact reconstruction from first-position states

Newton's finite-difference identity yields

```
F(m,l) = sum_{t=0}^l binom(l,t) F(m-t,0),
s(m,l) = sum_{t=0}^l (-1)^t binom(l,t) s(m-t,0).           (3a)
```

For a fixed exterior and module sizes 1..M, all M(M+1)/2 marked root-position
states can therefore be reconstructed from the M first-position states,
using (3a) and the known constants h_m and q_m. These constants are already
contained in the first-position records. This reduces expensive full graph
counts to M evaluations, followed by exact integer arithmetic. The contracted
states are a calculation table; their smaller H-values are **not** asserted
to be target-preserving physical repair steps. The preparation itself still
uses the fixed-size adjacent flips of Section1.

## 4. The neutral path is isometric

Give marked tournaments the arc-edit distance, allowing endpoint relabellings
only when they fix the named root. The root's indegree in T_l is a fixed
exterior contribution plus l. One arc reversal changes root indegree by at
most one, so any route from T_l to T_k needs at least |l-k| reversals, even
if intermediate graphs leave this family. The adjacent moves above attain
that bound while preserving h at every step. Thus the certified root orbit
contains an **isometric integer interval** in the exact-H arc-flip graph.

This does not claim that every pair in an H-fibre can be joined neutrally.
The native metagraph has a different edge operation; no native metagraph
connectivity conclusion follows without compiling each step into its edges.

### Canonical coordinate integration, with explicit prior credit

[THM-4463, homogeneous tournament contacts](../../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md),
Section2, already proves that the graph of homogeneous pairs is a disjoint
union of consistently directed paths, exactly the maximal transitive modules.
[Weighted pair-free cores](../../weighted_core_transport_20261020/core_count.md)
already proves the unique decomposition T=Q[w] into a pair-free tournament Q
and positive chain lengths w, modulo Aut(Q), including
Aut(T)=Stab_(Aut Q)(w). These structural and enumeration facts are inherited,
not new results of this session. A pair-free core can have nontransitive
modules; it need not be indecomposable.

The new response law can therefore be applied canonically, with no guessed
block boundary. For any fixed unrooted T=Q[w], its root positions are the pairs
(i,l), 0<=l<w_i, modulo the weight-preserving automorphisms of Q. Such an
automorphism preserves l because a finite transitive chain has only the
identity order automorphism. On each fixed i, the full response satisfies
(1)–(3) in its actual exterior, and its named-root transport is the isometric
path just proved. This classifies these coordinates for **all root choices**,
while not claiming a neutral path between different blocks. A singleton block
has no internal preparation. The quotient arcs and the other block lengths
are required sidecars for its numerical response.

Even for four vertices, homogeneous-pair flips are a strict subset of neutral
arc flips. In out-neighbour bitmasks on vertices0..3,

```
[8,1,3,6] -- reverse arc(0,2) --> [12,1,2,6]
```

has H=5 on both sides, although vertices0 and2 have different incidences to
vertex1. Thus the canonical block atlas does not exhaust H-neutral moves or
the changing-core repair routes.

## 5. Arithmetic selection defeats scalar extremization

The inherited physical return R^k N R has

```
t = 2^k d+q-D,                   b = 3t+2D.                     (4)
```

Its balanced output is primitive exactly when gcd(t,2D)=1. Native entry also
requires t divisible by3 and b not divisible by3, in addition to positivity.
Convexity optimizes D; it does not optimize those residue and gcd conditions.

A whole hostile family is already implicit in Session 298. Take
Q=C3 and m odd, m>=3, so h=2^m+1. At l=0, D=0 and b=3t, excluding primitive
entry. At l=m-1,

```
D = 2*3^(m-1)-2^(m-1)-2 == 0 (mod 3),
```

so b is divisible by3 for every k and native entry is impossible. The interior
l=1 has D=2+2^(m-2)==1 (mod3), and the cited Session 298 argument provides
infinitely many primitive native entries with even k. Thus an interior phase
works while both extreme positions fail. At the smallest case m=3:

| l | (h,s,q,F) | D |
|---|---|---:|
| 0 | (9,4,5,10) | 0 |
| 1 | (9,2,5,16) | 4 |
| 2 | (9,1,5,25) | 12 |

For l=1,k=2, (b,t)=(143,45); this passes the primitive native entry guard.
This example does not identify the ready-state defect with D and does not
assert that the final target equals the initial h. The **preparation** fixes h;
the subsequent return (4) has the stated new target b.

## 6. Hostile controls and exact reproduction

The smallest failure without module homogeneity has three vertices:
0->1->2->0, with candidate module {0,1}. Reversing 0->1 produces the
transitive tournament 1->2->0, with 1->0. The exterior vertex2 distinguishes
the two module vertices. The Hamiltonian count changes from3 to1 under that
single flip. A selective port requires its own enlarged skeleton/sidecar;
one cannot reuse this uniform-module theorem unchanged.

Run from this directory:

```
python -B root_phase.py
python -O -B root_phase.py
```

The explicit universe is all labelled contexts with 0..3 exterior vertices,
all uniform module incidences, module sizes 1..6, and all root positions; plus
16 deterministic contexts with 4 exterior vertices and module sizes 2..5.
There are 1,799 tournaments, 547 independent literal-permutation state checks,
1,285 single-arc steps, 1,087 strict even-difference checks, and 91 kernel
assignment checks. A further 2,945 exact contraction-difference identities
include 114 odd-order source equalities; all 1,799 marked states pass the
binomial reconstruction. Both four-vertex non-twin neutral states have separate
literal-permutation controls. These finite controls support but do not replace the
all-context combinatorial proof. No assert statement is used, so optimized
Python retains every check. Exact reports are checked in alongside the script.

## 7. Actual contribution to the open construction

The theorem supplies a reusable exact-target preparation **whenever a suitable
transitive module is already present**. It reduces the whole position scan to
positive kernels, supplies an exact neutral path, and separates convex scalar
optimization from arithmetic entry. It neither manufactures such a module at
every desired target nor proves that a good arithmetic phase always exists.
Changing the module size changes h in general. The remaining source-supply
and selective-port questions are genuine obligations, not consequences of
the convexity theorem.
