# Complete positive joint C2 preparation locus on all tournaments

**PROVED — [independent audit](JOINT_AUDIT.md) passed.** For every finite
tournament, exact positive preparation at a retained root forces
the auxiliary to have that root as its sole out-neighbour. The desired C2
digits2 and4 are therefore characterized by deletion-ending counts1 and2.
This classifies the complete joint input guard; it does not construct a return
of actual outputs to that class.

## 1. Inheritance and statement

The exact operator and all endpoint conventions are inherited from
[construction/RESULT.md §§1–3](../../multiport_defect_grammar_20261103/construction/RESULT.md).
The prior scalar-ready classification is
[READY.md §§1–3](../../two_entry_funnels_20261103/READY.md), read in converse.
That scalar classification includes exceptional two-exit auxiliaries with
sigma=2. The present joint statement removes those exceptional cases from
the `tau=0` interface; it does not retract their scalar readiness.

The closest mechanism is the cut-before-auxiliary injection underlying the
nonnegativity of tau. The canonical hostile is the two-exit scalar-ready
vertex in `(8,1,3,6)`, which has sigma2 but positive tau at every retained
root. The corrected near miss is outdegree one with the wrong exit: it can
have sigma=tau=0, so **positive** preparation is essential. The least-used
sidecar is the failed joining endpoint of the first free path. The live
concepts are actual arcs, paired endpoint slots, shortest directed access,
deletion ending counts, zero versus positive defect, and physical return.

Let T be a nonempty finite tournament, a its auxiliary and w its
retained root; a=w is initially allowed. In the square-zero vertex ring let
P include the empty path; A_a count nonempty paths starting a; W_w count
nonempty paths ending w; Q_aw count paths starting a and ending w. Put

    H=[P], F=[W_w P], S=[P A_a],
    G=[W_w P A_a]+[P²Q_aw], sigma=S−H, tau=G−F.      (1)

The proved C2 operator sends `(H,F)` to `(3H+2sigma,3F+2tau)`.
Every product below uses **ordered** factor slots. Write B=N^+(a), d=|B|,
and E_B=sum_(b in B) E_b, where E_b counts paths ending b.

**Theorem.**

1. On every finite tournament, with arbitrary coincident or distinct marks,

       tau=[E_B W_w A_a]+[E_B P Q_aw].                (2)

2. Whenever a reaches w,

       tau>=d−1 if a!=w, and tau>=d if a=w.           (3)

   If a does not reach w, tau>=sigma.

3. On **all finite tournaments**, the complete positive joint preparation locus is

       sigma>0 and tau=0
       iff a!=w, N^+(a)={w}, and e_w(T−a)>0.         (4)

   In that locus, sigma=e_w(T−a). In particular C2 has exact retained-root
   response `(3H+2,3F)` or `(3H+4,3F)` iff the auxiliary has sole exit w and
   the actual ending count in T−a is respectively1 or2.

No induced-extension, order, labelled-class or core-bank cutoff appears in
the theorem. The singleton obeys the formulas and has sigma=tau=0, so it
does not enter the positive locus. The later SCC table explicitly specializes
to nontrivial strong T (order at least3).

## 2. Exact decomposition by failed joins

Start with an F-cover `(W,V)`, W nonempty ending at w and V possibly empty.
If a lies in V, cut V immediately before a. This gives the ordered L-cover
`(W,P,A_a)`. If a lies in W, cut W immediately before a; its suffix is Q_aw,
and the ordered K-cover is `(P_before,V,Q_aw)`. The two destinations are
disjointly typed; restoring the cut recovers the input.

An L-cover lies in this image exactly when its P slot is empty or the last
vertex of P beats a. Otherwise its last vertex belongs to B. Thus the exact
complement of the image in L is `[W_w E_B A_a]`.
Likewise a K-cover lies in the image exactly when its first P slot is empty
or its last vertex beats a. Its complement is `[E_B P Q_aw]`. This proves
(2), including every empty-slot and coincident-mark case. For a=w, Q_aa is
the singleton path a, while the L term vanishes because W_a and A_a would
both use a.

This is a nonnegative combinatorial decomposition, not cancellation of
uncontrolled marginal counts. Its source is the full actual path-cover set;
its target records the failed joining endpoint and the paired a-to-w path.
Forgetting that endpoint loses the strictness mechanism below.

The same cut-before-a argument for a Hamilton path gives the exact scalar
failed-join identity

    sigma=[E_B A_a].                                  (2a)

In particular d=0 forces sigma=0. Suppose a cannot reach w. In every cover
of (2a), w must lie in the E_B path, since it cannot lie on a path starting
at a. It is not that path's final vertex b, because a→b would give access
to w=b. Cut the E_B path immediately after w. Its two nonempty pieces end
at w and b respectively; together with the original A_a path they give an
L-surplus cover in (2). Joining those pieces recovers the original cover,
so this is an injection. Consequently tau>=sigma in the nonreachability
case. This is the extra step needed for the all-tournament scope; it does
not assume an endpoint is reachable merely because the graph is named.

## 3. A shortest path pays a strict lower bound

For a!=w that a can reach, choose a shortest directed path Q from a to w.
It contains at most one vertex of B: if a later vertex also lies
in B, the direct arc from a to it shortens the path. Hence at least d−1
members b of B lie outside Q.

For each such b, place the singleton b in E_B, put Q in Q_aw, and put any
Hamilton path of the remaining induced tournament in P. The empty remainder
uses P's unit. Every tournament has a Hamilton path by the elementary
vertex-insertion argument. These are distinct K-surplus covers because their
first singleton endpoint is different. Thus tau>=d−1.

For a=w use Q=(a). All d exits lie outside it, and the same construction
gives tau>=d. Under positive sigma, (2a) gives d>=1, so a=w is impossible
when tau=0. If a!=w and sigma>0,tau=0, the nonreachability injection first
forces a to reach w; the bound then forces d=1.

The witness and the nonreachability injection together pay every tournament,
not only strongly connected inputs.

## 4. The remaining unique-exit case

Let N^+(a)={b}, a!=w, and U=T−a. On U write E_b for paths ending b and A_b
for paths starting b. Directly partitioning paths at a gives

    P_T=P+z_a(P−E_b+P A_b),
    W_(w,T)=W_w+z_a P Q_bw,
    A_(a,T)=z_a(1+A_b), Q_(aw,T)=z_a Q_bw.

Extracting (1) yields

    sigma=e_b(U), tau=[W_w E_b]_U.                  (5)

This rederives the needed identities here; no unpromoted companion theorem
is used as a dependency. If b=w, the product in (5) repeats w and vanishes.
If b!=w, cutting each HP of U ending b immediately after w gives an ordered
two-path cover ending w and b. Both slots are nonempty, and restoring the
join recovers the HP. Therefore tau>=e_b(U)=sigma. A positive sigma with
tau=0 forces b=w. Conversely b=w gives tau=0 and the stated sigma exactly.
Together with §§2–3 this proves (4) in its all-tournament scope.

Thus the exceptional two-exit scalar-ready language of READY can never
provide the full retained-root guard, even if w is freely chosen among
every actual vertex. This is an all-order joint consequence, not a search
failure for the first few exceptional words.

## 5. Equality boundary and physical controls

Both desired guards occur. Literal examples are

|Rows|Auxiliary a|Root w|H|F|sigma|tau|
|---|---:|---:|---:|---:|---:|---:|
|(8,1,3,6)|1|0|5|7|1|0|
|(20,1,2,7,14)|2|1|9|21|2|0|

In the first graph, auxiliary3 has two exits and sigma2, but its tau values
at roots0,1,2,3 are respectively7,1,1,3. Thus scalar readiness does not pay
the joint guard, and the d−1 bound can be sharp.

The positivity boundary is actual. Perform C2 at old auxiliary0 of the
triangle `(2,4,1)`, retaining old root1. Its output is

    (18,28,25,1,8).

The transported auxiliary3 has sole exit old0, not root1, yet sigma=tau=0.
The theorem therefore cannot replace “positive preparation” by “zero joint
defect.” This is precisely the inherited zero-channel behaviour, rather
than an abstract unattained hostile.

No statement here makes a preparation H-neutral, preserves F under a preceding
surgery, or says that either child lies in (4) again. A genuine reusable
two-branch grammar still needs the actual return operation, its changing
observer, and a closure certificate. The theorem identifies exactly which
positive C2 input interface such a return must create.

### Corollary: a fixed retained root has at most one positive auxiliary

Suppose a and a' were distinct positive joint-prepared auxiliaries for the
same actual retained root w. By (4), both would have w as their sole exit.
The arc between a and a' gives one of them a second exit, a contradiction.
Thus at most one auxiliary of an unchanged marked source can provide any
positive C2 preparation. In particular the same `(T,w)` cannot supply both
digits2 and4 by auxiliary choice alone. This leaves preceding source surgery
and a genuinely different marked source available. Strongness is not needed
for this uniqueness consequence.

### Corollary: the entire joint-zero locus has an exact SCC description

First, for **any** tournament U with distinct b,w,

    [W_w E_b]_U>0 iff {b,w} meets the terminal SCC of U.    (6)

Necessity follows because a directed path containing a vertex of the terminal
SCC cannot leave that component; one of the two path endpoints must be there.
For sufficiency, a Hamilton path can end at any specified vertex of the
terminal SCC: take Hamilton paths of the earlier SCCs in their unique order,
followed by a path of the terminal SCC ending at that vertex. Such an ending
path exists in a strong tournament by the endpoint-preserving insertion
argument in construction/RESULT §4. Cut that full Hamilton path immediately
after the other named endpoint. This produces the two required nonempty
paths and proves (6).

Return to a nontrivial strong T. By (3), tau=0 requires a!=w and a unique exit
b. In U=T−a, the vertex b belongs to the **initial** SCC. Indeed a simple path
from a to any vertex of U must begin a→b and cannot revisit a, so b reaches
every vertex within U. That is precisely membership in the initial SCC of
a tournament's totally ordered SCC decomposition.

If U is strong, e_b(U)>0; (5) and (6) give tau=0 exactly when w=b. If U is
reducible, its initial and terminal SCCs differ, so b is outside the terminal
one and e_b(U)=sigma=0. For w!=b, (5)–(6) give tau=0 exactly when w is also
outside the terminal SCC; for w=b the same conclusion holds since the product
repeats b and vanishes. The complete table is therefore:

|Actual case in strong T|Condition for tau=0|sigma|
|---|---|---|
|a=w|never|not needed|
|a!=w, outdegree(a)>=2|never|not needed|
|sole exit b, U strong|w=b|e_b(U)>0|
|sole exit b, U reducible|w outside the terminal SCC of U|0|

There is no omitted outdegree-zero case on a nontrivial strong T. This table
also proves that positive joint preparation has a strong deletion parent.
The inherited C2 zero-channel is visible structurally: deleting its newly
transported auxiliary x leaves `old T→y`. The old exit and retained old root
both lie outside the terminal singleton y, so the last row applies. No
Hamilton count is needed to recognize this zero-joint sidecar.

## 6. Standalone exact controls

The source `04-computation/nonzero_preparation_20261103/algebra_joint_guard.py`
imports no prior counter or receipt. It enumerates every labelled tournament
through order5 and every ordered `(a,w)`, including coincident marks: 26,705
marked cases in total. Directed simple paths are enumerated by DFS with both
endpoints, and full square-zero convolution compares the two sides of (2).
For every reachable shortest-path witness, the code checks an omitted exit,
the actual a-to-w path and a Hamilton path of its remaining vertices.

The all-strong subuniverse contains 2,24,544 graphs of orders3,4,5. It checks
(3) and both directions of (4) for every marked pair. The separate actual
transported-zero control and both positive digits are retained. These finite
checks test the proof; they do not replace its shortest-path argument.

```text
python3 -B 04-computation/nonzero_preparation_20261103/algebra_joint_guard.py --output /tmp/algebra_joint_guard.normal.json
python3 -B -O 04-computation/nonzero_preparation_20261103/algebra_joint_guard.py --output /tmp/algebra_joint_guard.optimized.json
```

Both modes pass **110,927 always-active gates** and produce identical
1,855-byte receipts. Frozen source SHA256:
`219e9da4aeea70cd95eda6fee4ad9a5b2704602e1dd9e75aa02a78cbfb2c7233`.
Frozen receipt SHA256:
`148a19cf1fd9bc54a75576782fc700030f3d53b27dcf3febb8c66e289224396a`.
The independent analytic/source/replay audit has passed.

The separate source `algebra_joint_scc.py` independently enumerates paths by
permutations, rather than importing the core source's DFS counter. Its finite
universe is every labelled tournament of orders2–5 and all 10,650 unordered
distinct ending pairs, plus all 14,002 ordered marked pairs on the strong
subuniverse. It checks (6), the complete table, the initial-SCC claim for each
sole exit, and uniqueness at every fixed root. Both normal and optimized runs
pass **26,851 gates** and produce identical receipts:

```text
python3 -B 04-computation/nonzero_preparation_20261103/algebra_joint_scc.py --output /tmp/algebra_joint_scc.normal.json
python3 -B -O 04-computation/nonzero_preparation_20261103/algebra_joint_scc.py --output /tmp/algebra_joint_scc.optimized.json
```

The sidecar source is 5,039 bytes, SHA256
`2fa1b5ee13f032f00e168ab41dfbe1f6785d441553b1ed9a27c727da7d096f11`;
its receipt is 581 bytes, SHA256
`d136b727d13ef78b877e613a9b20a8d02f778cc3cf60dfc6bf6ce44870e40052`.
The core source and its 110,927-gate receipt remain unchanged.

The all-tournament upgrade has its own independent standalone
`algebra_joint_weak.py`. It enumerates path objects by permutations and checks
the scalar failed-join identity, both directions of (4) on all 26,705 marked
cases, and the actual nonreachability injection by cutting/rejoining every
eligible path object. Thus this added scope is not inferred from the earlier
strong-only guard assertions. Reproduce by

```text
python3 -B 04-computation/nonzero_preparation_20261103/algebra_joint_weak.py --output /tmp/algebra_joint_weak.normal.json
python3 -B -O 04-computation/nonzero_preparation_20261103/algebra_joint_weak.py --output /tmp/algebra_joint_weak.optimized.json
```

Both modes pass **47,221 gates** and produce identical 555-byte receipts. They
check 3,212 unreachable marked pairs and 3,966 literal injection images, and
retain 312 positive weak guards. Frozen weak source is 5,106 bytes, SHA256
`fb56ac76151a3411c22728a9018f7d1aaa9a355140c53c3cc3a6d5d0e5c4d3cd`;
receipt SHA256 is
`3103157b1b6214be3326a83309af6d8d1ac9d1d8aa7c8ff24c60ccf80da19a7f`.
The core and SCC sidecar remain unchanged.
