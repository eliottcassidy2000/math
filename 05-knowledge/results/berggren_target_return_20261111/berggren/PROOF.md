# Berggren source charts and an exact neutral three-arc return

**PROVED in this note:** a completely classified infinite family of physical,
target-preserving three-arc rewrites after a `3H+2` exit; its full source and
return interfaces; its failure to restore small readiness. **INHERITED and
independently rechecked:** the one-apex classification, enumeration and
four-dimensional linear count carrier. **OPEN:** a repeatable target-preserving
ready return, coverage of the Chan--Pak spectrum, and the Chan--Pak conjecture.
No literature-priority claim is made.

The closest mechanism is Session256's [Berggren substitution lift](../../session256_chan_pak_triples_20261111/spinor/PROOF.md).
The hostile example is its three decorated H9 sources. The corrected near miss
is forgetting the apex before applying a branch. The underused sidecar is the
October [complete apex classification](../../session26_random_bridges_20261022/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md),
reinforced by the [fixed-H pointed/unpointed correction](../../session133_fractran_spectrum_refinement_20261101/path_counter/THEORY.md).
These are earlier results, not new classification claims from this session.

## 1. The actual objects behind 107, 73 and 17

For a binary word `w`, T(w) has a transitive base in word order and an apex
`a`, with `a->v` at O positions and `v->a` at I positions. The inherited
Berggren letters are

```
A: O -> OI,  I -> I
B: O -> OI,  I -> OIO
C: O -> O,   I -> OOI.
```

They act on decorated graphs, preserving the specified base order and apex.
The starting word is OI. The exact collision data are:

| Address | Word | (H,P,last bit,1) | Strong? | A-child H |
|---|---|---|---|---|
| B | OIOIO | (9,21,1,1) | no: the last O is a sink | 107 |
| C | OOOI | (9,14,0,1) | yes | 73 |
| AA | OIII | (9,8,0,1) | yes | 17 |

Here P reads O as 1 and I as 0. The last two rows are **the same ordinary
tournament in different apex charts**. Explicitly, with apex last, map the
vertices `(0,1,2,3,4)` of T(OOOI) to `(1,2,3,4,0)` of T(OIII). Every arc is
preserved. Both are the cyclic blowup of a transitive triple and two
singletons. The first row has strong core T(OIOI), which is a different H9
ordinary class. Thus there are three decorated sources here, but not three
different strong ordinary sources.

This is also the explanation of the infinite collision
`H(C^k)=H(A^(2k))`: the words `O^(2k+1) I` and `O I^(2k+1)` are the two
apex charts of the same ordinary cyclic blowup. The tree remains an honest
decorated construction; neither ordinary isomorphism nor H alone respects
its branch operation.

## 2. Four coordinates: exact meaning and previous work

The inherited append updates, in homogeneous coordinates `(H,P,e,1)`, are

```
O: (H,P,e,1) -> (H,2P+1,1,1)
I: (H,P,e,1) -> (H+P+e,2P,0,1).
```

The four-dimensional Hankel rank is already recorded in
[the G(1,0) dictionary, Section 3.2](../../spine_grammars_cutwidth_apexes_hankel_seams_20261104.md).
The growing bit-memory distinction and exact two-letter probes are older
[Session176 results](../../session176_tournament_information_20261104/gates/NOTE.md).
For a short self-contained rank certificate in the present convention, take
prefixes `empty,O,OI,OO` and suffixes `empty,I,II,OI`. The matrix of actual
Hamilton counts of their concatenations is

```
1 1  1 3
1 3  5 5
3 5  9 9
1 5 11 9
```

Its determinant is 16. Every homogeneous linear word realization of dimension
d factors this matrix through a d-dimensional space, so d is at least four.
The displayed carrier attains four. This is minimality for homogeneous
linear realizations of the **append-word series**, not a lower bound on the
number of nonlinear registers or the matrix-pair substitution compiler.
Nonlinearly, e=P mod2 and the constant is known, so `(H,P)` suffices for
append evaluation. For arbitrary words P alone also loses leading I's,
which matter to the B/C substitutions. On the Berggren image the initial O
removes that leading-zero ambiguity.

Use Q(w) for the reverse-complement binary value:

```
Q(w) = sum_(i=0)^(|w|-1) [w_i=I] 2^i.
```

The inherited concatenation identity, in these coordinates, is

```
H(pv)=H(p)+H(v)-1+P(p) Q(v)+e(p)[v begins I].       (1)
```

The empty suffix has H=1,Q=0 and beginning-I indicator zero. One proof is
to propagate the difference between the initial states for v: the coefficient
of incoming P is Q(v), while incoming e is used only when the first letter
is I. This is a coordinate version of the older one-apex concatenation
calculus, not a new general compiler theorem.

For nonempty w, the full rooted profile at its apex is

```
(H,E,S,U,V,d,K)
  = (H(w), [last I], [first O], P(w)+1, Q(w)+1, 1, 2^|w|).     (2)
```

Here E and S count paths ending and starting at the root; U and V count
ordered two-bin covers whose root bin respectively starts or ends there;
the other bin can be empty; d=H(T-root); K is the full two-bin cover count
of the deletion graph. For example, a root-starting path in a selected base
subset must begin at that subset's first O. Summing over subsets gives
P+1; reversal gives Q+1. This fixes the orientation convention explicitly.

## 3. The recovered complete class family, rather than a new census claim

The cited October classification proves that a strong word starts O and ends
I. At length at least three its admissible deletion apices are:

* one, when both letters occur at least twice;
* two, for `O^r I` or `O I^r`, r>=2; these two words represent one ordinary class.

C3 is the separate length-two case with three admissible apices but only
one pointed isomorphism class. Every strong one-apex graph of order at least
four is rigid. Indeed, an automorphism fixing an apex fixes its unique
transitive deletion order; in the two-apex case their outdegrees r and 1
are different. Every nonstrong graph in this family is a unique ordinal
sum `TT_a -> strong core -> TT_b`, or is transitive.

Consequently, for ordinary classes, the inherited formulas are

```
c_3=1;  c_n=2^(n-3)-1 for n>=4,
a_n=2^(n-1)-1-binomial(n-1,2) for n>=2;  a_1=1.    (3)
```

Only the n-2 classes with a C3 core have automorphism group of order three;
all other classes are rigid. Hence the labelled count for n>=3 is
`n! (a_n - 2(n-2)/3)`. The independent exhaustive controls here find 13,200
of the 32,768 labelled order-six tournaments in the one-apex class, forming
21 ordinary classes, seven strong. This organizes a completely understood
infinite subfamily of A000568. It does not give new terms of A000568 or a
classification of all tournaments.

## 4. Every strong C-image has a genuine ternary exit

Let w=O u, where u itself starts O and ends I. In T(w), the first base vertex
r=0 has only the apex as an in-neighbour. Deleting r gives T(u). Write
`q=H(u)` and `Q=Q(u)`. Equation (1) gives H(w)=q+Q.
The number of paths of T(u) starting at its apex is one. Deleting r from a
path starting at r gives all paths of T(u) except that one. Thus

```
indegree(r)=1,  H(T-r)=q,  s_T(r)=q-1,  F_T(r)-H(T)=1.        (4)
```

The last equality is the inherited indegree-one defect identity from
[the complete ready-root theorem](../../two_entry_funnels_20261103/READY.md).
Add x beaten only by r, then y beaten only by x. The inherited two-birth
law therefore gives H(G)=3H(T)+2. Applying the C substitution to any strong
apex word yields another strong word beginning OO, so it always supplies
such a literal exit.

Deleting the **old apex** from G gives

```
C3(r,x,y) -> TT_(|u|),    r->x->y->r.              (5)
```

This is an actual H3 deletion-core chart, with the old apex beating only r
inside the triangle, and seeing the tail through u. It connects the one-apex
carrier to the H3 member of the five-core program. It is not an identification
of the two interfaces.

There are two vertex-disjoint triangles in G: `(a,1,last base vertex)` and
`(r,x,y)`. Deleting `{a,r}` leaves a transitive graph. Thus G has feedback
vertex number exactly two, at every size in this family. Re-marking its
vertices cannot turn this same graph into a one-apex tournament.

The full profile at the retained old apex a is

```
(3H+2, 3, 1, 3(P(w)+1), 6Q(w)+9, 3, 6*2^|w|).   (6)
```

For clarity, the two-cover factors come directly from (5). For root-start
covers, assigning no core vertices to the root bin contributes
`3(P(u)+1)`; the three allowable nonempty core choices each contribute
`2^|u|`. This is `3(P(w)+1)`. For root-end covers, a nonempty selected tail
ending at an I contributes `12Q(u)`: 12 is the two-bin cover count of C3.
With empty selected tail the core contributions total nine. Since
Q(w)=2Q(u), this is `6Q(w)+9`. The deletion profile is d=3 and
K=`12*2^|u|`. The endpoint counts are S=1 and E=3.

The exit is not generally returnable to the one-apex spectrum at the same
target: C(OI)=OOOI has H9 and its exit has H29. The inherited strong-word
bound `H>=2^(m-2)+1` reduces any H29 core to m<=6; exhaustive words at those
lengths give none. Source/sink padding changes no H. This is an all-order
obstruction for this restricted target, not an obstruction to general
tournaments with H29.

## 5. A target-preserving three-arc rewrite, with exact eligibility

Now restrict the source word to

```
w = OOOO v,    v ends I.
```

Perform the two births of Section 4. Keep every vertex and every old arc.
Reverse exactly these three new arcs:

```
y->old0, y->old1, y->old2
     become
old0->y, old1->y, old2->y.                         (7)
```

The resulting graph has the same old apex and the transitive base order

```
old0, x, old1, old2, y, old3, old4, ... .
```

Its cut word is exactly `OIOOIO v`. This is a three-arc rewrite, not a claim
that each individual arc change preserves H, and not a native path-complement
move. Equivalently, delete only the newest y and reinsert it at this cut.

The two prefixes have states

```
OOOO   : (1,15,1,1)
OIOOIO : (15,45,1,1).
```

Equation (1) yields the exact difference

```
H(OIOOIO v)-3H(OOOO v)=14-2H(v)-2[v begins I].     (8)
```

Since H(v) is odd, equality with the standard target `3H+2` holds **if and
only if** v begins I and H(v)=5. Stripping leading I's and trailing O's, and
using the inherited finite length bound at H5, gives exactly the two reduced
words OOI and OII. Since v ends I, the complete eligibility language is

```
v=I^t z,    t>=1,    z in {OOI,OII}.              (9)
```

Thus (7) is a classified infinite count-neutral rewrite, with no parameter
fit or unbounded search assumption. The families are

| z | Source H | Returned H | First example |
|---|---|---|---|
| OOI | 75*2^t-9 | 225*2^t-25 | 141 -> 425 |
| OII | 105*2^t-9 | 315*2^t-25 | 201 -> 605 |

The 141 source here has nine vertices at t=1 and differs from the earlier
eight-vertex C4 source of the same count. The standard exit and rewritten
425 output have eleven vertices. No identification with that earlier
141/425 construction is assumed.

The graph changes from feedback vertex number two to one while preserving
the target and the complete old induced source. This is a physical return
to an earlier structural class, stronger than an equality between unrelated
integer expressions.

## 6. What the rewrite improves, and why it is not yet repeatable

Use the inherited positive future frame

```
q=(H,U,U-d,H+U-E)
```

for the rooted O/Z/J language from
[the dynamics theorem](../../h13_paired_descent_20261111/dynamics/PROOF.md).
Keep the **old apex** as root on both sides of (7). Formula (6) and (2)
give

```
q(standard exit)-q(rewritten exit)
  = (0, 2P(v)+2, 2P(v), 2P(v))
  = (0,14,12,12) for z=OOI,
    (0,10,8,8)  for z=OII.                       (10)
```

This is independent of t, since the leading I's are binary zeros. Every
common O/Z/J continuation therefore has no larger H after the rewrite,
by positivity of the inherited matrices. The current H is equal. This is
not a comparison under arbitrary root-amalgam contexts: the full seven
coordinates include endpoint and deletion changes that must be retained.

The returned graph has first base vertex of defect zero. Its second base
vertex, the retained x, has indegree one, but its defect is

```
65*2^t-4 for z=OOI;    91*2^t-4 for z=OII.        (11)
```

To prove (11), the second base vertex of an OI-start word has sole
in-neighbour its first vertex. Its defect is the number of paths starting
at that first vertex after deleting the second: this equals
`H(OOIO v)-1`. Since the prefix OOIO has state `(5,13,1,1)`, (1) gives
`H(OOIO v)-1=9+13Q(v)`. Now
`Q(I^t z)=2^t(Q(z)+1)-1`, with Q(z)=4 or 6.

All other base vertices have indegree at least three, as does the old apex.
The inherited lower bound on branching defect excludes their being 1- or
2-ready. Equation (11) excludes x, and the first base vertex has defect
zero. **The returned graph has no 1- or 2-ready root.** It cannot repeat the
standard small-defect ternary exit by changing the root.

The most direct fixed-suffix scalar repair also fails at the next step.
Multiplying the cut prefix P by three gives 15->45->135, but the prefix
Hamilton counts are 1->15->17. The second difference is
`17-3*15=-28`. Both prefixes still end O. By (1), reaching `3H+2` with an
unchanged suffix would require `H(v)+[v begins I]=-14`, impossible.
This stops that specific repeat mechanism, not every growing-defect or
multi-chart return.

There is a useful additional hostile bound. Delete old0 from the standard
exit, retain the old apex, and reinsert one new base vertex at any transitive
slot with either apex incidence. The remaining cut is `IIu`. Moving an
inserted O left only increases H, and moving an inserted I right only
increases H, by the inherited adjacent-inversion formula. Thus the maximum
is bounded by the two extremes

```
H(OIIu)=q+4Q+4,       H(IIuI)=q+P(u).
```

The first/last contribution in a strong word implies
`q>=Q/2+2` and `q>=P(u)/2+2`. Both extremes are strictly smaller than
`3(q+Q)+2`. So this entire natural old0-delete/reinsert attempt fails at
every size. The successful rewrite (7) changes the **newest** vertex's cut,
which is a different and precisely identified operation.

## 7. Reproduction and the next exact question

Run `python -B check.py` and `python -O -B check.py` in this directory.
`results.json` specifies the universe and retains literal rows, labels,
profiles, and positive/hostile controls. The checker uses fresh subset DP,
Kahn acyclicity, and degree-partition canonicalization; it imports no earlier
counting implementation. It checks every labelled tournament through order
six, all 511 binary cuts through length eight, all 15 OO-start strong cuts
through length six for the two-birth exit, the eight literal rewrite graphs
at t=1..4, and all eligibility suffixes through length nine. The infinite
claims rest on the proofs above.

The live next question is an exact compensation rule for the growing
indegree-one defect (11), or a new neutral rewrite changing the returned
prefix to a chart with a lawful ternary continuation. A fixed number of
response coordinates is available, but bounded registers do not bound their
integer values and do not establish coverage. This packet supplies the
source, target, physical map, preserved H, lost readiness, and necessary
sidecar for that next search.
