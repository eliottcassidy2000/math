# Unit joint entry: a complete grammar and Fibonacci enumeration

**PROVED, all finite orders.** Local theorem identifiers JE1--JE5 belong to
this file. Computational controls are separate evidence in `entry.json`.
This is a classification of one response stratum, not of all tournaments
or of the positive odd Hamilton spectrum.

Let H(T) count directed Hamilton vertex words, with H(empty)=1. Let
J_T(v,w) count ordered two-bin spanning path covers whose first bin is
nonempty, starts at v and ends at w. The second bin may be empty. The two
specified vertices may coincide. Every induced tournament has a Hamilton
path. Partitioning according to the first bin gives

    J_T(v,w) = sum_(simple directed paths P:v to w) H(T-V(P)).       (1)

For v=w the only simple path in (1) is the singleton. A tournament has
exactly one Hamilton path iff it is transitive: strong-component
decomposition and a Hamilton cycle in each nontrivial strong component
prove the nontrivial direction. These are the standard tournament facts
already used in [F1--F2, FUNNELS](../polarized_path_calculus_20261103/FUNNELS.md).

## JE1. Exact path-and-threshold grammar

For v=w, J_T(w,w)=1 iff T-w is transitive. For v!=w, J_T(v,w)=1 iff T
has the following presentation:

* A path P=(p_0=v,...,p_(L-1)=w), L>=2. Its consecutive arcs point forward
  and all its nonconsecutive arcs point backward.
* The remaining vertices have a transitive order z_1,...,z_m.
* Each z_i has a threshold k_i in {-1,0,...,L-1}: it beats precisely
  p_0,...,p_(k_i) and loses to the other path vertices.
* For every i<j, k_j-k_i<=1.

The presentation is unique: P is the unique simple directed v-to-w path,
and the residual transitive order is unique. It specifies every arc.

**Proof.** Every summand in (1) is positive. Hence J=1 means that there is
exactly one simple v-to-w path and its complement is transitive. A forward
nonconsecutive chord on P would give a shortcut. If an outside vertex loses
to p_a and beats p_b with a<b, it gives a shortcut through that vertex.
Thus its incidences must have the stated threshold form. If i<j and
k_j>=k_i+2, take a=max(0,k_i+1) and b=a+1. The legal path segment
p_a -> z_i -> z_j -> p_b replaces the corresponding segment of P, another
contradiction.

Conversely suppose the four conditions hold. A path first departing from
the P prefix cannot use a backward P chord, since that would revisit a
vertex. It must leave p_a through some z_i, so k_i<a. If it returns to an
unvisited p_b, then b>a and the last outside vertex z_j has k_j>=b. A
single outside vertex cannot do both; multiple outside vertices occur in
their transitive order, so i<j and k_j-k_i>=2, also impossible. It must
return to P to reach w. Thus P is the only path, and its transitive
complement contributes one to (1). QED.

**Algorithm and boundary.** Find a shortest v-to-w path by BFS, check its
chords, order the complement, read its thresholds, and check every rise.
This takes polynomial time (O(n^3) is a sufficient bound), without counting
all Hamilton paths. `entry.py` implements it. Failure does not imply J=0:
the next positive value may be greater than one. For instance the strong
four-vertex example with rows (8,1,3,6), w=1, has J(0,1)=2< H(T-w)=3.
The already rejected inequality J>=H(T-w) is not restored here.

## JE2. All-order class count and a weighted A000568 slice

Let a_n count isomorphism classes of triples (T,v,w) of order n with
J_T(v,w)=1, including v=w. Then

    A(x)=sum_(n>=1) a_n x^n = x/(1-3x+x^2),
    a_n=F_(2n),    a_1=1, a_2=3, a_n=3a_(n-1)-a_(n-2).          (2)

Here F_0=0,F_1=1. The first counts are 1,3,8,21,55,144,377,... . With
distinct ports only, subtract 2^(n-1). These are marked classes, not
Hamilton values, not ordinary unrooted class counts, and not a new range
of A000568 enumeration.

**Proof.** Diagonal ports contribute 2^(n-1): the transitive deletion has
one canonical order, and its arcs to w are arbitrary. For distinct ports,
fix path length L in vertices. Shift each threshold by one, obtaining a
word in {0,...,L} with no later rise of two. Record its running minimum r.
Below L there are two letters that retain r (r and r+1), and one letter
for each new lower minimum. At L there is only one retaining letter.
The word-count generating function from state r<L is therefore

    G_r(t)=(1-t)^r/(1-2t)^(r+1).

At initial state L it is

    G_L(t)=(1-t)^(L-1)/(1-2t)^L.                             (3)

These follow from (1-2t)G_r=1+t sum_(j<r)G_j, and
(1-t)G_L=1+t sum_(j<L)G_j. Hence the distinct-port series is

    sum_(L>=2) x^L G_L(x)
      = x^2(1-x)/[(1-2x)(1-3x+x^2)].

Adding x/(1-2x), the diagonal contribution, gives (2). JE1 guarantees
that these presentations neither omit nor duplicate marked classes. QED.

An automorphism fixing both ports fixes P pointwise and the ordered
complement pointwise; in the diagonal case it fixes the transitive
deletion. Thus all these marked structures are rigid. Writing j_1(T) for
the number of ordered pairs with J=1 gives the exact identity

    sum_[T of order n] j_1(T)/|Aut(T)| = F_(2n).              (4)

The number of labelled marked objects is n! F_(2n). This is a weighted
slice of the unrooted tournament class enumeration. It has the same
continuant denominator as the inherited unique-start funnel grammar:
the latter's strong rooted classes at order n+2 also number F_(2n).
Both proofs reduce to ordered transitive attachments and entry restrictions.
An operation-preserving bijection between these two families is not claimed.

For this stratum alone, the path/threshold code reconstructs the full
marked class. Ordering the diagonal bit words first, then L and the legal
threshold words lexicographically, gives an enumerative code of a_n
states, with fixed-length bit bound ceil(log2 a_n). This is minimal for
distinguishing all these marked classes at fixed n. Physical names require
their own transport sidecar. Closure under actual complements is not asserted.

## JE3. The unit dial's reentry condition, at every size

Let T be nontrivial strong and s_T(w)=1. Then

    (exists v, J_T(v,w)=1)  iff  H(T-w)=1.                  (5)

**Proof.** The reverse direction uses v=w. For the other direction use
the canonical deletion grammar from F1: T-w=U -> TT_t, where U is the
first strong component, w beats precisely one a in U, and s_U(a)=1.
Put q=H(T-w)=H(U). If q>1, U is a nontrivial strong tournament.

For v=w, J=q. For a tail vertex v, any simple v-to-w path stays in the
tail until w; its complement is U followed by a transitive sub-tail, so
every summand in (1) equals q. For v in U, start a directed Hamilton cycle
of U at v. Each prefix ending at a vertex other than a extends directly
to w. This gives |U|-1>=2 distinct simple paths, each with positive weight
in (1). Therefore J(v,w)>1 at every possible v. QED.

Now apply the inherited three-birth dial N^3 from [recovery NOTE, section3](../zero_ready_component_returns_20261103/recovery/NOTE.md).
It adds a,b,c with a->b->c->a, both b,c dominating T, and w->a as the only
old-to-new arc. For old v its new starting count is exactly J_T(v,w);
the new vertices have starting counts h,h+F-1,h+F, all at least three.
Consequently N^3 returns to a unique-start mark exactly when q=1.

If q=1, every possible unique-start output mark is an old vertex.
Deleting it still leaves the new directed triangle, so its deletion
count exceeds one. A second application of N^3 can be made from that
mark but has no unique-start output mark. This bounds an uninterrupted
sequence of legal unit-dial applications by two, even with arbitrary
mark changes between them. C3 attains two: after the first dial the old
mark has (s,q)=(1,3). There cannot be two successive *ready returns*.
Intervening source rewrites are outside this obstruction.

For the actual marked component underlying H170807, q=53. Thus JE3
explains, at all sizes, the earlier finite minimum-J=53 failure. It does
not classify all permissible repairs or rule out either ternary target.

## JE4. The only two-start escape is the marked four-vertex atom

Let T be nontrivial strong and s_T(w)=2. Then some J_T(v,w)=1 iff T has
four vertices. There is exactly one rooted isomorphism class in this
case: T-w is C3 on a->b->c->a, and w beats a,b and loses to c. Its
profile (h,s,q,F) is (5,2,3,6); every v in {a,b,c} has J(v,w)=1.

**Proof.** F1 again gives T-w=U->TT_t with U strong and
sum_(a in U,w->a) s_U(a)=2. Thus w has at most two out-neighbors in U.
If v is in U, Hamilton-cycle prefixes ending at its other vertices
give at least |U|-2 simple v-to-w paths. If |U|>=4 this is at least two.
The component U cannot be a singleton, and if |U|=3, w must beat exactly
two of its vertices. A nonempty tail also prevents J=1: besides a path
through U to its in-neighbor c of w, every v in U has the distinct path
v -> last_tail_vertex -> w. Tail vertices themselves have each path
weighted by H(U)>=3; v=w has that same deletion count. Hence U=C3 and
t=0 are necessary. In that four-vertex graph the only simple v-to-w
path follows C3 from v to c, then enters w. Its complement has at most
two vertices and is transitive, proving sufficiency. QED.

Its literal N^3 extension has H=31 and outgoing rows
(2,4,9,19,39,79,31). At old names0,1,2 its (h,s,q,F) profiles are
(31,1,17,163), (31,1,9,145), (31,1,5,87), respectively. Each is a valid
unit-start return, and each then meets JE3's obstruction to another unit
dial. These are actual incidences explaining a four-vertex/seven-vertex
transition. Equality of H=31 does not identify this graph with the older
31/95/97 source; the earlier source-reset packet retains its own graph
and observer. See [SOURCE_RESET](../source_reset_branches_20261103/README.md).

## JE5. Three starts are sufficient at every larger order, and minimal

For m>=1 let D_m consist of C3 on a->b->c->a, a transitive block Z of
size m, and w, with

    w -> Z -> C3,       w -> a,b,       c -> w.

This is strong and has order m+4. Deleting w leaves TT_m->C3, so q=3.
A Hamilton path starting at w must start with the whole Z block, then
traverse C3; otherwise it can never reach Z. Hence s_w=3. For each
v in C3 the only simple v-to-w path follows C3 to c and then w. Its
complement is a transitive block followed by at most two C3 vertices.
Thus J(v,w)=1 at all three core ports.

Consequently N^3(D_m,w) is a strong, source-preserving entry adapter:
it adds three vertices, retains all old arcs, and re-marks any of a,b,c
to obtain exactly one starting Hamilton path. Its inherited wrapper
counts and the N^3 identity give

    (h,s,q,F)(D_m,w)=(3+3*2^m, 3, 3, 3+3*2^m),
    H(N^3(D_m,w))=18*2^m+15.                              (9)

Indeed D_m is exactly Z_m(C3,c); the nonempty-prefix guard m>=1 is
essential. At m=0 the four-vertex boundary has H=5, not6, as in JE4.

Combining JE3 and JE4 proves a sharp threshold: among strong tournaments
of any specified order n>=5 with H(T-w)>1 and some J(v,w)=1, the minimum
possible starting count s_T(w) is exactly three. The minimum at order4
is two; order3 has no nontransitive deletion. This is an all-size
classification of the minimum entry budget, not just examples at larger n.

The adapter is uniform but not a closed ternary return: its restored
mark again has nontransitive deletion, and its scalar output is (9),
not a claimed 3H+2 or 3H+4. The useful new search target is a source
rewrite into such a three-entry interface with the required outer
response equation. Full joint responses and F49/F343 charts for m=0..6
are retained in `repair.json`; all-size claims follow from the proof.
