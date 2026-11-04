# Moving the feedback mark repairs an infinite blocked ternary branch

**PROVED scoped constructions and obstructions; FINITE-EXACT controls.**
Full positive odd spectrum coverage, including both prime lanes, remains
OPEN. This note supplies an actual infinite source repair and reusable
balanced children. It does not claim a target selector for arbitrary primes.
The numerical plus-two and15*2^m-1 families were already attainable by
inherited constructions. The additions concern the operation obstruction,
physical mark transport, and the explicit balanced attachment interface.

## Inheritance and the question changed by the new source

The closest proved mechanism is the marked clone formula
[CL1](../../return6325_completion_20261103/CLONE_LAW.md), together with
the rooted defect and two-birth reset in
[ADAPTIVE sections 2 and 3](../../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md).
The hostile is a long uniform clone block. The corrected near miss is
[STABLE S6](../../return18979_palindrome_20261103/STABLE.md): arbitrary
intermediate repairs cannot overcome a fixed exterior and retained clone
mark at a ready endpoint. The least-used sidecar here is the physical
choice of which old clone becomes a feedback vertex, with the other old
apex retained by name.

The live board compares (i) exact scalar target selection, (ii) the two-apex
word representation, (iii) positive clone-run coefficients, (iv) actual
feedback re-marking, and (v) balanced two-port reentry. A fixed word suffix
fails on the source below. Re-marking changes the prefix and repairs one
whole ternary branch; the repaired child then supplies a different prime
and seven-prime family through an old named root.

## 1 A positive cone for a unique-entry clone module

Let E be a finite tournament and y a vertex with at least one in-neighbor.
Use the actual rooted profile

    h=H(E), s=start_y(E), q=H(E-y), F=E_y(2),
    delta=F-h-q+s.

F counts two ordered spanning path bins, first nonempty ending at y,
second possibly empty. Adjoin an ascending clone module U of size m>=1,
with y->U and U->E-y. The inherited CL1 formula gives

    H(E[U])=A*2^m+B,
    A=F-q=h-s+delta,       B=q-delta.                 (C1)

Here the notation E[U] denotes this adjoining operation, not replacing y.
Every Hamilton path has at most two U-runs: one can start in U, and only
y can supply another entry. Two nonempty runs contribute 2^m-2; their
exterior cover count is F-q. This is precisely CL1's bijection.

**Cone theorem.** For every such actual E,y,

    A>=1,       2-A<=B<=A,       B odd.              (C2)

Proof. The inherited weighted-cut law gives delta>=0. For every Hamilton
path of E-y, insert y immediately after its last in-neighbor on that path.
The result is a path of E not starting at y. Restriction recovers the old
path, so h-s>=q>=1. Thus A+B=h-s+q>=2 and
A-B=h-s-q+2delta>=0. Finally q is odd and delta is even by the weighted-cut
law. QED. This derives the cone from positive actual path counts; it is not
an assumption about arbitrary integer profiles.

Consequently, for M=2^m, every attainable count in this operation lies in

    union_(A>=1) [A(M-1)+2, A(M+1)] intersect positive odds.   (C3)

When M>2A+1 this leaves a genuine interval before the next coefficient
layer. The bound is necessary, not a complete classification of pairs A,B.
For example (A,B)=(4,-1) satisfies the cone but is impossible.

**The missing pair.** If A=4,B=-1, then q-delta=-1 and
h-s+delta=4. Since h-s>=q, one obtains q=1, delta=2, h-s=2.
The deletion E-y is transitive, so s is zero or one. Oddness of h forces
h=3,s=1. The weighted-cut formula sums over ordered pairs of distinct
in-neighbors, so a positive defect forces y to have at least two
in-neighbors a,b. In particular indegree one has delta=0, also the direct
newborn N identity in ADAPTIVE. Since s=1, y beats the first vertex u of the deletion chain. The
induced tournament on y,u,a,b is strong, with its chain ordered u,a,b;
it has five Hamilton paths. Insertion monotonicity then gives h>=5,
a contradiction. This small nonrealizability witness is the extra
information needed beyond the numerical cone.

## 2 An infinite source on which every fixed-mark suffix fails

Let T_m have old names x,y,c_1,...,c_m and arcs

    x->y->c_1->...->c_m->x,
    y->c_i, c_i->x, c_i->c_j for i<j.

It is strong and H(T_m)=2^m+1. In the standard two-apex word convention,
feedback edge x->y and bit one meaning an apex beats a chain vertex,
its word is 2^m. A suffix v adds new vertices after every old chain vertex,
with arbitrary letters in {0,1,2,3}; the two apices remain x,y.
The old c_i still form a module with sole incoming exterior vertex y,
so C1--C2 apply to every suffix, however long.

**All-suffix obstruction.** For every m>=3, no such suffix has either

    H=3H(T_m)+2 or H=3H(T_m)+4.                      (C4)

For the first target, C1 would give A M+B=3M+5. If A<=3, B<=A
makes the left side at most3M+3. If A>=4, B>=2-A forces
(A-3)(M-1)<=6, impossible for M>=8. For the second target use3M+7:
the same argument excludes M>=16. At M=8 the sole remaining cone
solution is A=4,B=-1, excluded above. This proves C4 for arbitrary
suffix length, not merely for a finite tested bank.

The boundary is sharp for both branches. For m=1, suffixes20 and12 give
11 and13. For m=2, suffixes22 and02 give17 and19. These are literal
two-apex words, not abstract numerical profiles. A finite adaptive bank
of fixed-mark suffixes therefore cannot repair C4 either.

## 3 Re-mark one old vertex and recover the entire plus-two branch

For m>=2, select the ordered feedback pair (c_m,x), retaining the old
names. The remaining old tournament is the transitive chain
y,c_1,...,c_(m-1). Its actual word is

    2 0^(m-1).

Append two new chain vertices a,b with letters32. Equivalently:

* every vertex in the remaining old chain beats a and b;
* c_m->a, x->a, x->b, b->c_m, and a->b.

All old arcs and vertices survive. Call the completed tournament B_m.
It has order m+4, is strong, and

    H(B_m)=3*2^m+5=3H(T_m)+2.                        (C5)

Strongness follows because a reaches the old strong source through b->c_m,
while the old source reaches both new vertices. The first intermediate
birth is a sink and is not claimed strong.

For a direct all-size count, remove c_1,...,c_(m-1). The five-vertex
exterior is the fixed two-apex word232, rooted at its first chain
vertex y, with literal profile

    (h,s,q,F)=(9,3,5,11).

Each removed c_i loses only to y in this exterior. Apply C1 to the
module of size m-1 to get (F-q)2^(m-1)+(q-delta)
=6*2^(m-1)+5. The four fixed profile counts are independently checked
by permutation enumeration. This proves C5 without extrapolation.

The m=1 boundary must be separate: word232 has H9, not11. Use the
explicit suffix20 on word2 instead. The same actual selection rule with
that base case therefore gives the plus-two branch for every m>=1.

The numerical existence of this branch is already a consequence of the
inherited ready-root theorem: for m>=2, c_1 is one-ready in T_m, and c_m
is one-ready in the converse. The final C5 gadget is exactly the dual
two-birth ready reset at c_m, with the two fresh names exchanged. What
this section adds is its explicit compilation as two-apex suffix32 after
physical mark replacement, paired with C4's universal obstruction for
the original fixed marks. It does not claim a previously unknown family
of attainable values.

This is a precise operation-level distinction. The map from the original
marked source to B_m preserves its complete induced old tournament and
the old vertex names. It changes which old vertex is an apex of the
word representation. Forgetting that mark replacement loses the
construction: C4 proves that the original marks cannot perform it.

## 4 The retained old root supplies primitive seven-prime children

At the old physical vertex x in B_m, the exact two-slot count is

    E_x(2)=6*2^m-3,                 m>=2.            (C6)

Here is a short finite skeleton proof valid at every m. Put k=m-1 and
use the same five-vertex exterior. A two-bin path cover with first bin
ending at x has at most three clone runs: at most two starts of bins,
and at most one later entry through y. Grouping by r nonempty runs gives

    E_x(2)=d_1 Surj(k,1)+d_2 Surj(k,2)+d_3 Surj(k,3),

with nonnegative integer exterior skeleton counts d_r and
Surj(k,r)=r! S(k,r). The direct finite exterior counts at k=1,2,3
are21,45,93. Since these surjection rows are (1,0,0),(1,2,0),(1,6,6),
they determine (d_1,d_2,d_3)=(21,12,0). Thus
E_x(2)=21+12(2^k-2)=12*2^k-3, proving C6. The certificate enumerates
all vertex permutations with the cut immediately after x as an
independent count of each of the three required bin totals.

Now use the inherited two-birth reset at x: add a vertex beaten only by x,
then one beaten only by that newborn, retaining the first newborn as root.
Its balanced profile is (b,t,t,b), where

    b=H(B_m)+2E_x(2)=15*2^m-1,
    t=H(B_m)=3*2^m+5.                                (C7)

The final root and its child are named new vertices. All old names survive;
the reset intermediates are strong because the input is strong.
The pair is primitive exactly when m is not11 modulo12, since

    gcd(b,t)=gcd(t,26),
    13|t iff 2^m=7 (mod13) iff m=11 (mod12).

In particular m=3j, j>=1, always gives a primitive balanced pair with

    b=7 q_j,       q_j=(15*8^j-1)/7,
    q_1=17,        q_(j+1)=8q_j+1.                  (C8)

The first supplied quotients17,137,1097 are prime, giving actual strong
counts119,959,7679 at orders9,12,15. More generally, b is divisible by7
exactly when m=0 modulo3. No assertion of infinitely many prime q_j is
made. For reference, b is divisible by49 in this subfamily exactly when
j=5 modulo7, by expanding (1+14)(1+7)^j-1 modulo49.

The new child is a reusable balanced input to the inherited R/N algebra;
arbitrary primitive pairs are not automatically in a separately guarded
eight-state controller. The family itself has O(log X) outputs below X.
Neither this branch nor its scalar coprimality proves full prime/7p coverage.

The scalar family15*2^m-1 itself is also inherited:
[R2 of the prime-lane packet](../../prime_lane_ports_five_groups_20261102/realization.md)
already constructs it as V_(m,2), including119=7*17. C7 adds a specific
source-preserving route from T_m through B_m, the named balanced pair
(15*2^m-1,3*2^m+5), and its exact primitive selector. It is not an extension
of the previously known numerical spectrum.

Incoming [Session175's fixed clone library theorem](../../session175_tournament_spectrum_20261104/clone_obstruction/PROOF.md)
supplies a much broader prime-density obstruction, consistent with this
limitation. Its [asymmetric cyclic family](../../session175_tournament_spectrum_20261104/construction/README.md)
has counts27*2^r-17 and a different7-adic clock. Its arbitrary cyclic-module
kernel concerns uniform substitution. C5 instead singles out an old clone
as a feedback vertex before the new cuts, preserving every old arc but
splitting the old module's new exterior incidence. These are distinct
literal operations, even though both give exponential families.

## 5 The unclosed branch and exact controls

All possible feedback pairs of T_m, for m>=3, are the original pair or
one old apex and one clone. Deleting two clones leaves a directed triangle.
Their word representatives are

    2^m,
    2 0^i 1^(m-1-i),
    1^i 3^(m-1-i) 2,                 0<=i<m.         (C9)

The bit3 in the final family is essential: both marked apices beat each
later clone. The complete marked-pair generator in the certificate checks
C9 directly from the old relation and handles m=1,2 separately.

**FINITE-EXACT:** for every m=3,...,10, no re-marking in C9 followed by any
suffix realizes3H(T_m)+4. This is an all-suffix result for each of these
eight supplied sources: trailing3s are sinks and can be deleted; a remaining
suffix ends in a non3 and produces a strong tournament, so the cited
[Moon--Busch strong-order floor, Corollary1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf)
bounds its length. The maximum required
suffix lengths are2 for m=3,...,6 and3 for m=7,...,10. All such words are
enumerated. No all-m moving-mark obstruction is inferred.

The cheapest next obligation is consequently an operation that changes
the ordered chain relation or exits the two-apex word class on at least
the m=3 source, while retaining a sufficient return interface. A larger
finite menu of fixed-mark suffixes cannot help C4; merely allowing all
old feedback re-markings still fails at m=3 by the complete finite bound.

Run `python3 -B check.py` or `python3 -B -O check.py` in this directory.
Both must reproduce `check.json`. The script uses the inherited standard
library two-apex transfer, independently implements labelled subset DP,
and counts small examples directly by vertex permutations. No ordinary
assertions carry validation. Exact universes and check counts are recorded
in the JSON receipt. No new A000568 coefficient is claimed by this note.

