# The Hamilton-support image needs a realization guard

**PROVED I1–I5, independently audited [PASS](../audit_image/AUDIT.md).
FINITE-EXACT receipts.** Local labels I1–I5 belong
to this note. Finite results are controls, not the all-order proof.

Closest mechanisms: [two-step reconstruction G1](../../guard_reconstruction_20261026/README.md),
[elementary endpoint reconstruction THM-4807](../../../../01-canon/theorems/THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation.md),
and [joint rooted composition R1–R2](../../rooted_lines_frobenius_sheets_20261026/README.md).
The corrected near miss is to use a promised-image decoder as an image test.
The canonical hostiles are the transitive triple, the first strong unused
arc at five, and the new five-vertex false image below. The neglected
sidecar is an ordered spanning path cover with joint endpoints.

## I1. An exact image criterion with a positive defect polynomial

Let D be a loopless oriented graph on a fixed finite vertex set. U(G) is the
union of directed Hamilton-path arcs of G, for any digraph G; H(G) counts
its Hamilton vertex words. Include the empty graph with H=1 and U empty.

Run this partial decoder on D:

1. Its strong components must be totally ordered by reachability. Between
   consecutive components D must contain every forward arc; between
   nonconsecutive components it must contain none.
2. Keep its existing arcs and complete all cross-component pairs forward.
3. For each missing pair x,y inside a strong component, exactly one direction
   must have at least two D-paths of length two. Orient it that way. Reject
   when neither direction, or both, qualify.

Call the resulting tournament C(D), when all tests pass. This is the
inherited O(n³) decoder with its necessary support conditions made explicit;
it is not yet a recognition algorithm. Put M=E(C(D))\E(D), and define

    F_D(z) = sum over Hamilton words P of C(D) of z^|E(P)∩M|.

All coefficients are nonnegative integers. In particular F_D(0)=H(D),
F_D(1)=H(C(D)).

**I1 (PROVED).** D belongs to the Hamilton-support image of
tournaments if and only if:

    C(D) is defined;   U(D)=D;   F_D(1)=F_D(0).               (1)

Equivalently, every positive-degree coefficient of F_D vanishes.

**Proof.** If D=U(T), reconstruction gives C(D)=T. Every Hamilton word of T
lies in D and conversely, so U(D)=D and F_D is constant. Conversely D is a
subgraph of T=C(D). Equality of the counts, together with inclusion of the
Hamilton languages, forces their equality. Hence U(T)=U(D)=D. No signed
cancellation is used. The empty and singleton cases satisfy the same rule.
The uniqueness of T is inherited reconstruction. QED.

This is an exact all-order image characterization with an explicit
exponential finite evaluator, not a small local list of forbidden patterns.
It isolates the missing realization condition rather than assuming it.
The standard subset DP retains, for each visited set S and last vertex v,
a polynomial A(S,v); its update along v→w in C(D) multiplies by z if the
arc is missing from D and by1 otherwise. Starting from A({v},v)=1 and
summing at the full set evaluates F_D. U(D) is obtained by endpoint subset
DP or literal paths. Polynomial degree is at most n−1. Ordinary Hamilton
counting by this recurrence is classical; no improved complexity is claimed.

## I2. The defect is an ordered joint-path-cover transfer

For a nonempty vertex set, replace z by a variable z_e for each e in M. Split any
Hamilton word of C(D) at its M-edges. It becomes a unique ordered list
(P1,...,Pr) of nonempty directed paths in D, with pairwise disjoint vertex
sets covering all vertices. Consecutive blocks are joined by missing arcs
from the exit of Pi to the entrance of P_(i+1). Conversely every such list
with those actual missing arcs concatenates to exactly one Hamilton word.

Thus the coefficient of total degree r−1 is exactly the number of those
jointly compatible ordered r-path covers. This bijection proves I2. It
retains vertex incidence, all entrance/exit labels, block order and the
connector realization predicate. Multiplying separate endpoint totals
does not implement this map.

The empty vertex set contributes the unit1 separately.

The multivariate version also records which missing arcs cooperate. For
the directed four-cycle D=(2,4,8,1), complete it to
T=(6,12,8,1) by adding 0→2 and1→3. Every D-arc occurs in a Hamilton word;
adding either chord alone creates no new Hamilton word. Adding both does:

    F_D(z_02,z_13)=4+z_02*z_13,
    new word: (1,3,0,2).

This example is a general composition hostile, not an accepted local
decoder candidate: its missing pairs have one two-step route in each
direction. It shows why individual connector tests cannot replace joint
realization in an arbitrary partial graph.

## I3. The first local false image, and nonheredity

Exhausting all loopless oriented graphs through order5 gives the first
false image accepted by the local decoder and U(D)=D at order5:

    D=(10,20,16,6,1),   C(D)=(14,20,16,22,1).

The missing arcs are0→2 and3→4. D has seven Hamilton words; its completion
has nine, and F_D(z)=7+2z. Both new arcs become visible. Thus local decoding
and path-saturation alone do not characterize the image. Exactly480 of
1504 path-saturated decoder candidates at order5 fail the realization test;
the remaining1024 are precisely the distinct supports of all1024 labelled
tournaments. These counts are FINITE-EXACT.

The image is not hereditary under induced vertex deletion. The support
of the transitive triple is the directed path0→1→2. Deleting1 gives an
edgeless two-vertex graph, which is not a tournament support. Therefore
there is no characterization of this full image class by forbidden induced
subdigraphs. Its correct deletion operation is

    D \Downarrow v = U(C(D)−v),                              (2)

not induced deletion D−v. This operation commutes with relabelling and
successive deletions, because C and U are inverse on the image.

Insertion is not monotone in old visible arcs either. Inserting the middle
vertex of a transitive triple hides an old arc; extending that triple to
the strong tournament (6,12,8,1) makes the previously hidden0→2 visible.
An image recurrence must decode, attach, and recompute support, or retain
equivalent spanning-incidence information.

## I4. Exact complement-line incidence from two supports and a transport

Let D,E be valid supports on the same named set, and decode T=C(D), S=C(E).
For n≥3 they are the two ends of one literal complement line in some path
order precisely when their common directed arcs form one spanning directed
path P:

    E(T)∩E(S)=E(P).                                         (3)

Indeed K_P(T) keeps the path arcs and reverses every other pair, so (3)
is necessary and sufficient. The agreement graph has at most one candidate
path. Testing (3) after decoding costs O(n²); total decoding cost is O(n³).
Its blue guard is that the positional reversal of P is an anti-automorphism
of T. Otherwise the line is black. The guard has the same value at S.

For ordinary isomorphism representatives T,S, let b(T,gS) indicate that
(3) holds and its path is black. Then the number of black Hamilton words
of T whose complement target is isomorphic to S is exactly

    L_black(T,S) = (1/|Aut(S)|) sum_(g in S_n) b(T,gS).       (4)

Every target isomorphic to S has exactly |Aut(S)| such transports, and the
agreement path is unique, proving (4). The total/blue versions use their
respective guards. In particular L_black(T,S)/|Aut(T)| is symmetric in T,S.
The converse-merged graph further identifies S with S^op, retaining loops
and the declared line convention. This is an exact incidence description,
not a proof that its black support graph is connected. Forgetting the
transport g destroys the common vertex alignment needed by (3).

## I5. A framed common hub and the remaining connectivity obligation

Fix a tournament h with distinct ports s,t, and two Hamilton words P,Q
from s to t. Set C_P=K_P(h). The inherited hinge formula gives

    K_(a,P,b)(1_a→C_P→1_b)=e*[h;s,t]*e.                    (5)

The same labelled hub occurs for Q. Each padded path has exactly the
blue/black status of its inner path, since reflection swaps the two
singleton markers and restricts to the inner positional reflection.
Consequently two black paths with the same joint ports give a two-black-line
walk between their padded nonstrong sources. The hub is strong: its path
a,P,b closes by b→a. Distinctness of the sources is not guaranteed.

This is a useful exact hub family, derived from existing rooted composition.
It does not connect different port pairs automatically. The prior
[terminal-basin reduction](../../black_descent_20261020c/scc.md) shows that
black connectivity is equivalent to connectivity after contracting descent
forests onto nonstrong active classes. Equations(3)–(5) provide their joint
transport and hub guards, but not an all-order connection between every
terminal basin. Pure-blue exceptions and the earlier fixed-boundary hostile
remain load-bearing. Black-core connectivity remains OPEN.
