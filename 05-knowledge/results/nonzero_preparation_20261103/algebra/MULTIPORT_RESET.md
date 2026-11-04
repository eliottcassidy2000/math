# Exact multi-exit reset and its scalar growth obstruction

**PROVED — [independent audit](MULTIPORT_AUDIT.md) passed.** The proposed larger
two-birth operation is exactly evaluated below. For any tournament with H>=5,
using two or more old exits overshoots both `3H+2` and `3H+4`, even when the
root defect is allowed to change. Exact positive H offset with F tripling
also forces a singleton exit set at every order, including weak inputs.
A literal C3→13 boundary example genuinely supplies two scalar next letters;
it does not give a reusable larger-count class.

## 1. Object, inheritance, and operation

The closest mechanism is the single-exit C2 operator in
[construction/RESULT.md](../../multiport_defect_grammar_20261103/construction/RESULT.md).
The new variation makes the first newborn's old exit set arbitrary, retaining
the full incidence set rather than identifying it with one representative
port. The canonical hostile is loss of F tripling despite a scalar hit. The
corrected near miss is the genuine small positive C3→13: its two newborns
are scalar-ready, but its root debt is nonzero, and the same enlarged menu
cannot continue after the next reset. The least-used sidecar is the actual
position of each exit on a Hamilton path or in the two ordered cover slots.

The live concepts are: named exit sets, paired cover slots, scalar growth,
root debt, weak SCC order, and inherited reset-image nonreadiness.
All statements concern this fixed two-birth recipe. They do not exclude
larger surgeries, old-arc changes, deletions, or other ways to transport a
changing inner defect.

Let T be any nonempty tournament on V, r an old retained root, and B⊆V.
Birth x with x→b exactly for b∈B; every old vertex outside B beats x.
Then birth y with every old vertex beating y and y→x. Denote the result T'.
All old names and old arcs are retained; x,y are permanent fresh names.

If T is nontrivial strong, the first birth is strong iff B is a nonempty
proper subset. The final graph is strong iff B is nonempty: its cycle of
access is old T→y→x→B→old T. If B=V the first graph is weak but the final
one is strong. No strongness hypothesis is needed for the formulas or the
two lower bounds.

## 2. Full exact operator, including both correction terms

Use the square-zero vertex ring on V and full coefficient `[ ]`. P contains
the empty path, W counts nonempty paths ending r, and put

    A=sum_(b in B) A_b, E=sum_(b in B) E_b,
    Q=sum_(b in B) Q_br.

A_b counts paths starting b, E_b paths ending b, Q_br paths starting b and
ending r. These sums count each permitted path once, since its endpoint is
specified. Empty sums are zero. Unlike the singleton case, EA need not vanish.

The complete polynomials of T' are

    P'=P+z_x(P−E+PA−EA)+z_y P+z_x z_y(P+2PA−EA),
    W'=W+z_x(P−E)Q+z_x z_y P Q.                     (1)

To see the first x term, a path ending at x has old prefix ending outside B;
a path continuing after x has a suffix starting in B and its prefix, if
nonempty, must also end outside B. This gives P−E+(P−E)A. If both newborns
occur, either y ends the path, giving PA−EA, or y is followed immediately
by x, giving P+PA. An old-root-ending path containing both must use the
second alternative and continue from x to r. This proves (1) in every
empty-slot and root-in-B case.

Writing H=[P] and F=[WP], full extraction gives

    H'=H+2[PA]−[EA],
    F'=F+2([WPA]+[P²Q])−[WEA]−[PEQ].                (2)

Both cross-port corrections are necessary. Setting |B|=1 annihilates EA,
WEA and EQ by repeated use of its port and recovers the inherited C2 law.

The newborn start polynomials are

    A'_x=z_x(1+A)+z_x z_y A,
    A'_y=z_y+z_x z_y(1+A).

Consequently their exact scalar start gaps are

    sigma_x=[EA], sigma_y=H−e_B(T),                  (3)

where e_B(T)=[E]. Formula (3) shows that the generalized operation can
create scalar readiness; the obstruction below is its growth at the
required scale, not a blanket assertion that its newborns are never ready.

## 3. An all-order scalar lower bound

For nonempty B with m=|B|,

    H(T') >= (m+2) H(T).                             (4)

Proof. For each old Hamilton path P, construct these m+2 new paths:

* append y,x to P;
* insert x immediately before P's first B vertex, then append y;
* for each b∈B, insert the consecutive block y,x immediately before b.

The second path is legal because the predecessor of the first B vertex,
if any, is outside B and therefore beats x. In the last family, any old
predecessor beats y, and y→x→b. The first path ends at x, the second at y,
and the last m paths end at an old vertex; different b give different
block locations. Thus all m+2 outputs are distinct. Deleting x,y recovers
the old Hamilton path, so outputs from different old paths are also distinct.
This proves (4) without a strongness or an order restriction.

If H>=5 and m>=2, (4) gives H'−3H>=H>=5. Hence neither desired offset2 nor4
is possible. This conclusion permits arbitrary changes of F and of the root
defect; it is strictly stronger than an exact-F obstruction at these inputs.
In particular it applies to every actual prepared H1429 source and every
larger descendant, including weak variants of those counts.

For B empty, the actual final order is T→y→x, so H'=H and F'=F. That case
cannot supply a positive ternary digit either.

## 4. The rooted-cover lower bound and exact-tripling boundary

Let C_(B,r) count old F-covers `(W,V)` in which **both** paths meet B.
Here W is nonempty ending r; V may ordinarily be empty, but is nonempty
in this subcount. For nonempty B,

    F' >= 3F+2 C_(B,r).                              (5)

For every old cover, append y,x to its free path V. This is the first
lift. If V meets B, insert x before its first B vertex and append y to V
for a second lift; insert the block y,x there for a third. If V does not
meet B, W must meet B; instead insert x before its first B vertex and
append y to V, or insert y,x there and leave V unchanged.

If both slots meet B, the chosen second and third lifts used V. The same
two choices in W give two additional lifts. Their slots and the positions
of x,y distinguish all of these variants. Deleting x,y recovers the old
typed cover, including the empty V case. The legality checks are exactly
the first-B and block checks of §3, and W still ends at old r. This proves
the injection and (5).

If m>=2 and some b∈B reaches r, take a shortest directed path from the set
B to r. It contains exactly one B vertex: a later B vertex would give a
shorter such path. If r∈B, use its singleton path. Use this path as W and
any Hamilton path of its complement as V. The complement contains another
B vertex, so C_(B,r)>=1 and

    F'>=3F+2.                                       (6)

This argument does not require strongness; in a strong graph reachability
is automatic. It cannot simply be extrapolated to an unreachable root.

If no vertex of B reaches r, every B vertex lies strictly after r's SCC
in the total SCC order. Write T=A→D, where A is the prefix through r's
SCC and B⊆D. The two newborns also lie after A, so T'=A→D'. Direct ordered
path concatenation gives

    H=H(A)H(D), H'=H(A)H(D'),
    F_r=F_r(A)H(D), F'_r=F_r(A)H(D').                (7)

The first cover path cannot enter D and end at r; all of D belongs to the
free path. Since F_r(A)>0, equality F'=3F in (7) forces H'=3H, with zero
scalar offset.

Combining the empty case, (6), and (7), **for every tournament a positive
offset H'>3H with F'=3F forces B to be a singleton**. It is then exactly
the inherited C2 operation. By the proved complete
[joint classification](JOINT_CLASSIFICATION.md), its positive guard is
precisely an old auxiliary with sole exit r and positive e_r(T−a).
The reduction to a singleton is proved here independently; this last
identification retains its separate audit dependency.

## 5. Sharp small positives and their stopping boundary

Take T=(2,4,1), a directed triangle, and B={0,1}. Both births are strong.
The actual output is

    T'=(18,20,25,3,8), H=3, H'=13=3H+4.

At old roots0,1,2 the input F values are all3 and the output values are
respectively18,15,14. Thus this is a genuine scalar +4 with a nonzero
retained-root debt. Its five actual start gaps are

    (5,5,13,2,1).

The two newborns therefore supply both scalar next letters by singleton
C2 resets: auxiliary x=3 gives H43, and auxiliary y=4 gives H41. These
are literal strong outputs and the receipt retains their full matrices,
old-root responses and all start gaps. At old root0, the branches are
`(H,F)=(13,18)→(43,56)` and `(13,18)→(41,58)`; the root defects change from5
to13 or17. Their existence is not discarded
because the original F did not triple.

This small fork still does not close the reset menu. Each branch is now
a C2 image of a strong five-vertex input. The inherited all-auxiliary
reset-image theorem (construction/RESULT §4, routed to READY §4) excludes
scalar gaps1 and2 at every auxiliary. At H41 or43, a non-singleton reset
is separately excluded by (4), and the empty reset does not make a
ternary step. Thus another nonzero step requires an operation outside
this enlarged reset menu. This is a mechanism-level stopping result,
not a spectrum obstruction.

The H>=5 hypothesis of the scalar obstruction cannot be dropped. The
triangle above supplies H3; the transitive two-vertex input `(2,0)` with
B={0,1} supplies another literal boundary H1→H5=3H1+2. Neither is evidence
for a uniformly reusable operation at larger H. No global minimality or
classification of the H1 boundary is claimed.

## 6. Exact finite universes and reproducibility

The standalone source is
`04-computation/nonzero_preparation_20261103/algebra_multiport_reset.py`.
It reads no producer code or response receipt. It enumerates every labelled
input through order4, every B including empty and full, and every retained
old root. Directed simple paths are counted by literal permutations on the
source and actual two-birth target. It checks the full formulas (2)–(3),
every old arc, both strongness boundaries, and literal distinct image paths
and covers for the two injection proofs. It also checks the unreachable
SCC-factorization ratio and the all-tournament exact-tripling conclusion.

The separate cheap scalar probe covers **every labelled strong tournament
of orders3–5** and every exit set of size at least2, without a port or
isomorphism filter. It finds only the six labelled C3/pair choices for +4;
no +2, and no order4/5 hit. That finite failure does not prove the all-order
obstruction; the injection (4) does. The explicit H13 fork and H1 boundary
are counted again from their actual adjacency matrices.

```text
python3 -B 04-computation/nonzero_preparation_20261103/algebra_multiport_reset.py --output /tmp/algebra_multiport_reset.normal.json
python3 -B -O 04-computation/nonzero_preparation_20261103/algebra_multiport_reset.py --output /tmp/algebra_multiport_reset.optimized.json
```

Both normal and optimized runs pass **296,258 always-active gates** and
produce identical 2,958-byte receipts. The operator census contains 1,098
graph/exit cuts and 4,306 old-root marked cuts. The source checks 12,239
distinct Hamilton lift images and 79,757 rooted-cover lift images. The full
strong3/4/5 scalar probes have 8/264/14,144 multi-exit choices respectively.

Frozen source SHA256:
`b91ca738b48b3fe1aebb3cd11b15b3b52a7ba9b9cd6d5fbf0b6edf1e7d468d08`.
Frozen receipt SHA256:
`e3cdcb0cf11d288cd0361eafe0d8913006b66226286c4401cb90b7044fa4a30a`.
This new operator and its quantified stopping bounds passed independent analytic/source/replay audit. All earlier frozen artifacts are
unchanged.
