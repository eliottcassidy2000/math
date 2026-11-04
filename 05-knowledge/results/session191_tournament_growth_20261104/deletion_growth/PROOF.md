# Count-neutral padding and the first failure of Fibonacci deletion growth

**Session 191, 4 November 2026. Status: proof and finite-exact computation independently audited.**

For a tournament T, H(T) counts directed Hamiltonian vertex words, with H(empty)=1.
For a marked vertex v put Delta(T,v)=H(T)-H(T-v). This packet concerns an induced
vertex deletion, not changing arcs on a fixed vertex set. The all-order
{7,21} spectrum remains open.

## 1. Inheritance and scope

The inherited identities used here are:

* THM-002, `01-canon/theorems/THM-002-ocf.md`: the Grinberg--Stanley odd-cycle
  collection formula H(T)=sum_Gamma 2^|Gamma|, over sets of vertex-disjoint
  directed odd cycles of length at least three. A directed cycle is counted
  modulo cyclic rotation, not merely by its vertex set.
* THM-070, `01-canon/theorems/THM-070-claim-a-from-ocf.md`: its marked-vertex
  difference formula. THM-4114,
  `01-canon/theorems/THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature.md`,
  already proves the more general nonnegativity of presence-cube differences.
* Moon's strong-subtournament theorem, sourced and stated in THM-4097,
  `01-canon/theorems/THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md`,
  section 1: every strong tournament of order m>=4 has at least two vertices
  whose individual deletion remains strong. Source: J. W. Moon, *On
  Subtournaments of a Tournament*, Canadian Mathematical Bulletin 9(3)
  (1966), 297--301, Theorem 2, https://doi.org/10.4153/CMB-1966-038-7.
* The complete strong order-four tournament is the unique class with H=5,
  with two vertices whose deletion is C3. The script rechecks this directly
  over all 64 labelled tournaments; it does not import a class census.

The pinned incoming THM-4872 was read at commit
3192fab75fd27f262157a67c9e2aa0c17c47a36b, file
`01-canon/theorems/THM-4872-unit-start-plus-minus-two-port-pairs-and-the-z3-trivial-symmetry-split-of-three-cycle-tournaments.md`.
Its Fibonacci numbers count certain rooted/few-cycle isomorphism classes;
its +/-2 operation changes old arcs. Those are different invariants and
operations from the deletion increment here. The file is pre-audit at that
checkpoint. This packet does not claim its results or confuse its +/-2
pairs with one-vertex insertions.

## 2. Exact marked difference and hereditary nonnegativity

Write U(Gamma) for the union of the cycles in an odd-cycle collection Gamma.
The OCF gives directly

    Delta(T,v) = sum_{Gamma in T, v in U(Gamma)} 2^|Gamma|.       (1)

Indeed deleting v removes precisely those cycle collections whose covered
set contains v. Equivalently, by selecting the unique cycle through v,

    Delta(T,v) = 2 sum_{odd directed C containing v} H(T-V(C)). (2)

Each remaining tournament in (2) includes the empty-tournament convention.
This is the inherited THM-070 identity, restated to fix normalization.

If S contains v, every collection in the right side of (1) for T[S] is also
one for T. Therefore

    0 <= Delta(T[S],v) <= Delta(T,v).                          (3)

All these increments are even. Formula (3) is monotonicity of the **marked
increment**, stronger than merely H(T[S])<=H(T), and is also an immediate
specialization of inherited THM-4114. No equality of full response profiles
is asserted.

## 3. Exactly which vertex growth is count-neutral

**Corollary.** H(T)=H(T-v) if and only if v is a singleton strongly connected
component of T. More generally, for S subset V(T),

    H(T[S])=H(T)

if and only if every vertex outside S is a singleton strong component of T.

**Proof.** A singleton strong component lies on no directed cycle, so (1)
vanishes. Conversely, a vertex in a nontrivial strong component lies on a
directed triangle. To see this, take a shortest directed cycle through v,
written v,x1,x2,...,x(k-1),v. If k>3, the chord v->x2 would shorten it;
hence x2->v, and v,x1,x2,v is already a triangle, a contradiction to k>3.
The triangle contributes at least 2 to (1), so Delta>0. If a removed set
contains such a v, ordinary induced-subtournament monotonicity gives
H(T[S])<=H(T-v)<H(T). If all removed vertices are singleton components,
no cycle collection is lost. QED.

Thus an induced enlargement at the same H retains every nontrivial strong
component intact; it only inserts transitive singleton components among
those components. Order can grow indefinitely at fixed H for that reason,
whereas a strong one-vertex enlargement always raises H.

## 4. Simultaneous strong deletion

Call (T,v) a **strong deletion pair** if T and T-v are both strong. The
smallest order is four. The root v remains a named vertex.

**Lemma (joint strong-deletion descent).** Every strong deletion pair (T,v)
of order n>=5 has some w!=v such that both T-w and T-{v,w} are strong.

**Proof.** Let U=T-v, of order m=n-1>=4. Moon's theorem gives at least two
vertices w for which U-w is strong. The root v has at least one in-neighbor
and at least one out-neighbor in U, because T is strong and U is strong.
Deleting w can make v fail to have both kinds of neighbor only when w is
its sole in-neighbor or its sole out-neighbor. At most one of those two
neighbor classes is a singleton: they partition m>=4 vertices, so both
cannot be singletons. Hence at most one of the at least two candidate
vertices w is forbidden. Choose another. Then U-w is strong and v still
has both an in- and an out-neighbor in it, which makes (U-w)+v=T-w strong.
QED.

Together with (3), this proves that every strong deletion pair with
Delta<=B has a nested chain of such pairs down to order four, all with
Delta<=B. This is a complete construction theorem, not representative
selection by the scalar H or by Delta.

## 5. Exact bounded-increment census

Root label 0 is fixed. Start with all 64 labelled four-vertex tournaments,
retain strong deletion pairs and quotient by root-preserving isomorphism.
There are exactly two rooted classes, both with H=5 and H(T-0)=3.

For each retained rooted class on m vertices, add vertex m with every
nonconstant cut (2^m-2 choices). Because the old tournament is strong, each
such child is strong. Keep only children whose deletion of root 0 is
strong and whose exact Delta is at most 44. Quotient by root-preserving
isomorphism, not by count. Rooted canonicalization fixes 0 and exhausts
all permutations of the remaining vertices within their full-tournament
out-degree classes. Any rooted isomorphism preserves those degrees, so this
is a complete isomorphism test. Strongness is tested by reachability.

The direct path DP is a forward recurrence on subsets and their last
vertex. Its same table counts H(T-0) on the subset with bit 0 removed.
All arithmetic is integral; the largest examined child has order ten, so
even the crude 10! bound is safely within the signed 64-bit type.

The complete retained banks are:

| order | rooted classes with Delta<=44 | minimum Delta |
|---:|---:|---:|
| 4 | 2 | 2 |
| 5 | 14 | 4 |
| 6 | 140 | 6 |
| 7 | 456 | 10 |
| 8 | 434 | 16 |
| 9 | 175 | 26 |
| 10 | 0 | none |

At order ten, 175*510=89,250 extensions were checked. Exactly 88,900 also
preserve strong root-deletion. None has Delta<=44; the minimum among these
extensions is 46, with H=135 and H(T-0)=89.

By the descent lemma, a missing bounded pair cannot first reappear at a
larger order. Therefore:

**FINITE-EXACT theorem with all-order consequence.** The minimum possible
Delta(T,v), over strong deletion pairs of orders n=4,...,10, is

    2, 4, 6, 10, 16, 26, 46.                                (4)

For every n>=10, every strong deletion pair satisfies Delta(T,v)>=46.
The first six entries happen to be 2 F_(n-2), with F1=F2=1. At n=10 that
prediction is 42, and is false. No formula for the unrestricted minimum at
arbitrary n is proved here.

An explicit order-ten witness, with v=0, has outgoing bitmask rows

    (526,4,8,530,647,543,575,879,639,6).

The first occurrence in the recursion and both counts are generated in
`rooted_cap44/step10.log`. The unrooted interval-chain construction in the
companion note supplies a second description of the same minimum: the
chain with parts (3,2,4), deleting the middle vertex of its triangle block.

## 6. Reproduction and complete records

From the directory containing this note:

    g++ -O3 -std=c++17 rooted_ears.cpp -o rooted_ears
    python3 run_rooted.py 44 10

The generated `rooted_cap44/` contains the initial rooted representatives, every retained
canonical rooted bank, every raw retained child presentation before
quotienting, every extension log and the summary. `run_rooted.py` uses no
external data. The earlier order-eight bank sweep is exploratory corroboration,
not an input to this proof. The cap128 continuation is separate and is not
silently imported into the table above.

Independent review replayed the finite universe and checked the
root-preserving canonicalization, both strongness tests, the marked-root
bitmask after deletion, and the two cited inherited inputs. Normal and
optimized Python runs agreed. Counts of presentations are not counts of
isomorphism classes, and the all-order conclusion is only the explicit
Delta>=46 consequence, not global spectrum completeness.
