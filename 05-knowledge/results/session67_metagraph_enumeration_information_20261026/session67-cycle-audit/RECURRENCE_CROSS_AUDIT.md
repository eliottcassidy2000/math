# Cross-audit: prospective symmetry and canonical root selection

**Verdict: PASS. No mathematical correction required.** The infinite same-core
obstruction in `session67-recursion/REPORT.md`, equations(4)–(6), and its
canonical-root-selector proof are sound. This note also proves a stronger
single-child separation while keeping it distinct from the aggregate gap.
No producer files or shared repository files were changed.

Audited peer report SHA-256:
`356cad6322edb0541ce29a692a5cfa9b52eb076589d43af8950cb33b722639c2`.
Audited `verify_recursion.py` SHA-256:
`87d100964c6836de631e79524c7540f7926c788ec8b3daeca6a0d4b3ae070ab6`.

## 1. Exact mass and prospective stabilizers

For parent T of order n and a target set G of child isomorphism classes,

    K_n(T,G) = sum_(z:T+x_z in G) |Aut(T+x_z)| / ((n+1)|Aut(T)|)

is the correct contribution to the ordinary child count. It is a rational
incidence mass, **not a transition probability**; values above1 are legitimate.
The root-preserving stabilizer of an attachment is exactly the stabilizer of
its mask in Aut(T). Orbit–stabilizer therefore converts the raw-mask formula
to the rooted-class sum. Summing over all root orbits of any single child gives
one, including when its automorphisms act nontrivially on the vertices.

For `T=Q[l]`, let `H=Stab_Aut(Q)(l)` and `H_i=Stab_Aut(Q)(l+e_i)`.
The old maximal-chain theorem makes the internal chain orders rigid.
There are precisely `l_i+1` distinct clone-insertion masks on chain i.
All of them yield the weighted child `Q[l+e_i]`; every other mask strictly
increases core order. Thus the peer's equation(5) is exact:

    K_n(T, core order |Q|)
      = sum_i (l_i+1)|H_i| / ((n+1)|H|).

This formula sums attachment masks, not distinct children. If two different
indices reach isomorphic weighted children, both terms are retained as the
rooted incidence formula requires. No missing quotient or double-counting
correction occurs at that step.

## 2. Infinite family: proof checked

For every `k>=2`, write

    A_k=C3[3k-3,1,1],  B_k=C3[k-1,k,k],  n=3k-1.

Aut(C3) consists of the three cyclic rotations; a positive three-vector has
nontrivial stabilizer precisely when all three entries are equal. Thus A_k
and B_k are rigid. Both are strong, have actual canonical core C3, total
order n, the same one-term ordered SCC-size sequence `(n)`, and complete
current automorphism cycle index `p_1^n`.

- No single-coordinate increment of `(3k-3,1,1)` is uniform. All A_k's
  `n+3=3k+2` clone-insertion children are rigid.
- For B_k, incrementing its unique deficient chain gives `(k,k,k)`.
  That chain has length k-1 and therefore k insertion cuts. These k masks
  give child automorphism order3; all other clone masks give order1.

Consequently

    K_n(A_k, core order3)=(3k+2)/(3k),
    K_n(B_k, core order3)=(5k+2)/(3k),
    K_n(B_k, core order3)-K_n(A_k, core order3)=2/3.

The bound k>=2 is essential for the stated positive weights. This is a
failure of universal family transfer from the specified coarse state;
it does not prove that no scalar identity for the actual global sequence
can exist. The necessity claim for blockwise constancy is correctly limited
to transfers valid for arbitrary assigned parent masses.

## 3. Stronger single-child separation — proved audit addendum

Let

    U_k=C3[k,k,k].

Then for every `k>=2`,

    K_(3k-1)(A_k,{U_k})=0,
    K_(3k-1)(B_k,{U_k})=1.

**Direct proof.** Deleting any vertex from U_k leaves one chain of length
k-1 and two of length k. All three choices of the deficient chain are
isomorphic by core rotation. Thus B_k is the only unrooted parent class
of U_k, whereas A_k is not isomorphic to B_k (their chain-size multisets
differ). This already excludes all attachments from A_k to U_k.
Over a fixed B_k, exactly the k insertion cuts of its deficient chain
produce U_k. Each has child automorphism order3 and parent automorphism
order1. The raw-mask formula gives `k*3/(3k)=1`.

Equivalently U_k has k vertex orbits, one for each rigid chain rank,
each of size3. Its k rooted classes all lie over the parent B_k.

The two separations concern different target sets: **0 versus1 is for the
single child class {U_k}; the difference2/3 is for the union of all children
with canonical core order3.** They must not be interchanged in the final
statement. Both identify prospective symmetry as information omitted by the
current core/order/SCC/cycle-index state.

## 4. Canonical root selector — proof and implementation checked

A complete rooted canonical code c(U,v) has equal values exactly when the
rooted tournaments are isomorphic. For a fixed unrooted U, this means equal
values exactly on Aut(U)-orbits. Its minimum therefore selects one entire
root orbit and no second orbit.

There is exactly one rooted extension class over every pair consisting of a
parent isomorphism class and an attachment-mask orbit under its automorphisms.
Selecting precisely the new roots in the minimum-code orbit consequently
accepts every child once: existence follows by deleting a selected root;
uniqueness follows because all selected roots belong to the same rooted class.
The induction needs a complete preceding parent bank and one representative
per mask orbit, exactly as stated. A final class set can audit the result but
is not needed to repair duplicate accepted presentations.

The implementation's `root_key` puts the distinguished root first, then sorts
other vertices by the invariant pair `(relation to root, total outdegree)`
and minimizes over all permutations inside each equal-key cell. Every rooted
isomorphism preserves those cells, and equality of resulting full adjacency
codes gives a rooted isomorphism. Therefore this restricted minimization is
still a complete rooted code. It need not equal the numerical minimum over
*all* root-first permutations, only induce exactly the same equivalence.
The code uses the new vertex's position in the canonical ordering correctly.

The two-vertex source/sink example defeats unrooted-deletion-card minimization:
both cards are the singleton but the roots are in different orbits.
The triangle defeats dividing the rooted-class count by child order:
there are4 rooted classes at order3 and2 unrooted classes, not4/3.
Neither hostile invalidates the complete rooted-code selector.

## 5. Independent finite checks and frozen files

`audit_recurrence_cross.py` imports no producer automorphism or extension
routine. It constructs the weighted triangles directly and enumerates
automorphisms using exact arc-consistent bijection backtracking. The producer's
`root_key` alone is imported to compare its equivalence against unrestricted
brute root-first canonical codes.

Normal and optimized runs agree byte-for-byte and perform **3,195 checks**:

- All clone masks for both parents at every k=2..12 (506 masks), checking
  parent rigidity, canonical chain sizes, prospective stabilizers, the2/3
  aggregate gap and0-versus1 single-child gap.
- Every attachment mask for both parents at k=2,3,4: **4,672 instances**,
  explicitly verifying that core3 masks are precisely the clone masks.
  These independently recover the peer's complete k=2 raw and weighted
  histograms; they also test parent orders8 and11.
- Independent exhaustive labelled-orbit classification through order6,
  yielding76 nonempty classes total. Complete brute rooted codes and the
  producer's pruned codes agree on every within-class root equivalence
  (2,403 ordered root pairs); selected roots are exactly one actual
  automorphism orbit in each class.
- The rooted counts are `1,2,4,12,48,296` at orders1..6, matching the
  independently known incidence values.

The all-k conclusions are proved in sections2–3, not inferred from those
finite tests. The producer's full order8 run was inspected but was **not**
replayed in this cross-audit; no such replay is needed for the all-order proof.

Reproduce from this owned audit directory:

    python3 audit_recurrence_cross.py > RECURRENCE_CROSS_RECEIPT.json
    python3 -O audit_recurrence_cross.py > RECURRENCE_CROSS_RECEIPT.optimized.json
    cmp RECURRENCE_CROSS_RECEIPT.json RECURRENCE_CROSS_RECEIPT.optimized.json

Frozen deliverables are this note, the script, the two receipt JSON files,
and `RECURRENCE_CROSS_SHA256SUMS`. The earlier core-triangle audit and original
cycle-index audit remain separate, unchanged artifacts.
