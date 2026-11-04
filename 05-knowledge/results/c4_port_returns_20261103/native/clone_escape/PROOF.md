# The homogeneous-clone escape also excludes24/44

**PROVED using finite exact certificates and cited classical order bounds:**
no fixed marked strong tournament (R,x) has

    H(R[x := TT_t]) = 63·2^t−35  or  63·2^t−29

for every sufficiently large integer t. In fact the necessary two-run
coefficient63 never occurs at any marked strong base with H=91 or97.
This result allows arbitrary surrounding tournament shape. It does not
allow that shape or its incidences to change with t, or the cloned prefix
to lose its homogeneous-module property.

## Exact inherited mechanism and the target

The closest mechanism is the clone-run formula in
[`session182_physical_returns_20261104/defect/PROOF.md`, §3](../../../session182_physical_returns_20261104/defect/PROOF.md),
independently reviewed in that packet's `defect_audit/AUDIT.md`. For a fixed
marked tournament R, delete x to obtain U. Let I be x's old in-neighbors
and O its old out-neighbors. In the squarefree ring on U, let

    L = 1 + (all nonempty old paths ending in I),
    Q = 1 + (all nonempty old paths starting in O),
    M = (all old paths starting in O and ending in I),
    A_r = [U] L M^(r−1) Q.

Then, for t≥1,

    H(R[x := TT_t]) = Σ_(r≥1) A_r r! S(t,r).

Here S(t,r) is a Stirling number of the second kind. The coefficients are
nonnegative integers, A1=H(R), and (r−1)! divides A_r. Each intervening old
path starts in O and ends in the disjoint set I, so it consumes at least
two vertices. Thus A_r=0 for r>floor((|R|−1)/2)+1.

The corrected near miss is that preserving the old transitive prefix as a
module is much weaker than preserving the special M(TT_t;B,W) attachment
chart. A general marked quotient must be tested separately. The canonical
hostile is R=[20,1,2,7,14], x=0, with A=(9,10,2): its third-run contribution
gives H=2·3^t+4·2^t−5. The least-used sidecar is the full positive run vector,
which distinguishes a true binary clock from one contaminated by higher
exponential modes. The live concepts are marked modules / clone runs /
Stirling coefficients / strong deletion / actual-source acceptance.

## Why the search is finite, even for eventual identities

The expansion

    r! S(t,r) = Σ_(k=1)^r (−1)^(r−k) binom(r,k) k^t

shows that if r is the largest index with A_r≠0, the coefficient of r^t is
exactly the positive integer A_r. Distinct positive exponential bases are
independent on every infinite integer tail. Therefore either proposed
eventual identity forces A_r=0 for r≥3. The remaining expression is

    H(R_t) = A1 + A2(2^t−2).

Comparing coefficients gives A2=63 and A1=91 for intercept−35, or A1=97
for intercept−29. This reasoning does not assume the proposed identity
was initially stated at t=1.

**CITED order bound.** Busch's [primary paper, Theorem1 and Corollary1,
pp3–4](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/)
gives minimum125 at strong order ten and a nondecreasing minimum thereafter.
Hence a strong base of H91 or97 has order at most nine. A strong clone
family cannot hide a weak base: a nontrivial directed condensation cut of
R lifts unchanged to every homogeneous inflation R_t. Conversely a strong
nontrivial base remains strong after transitive cloning.

For completeness at order nine, every strong order-nine tournament has a
strong order-eight vertex deletion. This is Moon's strong-deletion theorem,
used with its corrected citation in
[`THM-1370-h-spectrum-omits-7-21-all-n.md`](../../../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md)
and the complete ear reductions of
[`THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`](../../../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md),
Proofs3–4. Insertion monotonicity gives H(R−v)≤H(R): insert v into each
Hamilton path, retaining an injective recovery by deletion. Thus only strong
order-eight parents with H≤97 need extending.

## Complete finite universe and receipts

The inherited `regeneration/catalog8_f126.json` contains all ordinary
tournament classes through order eight, with counts

    1, 1, 2, 4, 12, 56, 456, 6880.

All marks are scanned, so no root-orbit representative assumption is needed.
The relevant complete census is

| Order | Strong H91 classes | Strong H97 classes | Marks |
|---|---:|---:|---:|
| ≤6 | 0 | 0 | 0 |
| 7 | 15 | 6 | 147 |
| 8 | 12 | 10 | 176 |
| 9 | 2 | 7 | 81 |

The order-nine extension uses exactly135 eligible strong order-eight
parents and all256 cuts per parent:34,560 births, of which34,290 are strong.
Exactly39 births have H91 or97. Canonicalization after the count filter
leaves the nine classes above. It uses an explicit **64-bit** tournament
code, since nine vertices require36 pair bits; extending the old 32-bit
code would be invalid. The canonicalizer minimizes within invariant
outdegree groups and therefore identifies precisely the ordinary classes.

All404 marked strong cases fail A2=63. The lower-order scan also checks
42 weak order-eight classes (336 marks), without using those extra checks
to claim a complete weak-base search. Thus the total finite audit has740
marked cases, and no candidate even reaches the higher-run-vanishing test.

`check_clone.py` and `check_order9.py` compute the entire A vector by
squarefree convolution on the deleted old graph. A separate direct path
DP on the graph with two literal clones checks H(R_2)=A1+2A2 for every mark.
The base H is independently recomputed in Python, rather than trusting the
catalog or C++ filter. The C++ extension stage and Python clone stage use
separate implementations. Every stored source matrix, root, and A vector
is available in the receipts.

Positive control: the actual order-eight base

    rows=[190,60,184,144,160,8,63,98], x=6

has A=(89,63,0,0), giving H=63·2^t−37. It is the accepted22 component word
TT3 with ports{0,2}, followed by C3 with ports{3,4}, attached to one prefix
vertex and a root. This verifies that the desired slope63 is attainable
and is not accidentally rejected by the evaluator. The higher-run hostile
above independently checks the nonzero A3 branch. Both controls are checked
by literal substitution for t=1,…,5.

Reproduction from repository root:

    python 05-knowledge/results/c4_port_returns_20261103/native/clone_escape/check_clone.py
    python -O 05-knowledge/results/c4_port_returns_20261103/native/clone_escape/check_clone.py --output 05-knowledge/results/c4_port_returns_20261103/native/clone_escape/receipt_optimized.json
    python 05-knowledge/results/c4_port_returns_20261103/native/clone_escape/check_order9.py
    python -O 05-knowledge/results/c4_port_returns_20261103/native/clone_escape/check_order9.py --output 05-knowledge/results/c4_port_returns_20261103/native/clone_escape/order9_optimized.json

Both normal/optimized receipt pairs must be byte-identical. The order-nine
program requires g++ with C++20 support and compiles into a temporary
directory. No enormous order-nine catalog is generated.

## What must change next

This closes the fixed finite marked-quotient escape, even after abandoning
the old root-attachment chart and allowing arbitrary external arcs. A next
candidate must change a real hypothesis: split the prefix into differently
incident blocks, allow the surrounding marked quotient to grow or depend
on t, or use a clock transformation not represented by this fixed TT_t
substitution. The result does not say the integers themselves are absent
from the Hamiltonian spectrum, and does not exclude a t-dependent grammar
whose each individual instance is an actual tournament.
