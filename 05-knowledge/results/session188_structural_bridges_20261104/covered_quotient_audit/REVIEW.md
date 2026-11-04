# Independent audit of unequal transitive cyclic blocks

**PASS. No mathematical correction found.** The producer files were not changed. This audit distinguishes the unbounded proof from the finite execution checks.

## Exact tomography

The three displayed actual-query formulas have the stated endpoint convention. With a=e_u and c=0, any Hamilton path beginning outside the old source chain must enter x from u; y either starts the path or immediately follows x. Both positions supply the same terminal and two-run possibilities, explaining the factor two. With a=0, x is a universal source and the remaining count is a one-port insertion; changing just u->y loses precisely paths starting at u and gains precisely those ending at u or cut after u. Their linear combination extracts 2S_u exactly. It is not an inference from a deleted card, and it does not drop higher-degree response terms.

## Cyclic-run identity

For quotient 0->1->2->0, each maximal-run word beginning in block i is forced. Every length falls into exactly one of 3k, 3k-1 and 3k-2. The associated run counts are exactly those stated. An ordered path cover of each physical child assigns its paths to successive occurrences of that child's label, giving a bijection. Reversing the quotient exchanges the next and previous block, so the equal and doubly shortened cases cancel. The identity therefore holds for arbitrary nonempty child tournaments, not just transitive ones. Its zero extension at k=0 and beyond the child orders handles all boundary terms.

## Stirling minors and strictness

The polynomial expansion for b_k=k a_k+a_(k-1) is algebraically correct. For a log-concave row with consecutive support, a_k a_(k-1)>=a_(k-2)a_(k+1) follows by multiplying successive-ratio inequalities, and the zero boundary cases preserve it. This establishes row log-concavity inductively without appealing to an unverified citation.

The adjacent-order ratio k+s(m,k-1)/s(m,k) is strictly increasing on its positive support: its second term is nondecreasing by row log-concavity, and its first term strictly increases. The product of these positive ratios gives the order-b versus order-c ratio comparison for 2<=k<=c. The k=1 minor is zero; at k=c+1 only the nonnegative first product remains; above that both products vanish. The factorial factors cancel exactly. At k=2 the minor equals 2^b-2^c, so it is strictly signed whenever b differs from c.

For unequal a,b,c>=2, one block has two differently sized neighboring blocks, and its own two-path cover count is positive. Thus the k=2 contribution is strictly signed while every other contribution has the same weak sign. This proves a nonzero block-summed start difference. Quotient converse preserves H by whole converse followed by child converses. The resulting signed battery of 2|A_i| actual queries is therefore nonzero. The proof claims exact inequality, not a residue-two separator; this distinction is necessary.

## Scope and hostile control

The theorem concerns a C3 quotient with three transitive children, sizes at least two, and quotient reversal retaining child interiors and all physical names. It is a fully covered laminar example, not arbitrary module-family closure. Equal sizes require the separate argument explicitly cited by the producer. The four-vertex deletion hostile is valid and its minor is -8, so arbitrary child ratios cannot be substituted for Stirling ratios. No invalid deleted-card tomography is used.

## Execution

`audit.py` imports no producer helpers. It counts Hamilton paths by literal permutations and child covers by explicit assignments of vertices to ordered slots. Its receipt pins both reviewed producer sources, rejects source changes during a run, and records unchanged producer normal/optimized replays. All checks use explicit exceptions and remain active under optimized Python. Exact finite counts and hashes are in `receipt.json`; these support implementation and do not replace the proof above.
