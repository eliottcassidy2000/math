# Which Fermat gluing reaches 60,000?

The finite positive claim must name its gluing operation. “Onto at most 9” describes the **source** order, not the resulting tournament order. There are three materially different mechanisms in the repository.

## 1. Source lineage and the older non-Fermat alternative

- [THM-4744](../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md) gives Moon's interval chains: overlap the nearly transitive blocks at endpoints with the specified orientations, and their path counts **multiply**. [THM-4749](../../../01-canon/theorems/THM-4749-the-repunit-degree-law-phi-2-of-2-is-the-triangle-phi-3-of-2-is-the-gap-and-fermat-products-reverse-the-order.md) applies such products to Fermat and binary repunit counts. This is not a coverage certificate for all odd targets in an interval.
- The [Fermat/Monster census packet](../fermat_monster_hamiltonian_counts_prime_classes_20261102/README.md) pins the complete order-nine H-value list. By adjoining a source or sink, the path count is preserved, so its 1,520 values are precisely the counts on order **at most** nine; its maximum is 3,357. The replay keeps the 1,520 counts extracted from the pinned output and the needed witness, not an order-nine census rerun.
- An earlier exact alternative needs **no Fermat gluing**. The strong-count inventory in THM-4744 supplies 69..83; Busch's strong floor puts these realizers at order at most nine. [THM-4097](../../../01-canon/theorems/THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md) supplies the order-nine range from 85 through 2,881. [THM-4102](../../../01-canon/theorems/THM-4102-selected-order-ten-strong-ear-solid-interval.md) gives every odd 249..14,649 from selected strong order-nine parents, and [THM-4104](../../../01-canon/theorems/THM-4104-selected-order-eleven-strong-ear-solid-interval.md) gives every odd 429..80,265 via another selected one-vertex ear. Together these independently audited finite packets cover every odd 69..59,999 on at most eleven vertices. This is a selected realization certificate, not a complete census of the larger orders.

Searches of the pinned repository and searchable main did not locate a prior specific 69..60,000 Fermat-block coverage packet. The closest operational provenance is the S324 one-way Busch overlap below; the new full check using it is a separate finite result. Do not attribute this full check to THM-4744 or the old ear packet.

## 2. A complete multiplicative countercheck

On the 29,966 odd targets 69..59,999, the following exact closure computations start from the complete at-most-nine value set \(A\). Repeated factor insertion is allowed, even if one wants only a single block.

| multiplication grammar | misses |
|---|---:|
| one value of \(A\), times arbitrary actual Fermat-number factors below 60,000 | 22,013 |
| one value of \(A\), times arbitrary \(2^j+1\) factors below 60,000, \(j\ge1\) | 21,209 |
| arbitrary products of values in \(A\) | 10,849 |
| arbitrary products of values in \(A\) and all \(2^j+1\) factors below 60,000 | 10,849 |

Every row first misses **2,887**, a prime, which is absent from \(A\) and not any \(2^j+1\). The equality of the last two rows is structural rather than a fitted total: the only positive near-transitive counts below 60,000 not already in \(A\) are
\[
4097=17\cdot241,\quad8193=3\cdot2731,\quad
16385=29\cdot565,\quad32769=11\cdot2979,
\]
with both factors in \(A\) in each case.

Ordinary primes and \(7p\) give a sharper interpretation. There are 5,622 prime and 983 \(7p\) targets here that are absent from \(A\); no grammar in the table supplies any of them. For a prime, pure multiplication could only use that target itself as one factor. If \(2^j+1\) is prime then \(j\) is a power of two (use the odd-power factorization when \(j\) has an odd factor); all relevant possibilities here are already in \(A\). For \(7p\), the only nontrivial factorization is \(7\cdot p\), with \(7\) absent, and \(7\) divides no \(2^j+1\) because the powers of two modulo 7 are \(1,2,4\). This explains why a finite product bank is structurally wrong for this coverage claim, even though true Fermat gluing can be very useful after changing the interface.

## 3. The nonmultiplicative gluing genuinely changes the answer

The [S324 Busch inverse](../session324_busch_inverse_and_gap_structure_20261113/BUSCH_INVERSE.md) attaches a nearly transitive block of order \(j+2\) at a shared root \(v\), with all other cross arcs one-way. For the specified rooted source \(A\), write \(H=H(A)\), and \(L\) for the partition count on \(A-v\) into two directed paths, allowing empties and retaining the outgoing root guard on the second slot. The actual law is
\[
N=H+2^j L=(2^j+1)H+2^j(L-H).
\]
Here the last term is a nonnegative **boundary defect**, not a scaled Fermat factor. If \(j\) is a power of two, the isolated block count is an actual Fermat number. The glued count need not be divisible by it. In particular no prime target can be produced by \(D=L-H=0\) from a nontrivial strong source \(H>1\), and neither can \(7p\) for prime \(p\): zero defect gives an actual product of two realized factors, one of which would have to be the forbidden 7. Matching a divisor structure is not a replacement for selecting a positive defect at a precise source.

For a literal distinction, use the order-nine strong source with upper-triangle bits
\[
111110101111111111100111111101111111
\]
(row by row, 1 meaning the smaller index beats the larger), with root 8 under zero-based labels. Then \(H=911\), \(L=988\). Attaching just the triangle block (\(j=1\)) produces
\[
2887=911+2(988)=3(911)+2(77),
\]
the first multiplicative hole. The replay counts all paths of the completed eleven-vertex tournament by direct subset DP and checks strongness, independently of its evaluation of \(L\). This is a native source-to-target operation, not an extrapolation from the old ear counts.

No one-step chart with a fixed finite bank of at-most-nine rooted sources can prove asymptotic coverage. For each fixed pair \((H,L)\), unrestricted \(j\) supplies only \(O(\log X)\) targets through \(X\); restricting \(j\) to powers of two reduces this to \(O(1+\log\log X)\). These all-order obstructions do not contradict the new finite certificate and do not constrain a changing source bank or genuinely recursive sequence of gluings.

Replay: \`python3 replay_glue_scope.py\`. The bundled compact count snapshot and its historical source are bound by SHA-256 in the verifier; the original census came from \`origin/main\` commit \`2166f15c46d340aa7d929cd15552d2c9d6713183\`, path \`05-knowledge/results/fermat_monster_hamiltonian_counts_prime_classes_20261102/data/census_9.out\`.
