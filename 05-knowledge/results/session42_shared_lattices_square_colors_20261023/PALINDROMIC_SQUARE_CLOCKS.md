# Palindromic square clocks and the sharp order48 threshold

This extends [ordered-sum reconstruction](SQUARE_WORD_RECONSTRUCTION.md). All path clocks below retain the order of adjacent sums. Palindromic means the sum word equals its reversal, not that the vertex sequence repeats itself.

## Quotient theorem

For N>1, a palindromic adjacent-sum word on a Hamiltonian permutation requires even N=2m. Indeed reversing the path gives another realization of the same word; odd-order uniqueness rules this out. At even order the reverse is exactly its unique twin, so some d|m satisfies

x_(N-1-i)=tau_d(x_i),

where tau_d exchanges the two d-long halves of each 2d-block. Thus the first m vertices contain precisely one member of each tau_d pair.

Represent each pair by its lower member a in A_d. Define a simple quotient graph J_(N,d) on A_d, with distinct a,b adjacent when a+b+d is a square. Mark those vertices a for which 2a+d is a square. Then:

**A palindromic square-sum Hamiltonian path of displacement d exists if and only if J_(N,d) has a Hamiltonian path with at least one marked endpoint.**

Proof. The first half of the original path alternates A_d and B_d, so successive representatives a,b have actual sum a+b+d. Its final representative must be marked because the central edge joins a to a+d. Conversely orient a quotient path toward a marked endpoint; choose an initial side and alternate sides to lift its first half. Append the tau_d image in reverse order. Every vertex is used once, all noncentral edges have the prescribed square sums, and the marked endpoint supplies the central square. Reversal gives the twin. All choices of port position and side have been specified.

## Exact threshold: 48

No odd order N>1 is possible. The independent audit checks all 76 pairs (N,d) with even 2<=N<48 and 2d|N. Every candidate except (36,6) and (40,4) is excluded already by isolation, more than two leaves, or disconnection in the original cross-block graph.

For (40,4), the quotient has no marked vertex, so there is no possible central square edge.

For (36,6), delete quotient vertices {1,14}. Four components remain:

- {2,4,6,13,15,17,26,28,30}
- {3,16,27}
- {5,18,25}
- {29}.

A Hamiltonian path, after deletion of two vertices, has at most three path components. Since it would span every graph component, four components are impossible. Thus this quotient has no Hamiltonian path at all, regardless of its endpoint marks. This structural cut certificate replaces reliance on a search failure.

At (48,8), the following actual path has a palindromic square clock:

3,46,18,31,5,44,37,27,22,42,7,9,40,41,8,28,21,43,6,10,39,25,24,12,4,32,17,47,2,14,35,29,20,16,33,48,1,15,34,30,19,45,36,13,23,26,38,11.

Its central pair is 12,4 with sum16; their lower representative4 is marked. The second half is exactly tau_8 applied to the first half in reverse order. Every edge, vertex and reflection identity is independently checked.

Consequently the first square-clock twin occurs at36, but the first palindromic square clock occurs at48. These are distinct reconstruction phenomena. The numbers do not establish a Bott-periodic law or a uniform family at all subsequent orders.

## A useful hierarchy of retained information

An ordered path determines an ordered square clock; the clock determines at most two paths. A palindromic clock is exactly the case where its two paths are each other's reversals. Forgetting the clock order loses substantially more information and is outside these bounds.

The quotient proof gives an exact finite reduction: retain the divisor-block coloring, pass to pair representatives, and keep the marked central endpoint. A bare quotient graph without that mark would wrongly permit the (40,4) case. This is a direct example of the shared-port principle in a different graph family.

For unrestricted permutations, the independent [enumeration proof](audit-square-words/REVIEW.md) gives tau(m)m! palindromic words and ((2m)!−tau(m)(m!)²+tau(m)m!)/2 word classes modulo reversal. These factorial formulas do not count square-restricted paths.
