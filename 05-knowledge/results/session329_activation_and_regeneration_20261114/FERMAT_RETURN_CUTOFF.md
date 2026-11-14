# Two-port Fermat returns miss almost all prime targets

**Theorem.** Start from any finite bank of rooted tournaments. At each step use a one-way ordered seam with (R_j), (j=1,2,4,8,\dots) (the Fermat exponents), and choose either the old root or the new terminal as next root. There is an (s<1) for which this entire grammar supplies (O(X^s)) counts at most (X) that are prime or seven times a prime. In particular no cofinite Chan–Pak cutoff follows from these two ports, even allowing arbitrary switches and infinitely many Fermat block lengths.

## Transition and the obstruction to short returns

Use the inherited [ordered Busch overlap](../session324_busch_inverse_and_gap_structure_20261113/BUSCH_INVERSE.md) and [rooted Fermat clock](../session328_fermat_gluing_and_phase_repair_20261114/FERMAT_CLOCK_MEMORY.md). Put (H=H(T)), (h=H(T-v)) (the empty tournament has count one) and (L) equal to the relaxed two-piece cover count: the suffix must be empty or start at an out-neighbor of (v). Always (L\ge H\ge h>0), the second inequality by path insertion. At a seam the new total is (H'=H+2^j L\ge3H). At the retained root,
[
(h',L')=(h,2^{j+1}L-h),
]
while at the new terminal,
[
(h',L')=(H,H').
]
For the latter, deleting the terminal leaves the source followed by a transitive chain, and the new terminal's only exit is to the old root; the zero-defect criterion gives (L'=H'). The singleton boundary also follows directly from the pointed (R_j).

Call a move to the new terminal a reset. Just after a reset (L=H). One further letter that itself resets, starting at count (H>1), produces
[
(H',h',L')=((1+2^j)H,H,(1+2^j)H).
]
Subsequent arithmetic at either port preserves their common divisor (H). A prime output is therefore impossible. For an output (7p) with (p) prime, the only possible proper divisor (H>1) other than the globally forbidden count seven is (p). The short block itself has output ((1+2^j)p\ne7p), since (2^j\ne6). One more move would already give at least (9p>7p). Hence to reach a prime or (7p), every block between successive resets must use at least two letters. The initial prefix up to the first reset and the terminal suffix after the last reset are not subject to this restriction. After any first move (H\ge3), so (H=1) needs no later exception.

## Count every admissible itinerary, not just a chosen strategy

For a block of (k\ge1) moves with no interior reset, each of its first (k-1) steps multiplies (L) by at least (2^{j_i+1}-1). The last step takes the total to at least (2^{j_k}L). Since initial (L\ge H), every block increases (H) by a factor of at least
[
g(j_1,\ldots,j_k)=2^{j_k}\prod_{i<k}(2^{j_i+1}-1).
]
Define (A_s=\sum_j(2^{j+1}-1)^{-s}) and (B_s=\sum_j2^{-js}) over Fermat exponents. These are continuous near (s=1), by uniform geometric domination. The weight sum of admissible internal blocks of length at least two is bounded by (D_s=B_s A_s/(1-A_s)). At (s=1) even replacing all future Fermat exponents by every integer (j\ge8) gives
[
A_1\le 1/3+1/7+1/31+1/128<0.517,\qquad
B_1\le 1/2+1/4+1/16+1/128<0.821,
]
so (D_1<1). Fix (s<1), positive and close enough to one that (A_s,D_s<1). The unrestricted first and last segments each cost at most (1+B_s/(1-A_s)), counting the empty choice. An all-retained word may be assigned to one such segment. The allowed whole-word weight is therefore bounded by
[
\sum_{\text{words}} \prod_{\text{blocks}} g^{-s}
\le \frac{(1+B_s/(1-A_s))^2}{1-D_s}.
]
A word ending at (H\le X) has (prod g\le X/H_{seed}), so it contributes weight at least ((H_{seed}/X)^s). Sum over the finite seeds. Counting itineraries bounds distinct counts even when many words give the same integer. The elementary Chebyshev lower bound (\pi(X)\gg X/\log X) then implies almost all primes, and almost all seven-times-prime counts, escape the grammar.

This is not a global spectrum exclusion. A third root, a growing or target-chosen source bank, or a seam with a second simultaneously active module is outside the theorem. Its point is that arbitrary switching is not enough when a reset can permanently lock all interface coordinates to a common divisor.
