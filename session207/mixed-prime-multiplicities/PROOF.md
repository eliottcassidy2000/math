# Session 207: two coefficients and a metagraph degree law for an infinite family

Status: PROVED relative to inherited fixed-count SCC enumeration, the complete small strong inventories, and the explicitly named paired-black metagraph degree theorem. The verifier checks the polynomial algebra and finite controls; it is not a proof by numerical extrapolation. No global Hamiltonian-count completeness or exact infinite minimum-order formula is claimed.

## Inherited inputs and exact object

Let N_n(H) count ordinary isomorphism classes of n-vertex tournaments with exactly H directed Hamiltonian vertex-word paths. These are not labelled counts and not converse-merged counts.

The mixed-prime leading coefficient 1/(b!c!) and maximal-SCC argument are already explicit in spectrum_basins_20261020e/spectrum.md section4; they are not claimed new here. The additional deductions are the next coefficient, the quantified first-order concentration defect, and the mixed-prime paired-black degree law. The rational fixed-H series and eventual polynomial law are inherited from `05-knowledge/results/spectrum_basins_20261020e/spectrum.md`, section4, and `05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md`. An ordered word of k nontrivial strong components of total order t contributes binom(n-t+k,k), once n is sufficiently large. The k+1 gaps contain ordered singleton components. Every fixed H has only finitely many strong atoms, by Busch's bound.

The required complete strong inventory is

    P_3=z^3, P_5=z^4, P_9=2z^5,
    P_15=2z^5+2z^6, P_25=2z^6+z^7.

This is inherited from `05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/low_H_formulas.json` (current GitHub blob a10db7b891aad9ad50334ee6e428f194787ed73f; all bundled controls matched the remote source) and the audited strong inventories in THM-4744 / `05-knowledge/results/h_spectrum_multiplicities_20261021.md`. Busch excludes larger orders: f(6)=15>9, f(7)=25>15, f(8)=45>25. The unique H3 atom is C3; the unique H5 atom is the strong four-vertex tournament S4.

## 1. Explicit first two coefficients for every 3^b5^c

Fix b,c>=0 with k=b+c>=1, and put H=3^b5^c. Then, as n tends to infinity,

    N_n(H) = 1/(b!c!) * [n^k - (b^2+4bc+c^2+3b+5c)/2 * n^(k-1)]
             + O(n^(k-2)).                                           (1)

For k=1 the displayed linear polynomial is exact and the remainder is zero. Parameters b,c are fixed, not growing with n. Equation(1) identifies the two highest coefficients of the eventual exact polynomial, not merely a numerical approximation.

Proof. A nontrivial SCC contributes at least one prime factor to H. The maximum k components is attained exactly by b copies of C3 and c copies of S4. Their total order is t=3b+4c, and there are k!/(b!c!) distinct ordered words. Their exact contribution is

    T_n = k!/(b!c!) * binom(n-t+k,k).                                (2)

Its coefficient of n^(k-1) is [k(k+1)/2-kt]/(b!c!).

A word with k-1 components must have precisely one component containing two prime factors, and all others prime atoms. Thus its special atom is H9, H15, or H25. Using the complete totals P_9(1)=2, P_15(1)=4 and P_25(1)=3, their respective leading contributions at degree k-1 are

    2/((b-2)!c!), 4/((b-1)!(c-1)!), 3/(b!(c-2)!),

with terms omitted when a factorial argument is negative. All other words have at most k-2 components and contribute only O(n^(k-2)). Adding these coefficients and simplifying gives(1). Positivity prevents cancellation of the leading coefficient. This argument exhausts every way to lose exactly one SCC and leaves no assumption about higher strong atoms. QED.

Examples reproduce N_n(9)=(n-4)(n-1)/2, N_n(15)=(n-4)(n-3), and N_n(25)=(n^2-7n+10)/2 in their eventual ranges. For H27, (1) gives n^3/6-3n^2/2+O(n), agreeing with the inherited n(n-5)(n-4)/6. These are consistency controls, not the proof.

## 2. Quantified structural concentration

Let a uniformly selected ordinary class from this fixed-H stratum be called a prime stack if its nontrivial SCCs are exactly b triangles and c copies of S4. By(2) its count is T_n. Define

    D(b,c)=2b(b-1)+4bc+3c(c-1).

The preceding exhaustive k-1-component classification gives

    Pr(not a prime stack) = D(b,c)/n + O(n^(-2)).                    (3)

For k=1 this probability is identically zero. Equation(3) is stronger than the inherited qualitative top-pole concentration: it supplies the exact first-order defect for every member of an infinite two-prime family. It measures failure to be a prime stack, not failure to have a particular metagraph degree; some exceptional cores can have the same degree.

## 3. Direct bridge to the paired-black metagraph

Use exactly the inherited fundamental paired-black converse-merged metagraph Q of THM-4731 (`01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`), called B there. Its non-self-converse node degree is H/|Aut(T)|; the self-converse correction is a separate term, and loops use the inherited multigraph convention. This is not a statement about an arc-flip graph, the unpaired doubled network, or every object called a tournament metagraph.

Every prime stack has automorphism group of order3^b: the unique SCC order prevents component permutations; C3 contributes3, S4 contributes1, and singleton components contribute1. S4's automorphism order is1 because it divides both its five Hamiltonian paths and4!, while C3's cyclic order3 is immediate. Therefore its presentation count is

    H/|Aut(T)| = 5^c.                                                 (4)

For fixed b,c, self-converse classes have strictly smaller polynomial order than N_n(H). For a fixed symmetric word of r<=k nontrivial SCCs, its r+1 singleton gaps are paired by reversal; at fixed total their solution count is O(n^floor(r/2)), hence O(n^floor(k/2)). There are only finitely many possible cores. Thus self-converse classes are o(N_n(H)); converse-merging divides the dominant ordinary counts by two.

Consequently, as n grows, a fraction1-O(1/n) of the Q-nodes in the fixed raw-count stratum H=3^b5^c have degree exactly5^c. For c=0 this recovers the inherited almost-all-leaves theorem for powers of3. For c>0 it gives the mixed-prime degree law5,25,125,... independently of b. No claim is made that all nodes have that degree, or that(3) is the exact degree-exception probability.

## 4. What remains open

These results classify multiplicity growth and typical paired-black degree for an infinite family. They do not give a closed formula for the first occurrence m(3^b5^c), which can be governed by large efficient strong atoms rather than the asymptotically dominant padded stacks. In particular a prime stack is typical at large fixed-H order but need not minimize order. The h9 counterexample, m(9)=5 versus six vertices for two triangles, is the essential distinction.
