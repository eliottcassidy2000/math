# The exact H27 leaf-majority threshold is order23

Status: PROVED corollary of the complete H27 SCC inventory and paired-black degree convention in PROOF.md. This is a population threshold inside ONE fixed Hamiltonian-count stratum. It is not an eventual threshold for realizing all admissible counts, and does not settle spectrum completeness or graph connectivity. No priority claim is made.

Let M_n be the number of converse-merged classes with H=27, and L_n the number of those nodes that are leaves in the FULL fundamental paired-black graph Q_n. For every n>=6,

    N_n = n(n-5)(n-4)/6,
    S_n = [n odd](n-5)/2,
    M_n = (N_n+S_n)/2.

Here N_n and S_n count ordinary and self-converse classes respectively. The previous global strong-atom inventory determines these identities, not extrapolation.

Every H27 leaf is a nonpalindromic word in singleton and C3 components with exactly three C3 letters. Indeed the complete inventory gives automorphism order27 only for these three-triangle stacks. The other H27 SCC types have automorphism orders1 or3, hence non-self-converse degrees27 or9. A self-converse mixed leaf cannot occur at H27 by the inherited Session15 leaf classification: its raw count would be5*9^a. The self-converse rigid words instead have degree0.

There are

    A_n = binom(n-6,3)

three-triangle words: their total number of letters is n-6, of which three are C3. The binomial is zero when its upper argument is less than3. A palindromic such word must have C3 at its center, one C3 on each side, and an even number n-9 of singleton letters. Thus its count is

    P_n = [n odd](n-7)/2

for n>=7, with P_6=0. Equivalently the four singleton gaps around three ordered triangles must satisfy g_0=g_3 and g_1=g_2, so g_0+g_1=(n-9)/2. This gives(n-7)/2 choices for odd n>=9; at n7 there are none. Nonpalindromic words pair freely under reversal, giving

    L_n = (A_n-P_n)/2.

The exact fractions at the relevant orders are

- n14:28 leaves among105 merged nodes, fraction4/15;
- n22:280 among561, strictly less than half;
- n23:336 among660, fraction28/55, strictly greater than half;
- n24:408 among760, fraction51/95, also greater than half.

The claim is that n23 is the FIRST majority, and every later order remains a majority.

## Proof of the all-order threshold

The inequality L_n>M_n/2 is equivalent to2A_n-2P_n>N_n+S_n. Multiplying by6 gives

    E(n)=n^3-33n^2+272n-672 >0                for even n,
    O(n)=n^3-33n^2+263n-615 >0                for odd n.

These formulae are valid for n>=6 with the zero-binomial convention above. Check the finite n6,...,22 values; all have the corresponding polynomial negative. The final ones are E(22)=-12 and O(21)=-384. At the first positive orders, O(23)=144 and E(24)=672.

Their two-step differences are

    E(n+2)-E(n)=6(n-10)^2-180,
    O(n+2)-O(n)=6(n-10)^2-198.

Both are positive for all n>=23. This proves positivity throughout both subsequent parity progressions, and hence the exact first-majority threshold23. None of these thresholds concerns existence of raw H-values.

Run `python3 verify_majority.py` and `python3 -O verify_majority.py`. Besides exact arithmetic, the checker directly enumerates every singleton/C3 word with exactly three C3 letters at orders6 through40 and forms its reversal orbits, independently checking the leaf and palindrome counts. It verifies the finite initial sign range and polynomial difference identities; the displayed positivity argument establishes the unbounded tail.
