# Independent orbit-counting audit and recurrence obstruction

Session 175, 2026-11-04. **PROVED** statements below are elementary deductions from the cited Burnside formula. **FINITE-EXACT** computations are separately recorded. No historical-priority or polynomial-runtime claim is made.

## Sources and inheritance

The primary definition [OEIS A000568](https://oeis.org/A000568) is the number of unlabeled tournaments on n vertices, with offset zero. The page explicitly supplies Davis's odd-partition formula. The repository already contains the exact future-response quotient in `05-knowledge/results/tournament_future_boundary_20261020/README.md` and the odd-factorial packed transfer in `05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/local_scale.py`. The one-tail/core-kernel split is likewise inherited. The ordinary residual DP below is an independently implemented backward form of that existing quotient, not a new theorem.

The companion Session175 implementation in `../a000568/projective.py` reproduces the fixed-target projective quotient already proved in [Session98, Sections3–4](../../session98_tournament_growth_20261029/recurrence/PROJECTIVE_RESPONSE.md); [Session105, Sections2–3](../../session105_metagraph_recursion_20261029/enumeration/INTEGER_PROJECTIVE.md) already supplies its odd-local integer realization and packing. The present folder supplies an independent backward implementation, a different normalization, a signed-edge permutation oracle, and a recurrence-class obstruction.

## 1. Fixed-target projective equivalence is exact and minimal wordwise

At a fixed processing stage let J be the allowed future odd cycle lengths, including 1. Let M be the exact remaining mass, and retain only j<=M. A processed history has response R_j=sum over past cycles a of gcd(a,j). A future cycle multiplicity vector k has sum j k_j=M and contributes a history-dependent cross factor 2^(sum k_j R_j).

Two responses R and R' induce proportional factors on **every exact-mass completion** if and only if

    R'_j - R_j = j c for every retained j,

for an integer c. Necessity: the all-singleton completion gives cross difference M c, with c=R'_1-R_1. Compare it to one j-cycle plus M-j singletons. The difference is constant only when R'_j-R_j=j c. Sufficiency follows immediately by summing against k and using its mass equation.

Thus Q_j=R_j-jR_1, j>=3, is a complete invariant for this wordwise projective equivalence. The eliminated factor is 2^(M R_1), which must be absorbed into the historical weight; discarding R_1 without this compensation is incorrect. This is minimal only for preserving each individual allowed completion up to one scalar, not for all conceivable symbolic aggregate algorithms. M=0 has no response and is a separate trivial terminal case. The singleton assumption is essential.

A concrete nontrivial merge occurs for past histories (27) and (9,9,9), both mass 27, when only future lengths 1 and 3 remain. Their responses are (1,3) and (3,9), hence both Q_3=0. At target 30 the cross factors on completions (1,1,1) and (3) are constant within each history, and differ between histories by 2^6. Their original Burnside weights still differ and are added only after compensation.

## 2. Independent integer backward recurrence

Let F(R,L,s) be R! times the remaining Burnside sum for residual mass R, allowed odd lengths <=L, and past response s_j for odd j<=L. Let f_i=i!. For c new L-cycles, set r=R-cL and

    E = c s_L + L c(c-1)/2 + c(L-1)/2,
    C = f_R / (f_r L^c c!).

Then

    F(R,L,s) = sum_c C 2^E F(r,L',s'),
    L' = largest odd integer <= min(L-2,r),
    s'_j = s_j + c gcd(L,j).

The coefficient C is integral: it counts ways to select and arrange c disjoint L-cycles on R labeled positions. Terminals are F(0,...)=1 and F(R,1,s)=2^(R(R-1)/2+R s_1). Therefore a(n)=F(n,n_odd,0)/n!. The source code checks all divisions with explicit exceptions, which remain enabled under python -O.

For independent projective normalization choose

    b = min_j floor(s'_j/j),
    u_j=s'_j-bj.

All u_j remain nonnegative, and F(r,L',s')=2^(br) F(r,L',u). This canonical representative is invariant under integer mass-direction shifts. It avoids the companion forward algorithm's signed exponents, trading additional tuple normalization for a simpler integer proof. It is not claimed faster.

## 3. No finite-order polynomial-coefficient recurrence

A000568 is **not P-recursive over the rational or complex numbers**, and its ordinary generating series is not D-finite. This does not rule out recurrences with exponential coefficients, growing order, nonlinear recurrences, or the finite-horizon state algorithms above.

Proof: there are 2^(n(n-1)/2) labeled tournaments and each isomorphism class contains at most n! labeled copies. Consequently

    a(n) >= 2^(n(n-1)/2)/n!.

Suppose sum_{i=0}^r p_i(n)a(n+i)=0 for fixed polynomials p_i, not all zero. Choose the largest r with p_r nonzero. For sufficiently large integer n, p_r(n) is nonzero, and every ratio |p_i(n)/p_r(n)| is bounded by C n^D for constants C,D>=0. The maximum absolute value of the sequence in a sliding window therefore grows by at most max(1,r C n^D) per step. Iteration yields |a(n)|<=exp(O(n log n)). This contradicts the displayed lower bound, whose logarithm is (log 2)n^2/2-O(n log n). Order zero is immediately impossible since a(n)>0 and a nonzero polynomial has finitely many roots. The OGF equivalence with P-recursiveness can also be obtained directly by comparing coefficients of a hypothetical polynomial-coefficient differential equation.

## 4. Verification and boundaries

Run `python independent_audit.py` from this directory. It compares six independent-direction/normalization representations for every n=0..30, including both companion tuple and packed modes. For n=0..7 it independently loops over **all permutations**, follows their signed action on unordered edges, rejects an orbit when its orientation flip parity is odd, and averages 2^(number of consistent edge orbits). This oracle does not invoke the odd-cycle restriction or gcd formula. It also tests the projective-equivalence criterion on 600 arbitrary response-difference vectors and all exact-mass completion vectors through M=12.

Larger comparisons at n=50,100,150,200 record hashes and timings in `independent_results.json`. Timing comparisons concern these implementations in this environment, not external best-known algorithms. The backward cache counts are not directly comparable to forward peak-live-state counts. All algorithms retain an unproved asymptotic state-complexity bound; no polynomial-time theorem follows from finite state compression. Canonical augmentation produces representatives and is output-size bounded below by a(n), so it cannot be competitive with counting-only arithmetic at large n when the desired output is a single integer.

## 5. Final packed-kernel optimization audit

The final companion kernel uses field width `max(1,n.bit_length())`, not the looser quadratic bound. This is safe: at any descending-length stage every past cycle length a exceeds every retained future length j. Hence each summand j-gcd(a,j) is nonnegative and strictly less than a, so S_j=sum_a(j-gcd(a,j))<d<=n for nonempty history (zero for empty history). New-cycle updates obey the same bound on remaining coordinates, excluding carries between packed fields. The ungauged responses also satisfy R_j<=d<=n.

The transfer multiplier is odd, so exact divisibility by a power of two can be checked and executed on the historical integer before multiplication. Checking the low-bit mask against zero is equivalent to a zero remainder modulo that power. Final normal and optimized independent receipts include the SHA-256 of the exact companion source tested.
