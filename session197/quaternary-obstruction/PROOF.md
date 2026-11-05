# A general unit-clock reduction and a quaternary obstruction

Session 197. **PROVED:** every eventual transitive-root-clone clock A b^m+epsilon, with real A>0, b>1 and epsilon in {-1,1}, has integer b>=2, positive integer A divisible by (b-1)!, and only singleton SCCs outside its nontrivial strong root core. In particular, a quaternary clock A 4^m+epsilon requires A>=60, without a symmetry assumption. This is a necessary bound, not an assertion of existence at 60 or any larger amplitude.

## 1. The factorial and SCC reduction

For a rooted finite tournament (T,r), replace r by a transitive module of m>=1 vertices, retaining every external incidence. Let H_m be its number of directed Hamilton paths. Splitting a Hamilton word into its maximal nonempty clone runs gives

    H_m=sum_(j=1)^K a_j j! S(m,j),

where a_j are nonnegative integer skeleton counts and a_K>0. The maximum number K is finite. All old fragments between consecutive clone runs are nonempty paths starting in O=N+(r) and ending in I=N-(r). Holding the two exterior fragments fixed, permuting the j-1 internal fragments always gives another valid skeleton. This is a free permutation action: disjoint nonempty vertex words distinguish the fragments. Thus

    (j-1)! divides a_j.

Since j! S(m,j) is a linear combination of 1^m,...,j^m with top coefficient 1, the greatest exponential base of H_m is K, with coefficient a_K>0. An eventual identity H_m=A b^m+epsilon with A>0 and b>1 therefore forces b=K and A=a_K. These are integers, and (b-1)! divides A. Subtracting the putative identity leaves a finite combination of distinct positive integral bases vanishing eventually; an invertible Vandermonde system makes every coefficient zero. Hence the identity holds for every positive m. No value at m=0 is asserted.

Let g be the gcd of all positive-m values. Subtracting bH_m from H_(m+1) shows g divides b-1. The factorial condition implies b-1 divides A, including the trivial b=2 case. Consequently H_m is epsilon modulo b-1, and g=1.

The SCCs of a tournament form a total order, and a Hamilton path traverses them in that order. Its count is the product of their Hamilton-path counts. Nonempty transitive cloning preserves this ordered SCC structure, except for replacing the root vertex within its component. Every component outside the root component contributes a fixed positive factor to all H_m. Since their product divides g=1, every such factor is 1. A nontrivial strong tournament contains a directed triangle, whose three Hamilton paths extend injectively under induced vertex insertion. It therefore has H>=3. All exterior SCCs are consequently singleton SCCs, forming independent transitive prefix and suffix padding.

The root component cannot itself be a singleton, since cloning a singleton SCC leaves the total count constant. Removing the exterior singleton SCCs leaves a nontrivial strong rooted core with exactly the same H_m at every positive m. This reduction covers arbitrary weak original sources and arbitrary padding lengths; no strong-order bound is applied to their total padded order.

## 2. Quaternary finite reduction

For b=4, Section 1 gives 6|A. If A<60 then A is one of

    6,12,18,24,30,36,42,48,54.

The strong core has H_1=4A+epsilon<=217. Busch's strong-tournament lower bound gives 5^(n-1)<=H^3. Since 5^11>217^3, orders n>=12 are excluded. The exact order-eleven branch gives H>=9*5^2=225>217, excluding n=11 separately. Thus every candidate core has order at most ten. The sharp order-eleven step cannot be replaced by the weaker uniform bound alone.

We use the already independently reviewed complete bank of all strong tournaments with H<=217 through order ten. It has class inventory

    n=3,...,10: 1,1,6,35,353,1144,726,190,

for 2,456 classes in total. Its construction and isomorphism completeness are inherited from the published equality packet below, rather than claimed as new in this packet. The bank includes every such strong graph, not merely graphs or roots passing that packet's positive-ternary filter.

## 3. Two independent finite exclusions

The first verifier constructs actual clone graphs and counts their Hamilton paths by ordinary endpoint subset DP. It independently recounts all 2,456 original H values. Among the 410 graph/sign choices for which H_1=4A+epsilon with a listed A, it tests every one of the 3,426 roots, without a matching-rank filter. A necessary quaternary identity is

    H_2=16A+epsilon, H_3=64A+epsilon.

Exactly two roots pass the first equation and none pass the second. Thus no candidate core exists. The optional m=4 test is reached by no root and is not needed for exclusion.

A second counter independently computes root-run coefficients by an exterior-word gap polynomial. For each old-vertex ordering, a legal ordinary edge contributes 1 at its gap; a legal insertion through r contributes z; both choices are retained when both are legal. The two exterior gaps are included. This is evaluated by subset DP on the old vertices, not by counting actual clone graphs. Matching rank three is computed by augmenting paths; all 2,410 such roots in the bank are tested. The exact unit-quaternary coefficient target is

    (a_1,a_2,a_3,a_4)=(4A+epsilon,6A,4A,A).

There are no hits. The gap implementation caps each coefficient at 1,024. Every coefficient in the target range A<=54 is below that cap; componentwise nonnegative addition and shifts make saturation exact for this comparison. Five 12-bit lanes hold the coefficients: a pre-cap sum is at most 2,048, so no carry crosses a lane. The maximum rank bounds all higher coefficients by zero. Its broader A<=125 comparison is only an extra control; the all-order theorem here is the A<60 exclusion.

Together with Section 2, the finite exclusion proves A>=60 for every eventual pure unit quaternary clone clock of either sign, including weak sources. We do not extrapolate the finite bank to larger amplitudes, other bases, or arbitrary graph-changing constructions.

## 4. Dependencies and reproduction

The complete strong bank is copied byte-for-byte as `strong-bank.jsonl`, SHA256 `3d4ddbe02807c7c6bf962147c582e0c6707eb4d221134a715cdd8a84bf730951`. The published source is:

https://github.com/eliottcassidy2000/math/blob/a7d7c95bb3532507b63025b0f74d702ced0e692c/session197/minimum-positive-clock/equality/PROOF.md

Its proof SHA256 is `02363cc76958663bae028df5bfc0d15f03c138866cd5881d3cbebb2d9b590934`. The adjacent independent `REVIEW.md` has SHA256 `6daff2578b8926effb4c6890b7a1003fa51dbaf2c89a39c2aaff80b14d36b161`. The copied bank, strong-order cutoff and complete-generation scope were verified against those named sources.

The external bound is Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC 13 (2006), N3, Theorem 1, DOI 10.37236/1141. Its exact branches are already checked in the bank's proof and independent review.

The run-skeleton and factorial mechanisms are inherited from the session's reviewed matching and negative-base rigidity work; Section 1 states and proves the general unit-offset consequence needed here. There is no literature-priority claim. This packet's new bounded conclusion is the unrestricted quaternary amplitude obstruction using the complete H<=217 bank.

Run `python3 verify.py` and `python3 -O verify.py`. They compile both included standard-library C++ counters in a temporary directory, verify the bank hash/inventory and exact numerical results, and check the strict arithmetic cutoff. No external data fetch or probabilistic test is used.
