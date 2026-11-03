# An overlapping-cycle tournament, transfer matrix, and a fixed FRACTRAN program

Status: elementary proofs, independently controlled by direct Hamilton-path dynamic programming and small cycle enumeration. The transfer is a three-state analogue of a continuant, not a claimed rational-tangle invariant. This is a family of strong tournaments, not an exact-spectrum classification. No literature novelty is claimed.

## Construction and result

For integer m >= 0, take vertices 0,...,2m+1 and orient each edge from smaller to larger except for the m reversed edges (2i+3) -> (2i), for 0 <= i < m. Let H_m count its directed Hamilton paths. For m >= 1 the tournament is strong: the directed 4-cycles

    2i -> 2i+1 -> 2i+2 -> 2i+3 -> 2i

overlap in two vertices and jointly span the tournament. For m=0 the tournament consists of a single directed edge.

The exact ordinary generating function is

    Σ_{m>=0} H_m z^m = (1+z)^2 / (1 - 3z - z^2 - z^3).

Equivalently, H_0=1, H_1=5, H_2=17, and H_m=3H_{m-1}+H_{m-2}+H_{m-3} for m>=3. The first terms are

    1, 5, 17, 57, 193, 653, 2209, 7473, 25281, 85525, 289329.

Thus 5,17,57,... have explicit strong witnesses. In particular, the apparently plausible order-two recurrence H_m=H_{m-1}+8H_{m-2} fails at m=5: it predicts 649 instead of 653.

The transfer matrix is the unimodular nonnegative matrix

    (H_m,H_{m-1},H_{m-2})^t = [[3,1,1],[1,0,0],[0,1,0]] (H_{m-1},H_{m-2},H_{m-3})^t.

Its determinant is 1. The numerator and denominator above are coprime, so this sequence satisfies no homogeneous constant-coefficient recurrence of order <=2 over the rationals. In that precise sense the simplest constant two-state Conway-continuant analogy does not suffice; the three-state transfer is exact.

## Odd-cycle expansion and proof

We use the identity, valid for every finite tournament T,

    H(T) = Σ_F 2^|F|,

where F ranges over all collections of pairwise vertex-disjoint directed odd cycles of length at least 3; the empty collection contributes 1. A self-contained derivation follows. Let A be the tournament adjacency matrix, X=diag(x_i), and work in Q[x_i]/(x_i^2). The full path generating series is

    P(x) = 1 + 1^t X (I-AX)^(-1) 1
         = det(I - AX + 11^t X) / det(I-AX).

The coefficient of a squarefree monomial is the number of directed paths on its variables, since repeated vertices vanish. As A+A^t=11^t-I, expansion by principal minors in the squarefree ring gives

    P(x) = (Π_i (1+x_i)) det(I+AX)/det(I-AX).

Indeed, the numerator coefficient of x_S is det(I+A_S^t)=Σ_{U⊆S}det(A_U), the coefficient of x_S in the displayed product. Taking logarithms of the determinant ratio gives 2 Σ_{k odd>=3} tr((AX)^k)/k. In the squarefree ring, tr((AX)^k)/k counts directed k-cycles once. Exponentiating selects vertex-disjoint odd cycles, each with weight 2; the product Π(1+x_i) supplies all remaining singleton vertices. The full coefficient yields the identity.

Now label reversed edges by i, with interval [2i,2i+3]. Every directed cycle must use a reversed edge. Its used reversed-edge indices form a consecutive interval [a,b]: if two used indices have a missing index between them, cut between the two neighboring reversed-edge intervals on opposite sides of that missing index. The missing edge is the only backward arc crossing that cut, so a directed cycle using arcs on both sides must use it, a contradiction. The minimum and maximum vertices of the cycle lie within [2a,2b+3], by the same cut argument. The k=b-a+1 used backward arcs have disjoint endpoints; these 2k vertices fill the interval except for the two boundary vertices L=2a+1 and R=2b+2. An odd cycle must use exactly one of L,R.

For each choice there is exactly one cycle. Contract the kth backward edge as a directed block (2k+3,2k), for k=a,...,b. Between two distinct blocks i,j an arc is possible exactly when j>=i-1. If the cap is L, only block a can precede it, so the block immediately before a must be a+1, and inductively the unique block order is b,b-1,...,a. If the cap is R, only block b can follow it and the same descending order is forced. Both orders satisfy all arcs. Therefore every nonempty consecutive block of backward-edge indices carries exactly two odd cycles, distinguished by their boundary cap.

Collections use disjoint backward-edge blocks. Two consecutive cycle blocks, with no unused index between them, are vertex-disjoint only if the earlier block takes its left cap and the later block takes its right cap: the other cap choices collide with the neighboring block's endpoint. Thus three cycle blocks cannot be mutually adjacent. If at least one index is unused between two blocks, all cap choices are compatible.

For a maximal occupied run of k>=1 backward-edge indices, there are either one cycle block, in two cap choices and weight 2, or two adjacent blocks, with k-1 split choices, one compatible cap choice, and weight 2^2. Its total weight is 4+4(k-1)=4k. Let C(z)=Σ_{k>=1}4kz^k=4z/(1-z)^2. Runs are separated by one or more unused indices, so the sum over binary runs is

    (1+C(z))/(1-z(1+C(z))) = (1+z)^2/(1-3z-z^2-z^3).

This proves the generating function and recurrence.

## A family-specific 7 obstruction

Modulo 7 the recurrence state is invertible because the transfer matrix has determinant 1. Starting from (H_0,H_1,H_2)=(1,5,3), its exact period is 24, and one period of H_m modulo 7 is

    1,5,3,1,4,2,4,4,4,6,5,4,2,1,2,2,2,3,6,2,1,4,1,1.

No entry is zero, so 7 divides none of these H_m. This is an obstruction for this one explicit family, not a nonexistence claim for other tournaments.

## A fixed 13-fraction FRACTRAN program

Run the following ordered list under the ordinary first-applicable FRACTRAN rule:

    17/143, 133/85, 17/19, 23/17, 1015/69, 23/29, 31/23,
    38073/62, 31/37, 41/31, 86/287, 41/43, 13/41.

Start from N_m=13·11^m·2^17·3^5·5. The machine halts at

    13·2^H_{m+2}·3^H_{m+1}·5^H_m.

Here m may be any nonnegative integer. This single fixed program computes this specific recurrence; it is not claimed to be universal or to reproduce Conway's prime-generating game. Its fixed list has 13 fractions, related in spirit but not behavior to the well-known 14-fraction program.

Proof and phases: let the exponents of 2,3,5 be A,B,C, of 7 be scratch X, and of 11 be the iteration clock. At the entry control prime 13, the first fraction consumes one unit of the clock and changes to control 17. There are four guarded transfer loops. The first repeatedly consumes C and adds one to X, with 17 and 19 alternating, then falls through to 23. The second consumes B, adding one to each of X and C, with 23 and 29 alternating, then falls through to 31. The third consumes A, adding three to X and one to B, with 31 and 37 alternating, then falls through to 41. The fourth consumes X and adds one to A, with 41 and 43 alternating, then returns to 13. At each phase its transfer fraction precedes the fallback, so the fallback fires exactly once the source exponent is zero. Thus the net transformation is

    (A,B,C, X=0) -> (3A+B+C, A, B, X=0).

All registers are nonnegative, each loop exhausts a finite source, and each pass consumes one clock unit. Once the clock is zero, the first fraction fails and no other fraction is applicable at control 13. The initial triple is (H_2,H_1,H_0), so after m passes the stated output follows by induction. The output integer grows rapidly; no resource-efficiency or universality claim is made.

`python3 fractran.py` executes the actual integers and checks m=0,...,4. It also works under `python3 -O fractran.py`; fractran.py verifies the expected output without relying on assertions. The step counts are 0,165,714,2559,8788. In particular m=1 yields 13·2^57·3^17·5^5.

## Reproduction and scope

Run `python3 check.py` and `python3 -O check.py`. Independent subset/endpoint DP checks H_m for 0<=m<=6; explicit simple-cycle enumeration and independent packing enumeration check the cycle description for 0<=m<=4; independent reachability checks strength for 1<=m<=6. The script checks the entire 24-state mod-7 orbit. These finite calculations are controls, not the general proof. It is deliberately not asserted that the cycle expansion, recurrence, or resulting counts are new in the literature.
