# Gold certificates under attachments through an entire chain prefix

Date: 2026-10-19. Status: **PROVED**, with independent exact controls. This extends the frozen two-chain/fork packet without changing it. The general Gold Partition and 1/3–2/3 conjectures remain **OPEN**.

## Main result

A fixed two-comparison Gold certificate survives arbitrary side posets attached through a prefix of one shuffle chain, provided an explicit quadratic inequality holds. The proof uses a coordinatewise interleaving injection and a monotone coupling of uniform subsets. The permitted prefix can occupy an asymptotic fraction phi^(-2)=(3-sqrt(5))/2 of the chain.

The parameter range and universal normalized margin are sharp for this fixed strategy over the stated structural class. Failure outside the range means that some member defeats this strategy, not that the member defeats the Gold conjecture.

## 1. Construction and exact statement

Fix integers m,n>=2 and 1<=j<m. Take chains

    a_1<...<a_m, c_1<...<c_n.

An arbitrary finite cap H may sit above both chains. Take any nonempty finite side poset D and a maximal element p of D. Choose nested upsets

    W_1 superset W_2 superset ... superset W_j

of D, allowing empty sets. The only cross-relations between D and the backbone are

    a_i<w whenever w in W_i, 1<=i<=j.                       (1)

Each W_i being an upset, and their nesting in this direction, ensure transitivity. The boundary ports a_1,p are respectively minimal and maximal, and distinct.

Use the same fixed comparisons as in the earlier packet:

    first: c_n versus a_m;
    second: c_(n-1) versus a_m.                             (2)

The increasing answers select the first row and column. Assign core extensions any nonnegative flag weight

    alpha+beta I+gamma J+delta IJ,
    I=1{a_1 first}, J=1{p last},

with coefficients nonnegative and not all zero.

Set

    r=m-j, T=m+n, E=(r+n)(r+n-1),
    Q_j(m,n)=n^2-(r+1)n-r(r-1),

    Delta_j(m,n)=min{(n-m)/T, r(r-1)/E, -Q_j(m,n)/E}.       (3)

**Prefix attachment theorem.** Every member of this structural class and every allowed flag mixture satisfies

    (t_0-t_1-t_2)/t_0 >= Delta_j(m,n)                       (4)

on each realized answer path of strategy (2). Its first increasing-answer probability lies in

    [r/(r+n), m/(m+n)].                                    (5)

The constants in (3) and endpoints in (5) are sharp uniformly over the class. Here sharp means an infimum or supremum; a lower extremum may require a sequence of growing side posets rather than one finite member.

Consequently the strategy is Gold for every such member and every flag mixture if and only if

    m<=n and Q_j(m,n)<=0.                                  (6)

Under (6), the first comparison is at least r/(r+n)>=1/3 balanced. The same conclusions hold after arbitrary exterior terminal-upset and initial-ideal attachments at the two ports: their exact restriction weights belong to the flag cone, by `gold_unbounded_flag_context_20261019.md`.

Positive uniform margins occur, for example, at

    (m,n,j)=(7,8,2): Delta_j=1/39;
    (m,n,j)=(10,11,3): Delta_j=1/34.

These bounds permit arbitrary internal structure and size in D and H, arbitrary nested upsets, and arbitrary allowed exterior contexts. No strict boundary claim is made in general: empty W_i give constant fibres, which can realize zero slack on m=n.

## 2. A coordinatewise decreasing restriction fibre

Fix an extension of H and an extension tau of D. A backbone extension is a shuffle of the a-chain and c-chain, followed by the chosen cap order. Write

    t_i=position of a_i in the shuffle, 1<=i<=j.

All cross-constraints (1) depend only on the prefix arrival vector (t_1,...,t_j), not on the rest of the shuffle.

Represent an interleaving with tau by a binary word whose symbols select the next backbone or D element. If the prefix arrival vector is moved coordinatewise earlier, keep tau and this binary word unchanged. Every constrained a_i then appears no later, so all inequalities a_i<w remain valid. Thus the allowed interleaving words at a later vector inject into those at an earlier vector.

This injection preserves the p-last flag: that flag requires tau to end in p and the binary word to end with a D-symbol. It also preserves the a_1-first flag when t_1=1 is fixed: the binary word must begin with a backbone symbol, a restriction independent of the other coordinates.

After summing over tau, every relevant restriction fibre is therefore coordinatewise nonincreasing in the arrival vector. It is positive on every permitted vector: placing all of D after the backbone provides extensions, with an extension of D ending in p when that flag is required. The cap contributes a common factor e(H) and its size to the binary-word length, without affecting the monotonicity.

## 3. Averaging over earlier prefix positions

The coordinatewise statement must be converted into one-coordinate monotonicity. Without this step, decreasing dependence on several arrival positions would not by itself justify the later covariance argument.

Under the uniform two-chain shuffle, conditional on t_j=t, the earlier positions t_1,...,t_(j-1) form a uniform (j-1)-subset of {1,...,t-1}. Given that subset, the remaining tail shuffle is independent and uniform.

**Uniform-subset coupling.** Uniform k-subsets of {1,...,q} and {1,...,q+1} can be coupled so that their ordered coordinates never decrease. From a uniform subset S of size k, retain it with probability (q+1-k)/(q+1). Otherwise replace a uniformly chosen member by q+1. The resulting subset is uniform, and replacing an element by the new largest element can only increase its ordered coordinates.

Apply this coupling as t grows, and increase the final coordinate t_j itself. The coordinatewise decreasing fibre has a nonincreasing conditional average w(t).

For the a_1-first flag, t_1=1. If j>=2, the remaining earlier positions are a uniform (j-2)-subset of {2,...,t-1}. The shifted version of the same coupling works while preserving t_1=1. If j=1, the flag simply fixes t=1. The p-last condition changes only the fibre, not the underlying uniform-subset law. Therefore the conditional-average fibre is nonincreasing in t_j in each of the four flag measures, relative to its unweighted reference law.

## 4. Conditional tail law and Gold inequalities

Conditional on t_j=t, the remaining tail contains

    r=m-j a-elements, s=n-t+j c-elements.

The earlier prefix contains the first t-j c-elements, so it does not affect how many c-elements follow a_m. The tail remains a uniform shuffle after all allowed reweightings, because their fibres depend only on the prefix arrival vector.

The joint comparison table is triangular:

    [[A,0],[B,C]],

where A means zero c-elements after a_m, B means exactly one, and C means at least two. After normalization, its realized Gold deficits are

    1-2p_A, p_A-p_B, p_A-p_C.                              (7)

For the conditional r-by-s shuffle,

    p_A=r/(r+s),
    p_B=rs/[(r+s)(r+s-1)],

with p_B=0 when s=0. Both p_A-p_B and p_A-p_C are nondecreasing as s decreases. The first follows from r(r-1)/[(r+s)(r+s-1)], interpreting the one-letter case directly. For the second, assign independent continuous labels to all elements. Adding a c-element can only increase the number of c-labels above the largest a-label; the function taking values 1,0,-1 for counts zero, one, and at least two is decreasing.

Since s<=n, both lower deficits are bounded by their values at t=j, namely

    p_A-p_B >= r(r-1)/E,
    p_A-p_C >= [r(r-1)+(r+1)n-n^2]/E=-Q_j/E.               (8)

These inequalities do not need decreasing fibres; arbitrary mixtures over prefix vectors preserve them.

For the upper deficit, the conditional first-answer probability

    f(t)=r/(T-t)

is increasing in t. Each averaged restriction fibre w(t) is nonincreasing by Section 3. Reverse covariance therefore gives

    E[w f]/E[w] <= E[f].                                   (9)

One proof takes independent copies t,t' and observes that (w(t)-w(t'))(f(t)-f(t'))<=0 pointwise. The expectation is twice the covariance.

For the unflagged or p-last reference law, E[f]=m/T. For the a_1-first reference law, E[f]=(m-1)/(T-1)<=m/T. Hence every flag measure has

    1-2p_A >= (n-m)/T.                                    (10)

Equations (8) and (10), followed by convexity under normalized flag mixing, prove (4). The conditional probability f(t) is at least r/(r+n), and (9) supplies the upper endpoint, proving (5). Since the middle candidate in (3) is nonnegative, (6) is sufficient.

Finally, Q_j<=0 forces n<=2r: the positive root is at most 2r because Q_j(2r)=r(r-1)>=0 and the quadratic is increasing beyond its positive root. Thus r/(r+n)>=1/3.

## 5. Sharpness and genuine finite obstructions outside the region

Choose W_i empty. Then the side poset is independent of the backbone, and the unflagged shuffle law is uniform. The first-answer probability is exactly m/T, so the upper deficit (n-m)/T is attained. In particular m>n gives a genuine finite failure of this fixed strategy.

For the lower extrema, choose D to be a chain of length L, entirely above a_j. Thus W_1=...=W_j=D, and p is the top of D. Let N=T+|H|. For a backbone extension with t_j=t, the restriction fibre is exactly

    w_L(t)=binom(N+L-t,L).                                 (11)

Indeed no side element can occur before a_j, and after it one freely shuffles the remaining N-t backbone elements with the L side elements. Consecutive weights obey

    w_L(t+1)/w_L(t)=(N-t)/(N+L-t) -> 0.

As L tends to infinity, the induced law concentrates on t_j=j. Both lower deficits in (8), and the lower probability endpoint in (5), are approached. This proves that (3) is the exact universal infimum, not merely a sufficient estimate.

If Q_j>0, this limit has negative A-C. Therefore some finite L already defeats the fixed strategy. This establishes necessity in (6) using ordinary finite posets, without relying on an unrealizable abstract measure.

For a concrete example, m=n=4,j=2 has Q_j=2. A side chain of length 17 above a_2 gives

    [[779247,0],[602262,781641]].

The path in the lower-right cell has Gold deficit A-C=-2394. Every smaller positive side-chain length succeeds for this particular strategy; that finite first-failure claim is recorded by exact counts. No counterexample to the Gold conjecture is claimed.

## 6. Exact feasible prefix length and a Fibonacci-square boundary

For fixed m,j, some n>=m satisfies (6) exactly when

    j(m-1) <= (m-j)^2.                                    (12)

To see this, Q_j is increasing in n for n>=m, and

    Q_j(m,m)=j(m-1)-(m-j)^2.

Equivalently the largest admissible integer prefix length is

    floor([3m-1-sqrt(5m^2-6m+1)]/2).                       (13)

The asymptotic maximal fraction j/m is phi^(-2). More generally, if j/m tends to lambda and n/m tends to q, the limiting admissible band is

    1 <= q <= phi(1-lambda).

The curved boundary is exactly the old Pell curve in an effective width. Set

    u=m-j+1.

Then Q_j(m,n)=Q(u,n), where Q(u,n)=n^2-un-(u-1)(u-2) from the frozen fork theorem. Its positive integer equality pairs are

    (u,n)=(2,2),(7,10),(41,65),(274,442),...

For each such pair, every 1<=j<=n-u+1 gives a prefix-boundary member with physical chain length m=u+j-1<=n. This is a pullback of the same equality curve; different j generally define different posets, not identical dynamics.

The maximal-prefix diagonal cases have an especially simple Fibonacci form. Let

    d=F_(2k+1), b=F_(2k+2), e=F_(2k+3)=d+b.

Take

    m=n=d e, j=d^2.

Cassini gives b^2=d e-1. Therefore

    (m-j)^2=d^2 b^2=j(m-1),

so (12) is exact equality. The first pairs (m,j) are

    (2,1),(10,4),(65,25),(442,169),(3026,1156),...

and j/m=d/e tends to phi^(-2). Thus the prefix-size threshold is attained along an explicit Fibonacci-square sequence, rather than merely approached through an unspecified real asymptote.

There is an elementary reason the squares are forced. In any integer saturation case j(m-1)=(m-j)^2, write

    j=s a^2, m-1=s b^2

with common squarefree s>=1. Since m-j>0, it equals sab. Now m=s(a^2+ab), while m-1=sb^2, so

    s(a^2+ab-b^2)=1.

Hence s=1: both j and m-1 must be perfect squares, and

    j=a^2, m-1=b^2, m-j=ab, a^2+ab-b^2=1.                 (14)

Equivalently, the positive rank-one integer matrix

    [[j,m-j],[m-j,m-1]]

is exactly the outer product (a,b)(a,b)^T; the adjacent-integer normalization eliminates every nontrivial squarefree factor.

These unit-norm solutions are precisely (a,b)=(F_(2k+1),F_(2k+2)). A short descent proves this directly. Apart from (1,1), the norm equation gives a<b<2a and, for a>=2, b>=3a/2. Therefore

    (a,b) -> (2a-b,b-a)

preserves the norm, retains positive ordered coordinates, and strictly decreases a. It terminates at (1,1); reversing the map gives (a,b)->(a+b,a+2b), the square of the Fibonacci matrix. Consequently the displayed Fibonacci-square sequence exhausts all integer diagonal saturation cases, not just an example subsequence.

This structural reduction connects to the earlier Pell/Fibonacci theorem through u. No finite-field operation or Collatz dynamics is needed to establish the extension-count inequalities.

## 7. Independent exact controls and limits of the result

`gold_prefix_side_poset.py` checks:

- 2,074 actual nested-upset side-poset/cap configurations, including j=2 and j=3, by ideal-DP extension counts of all four flag tables;
- 894 individual interleaving-fibre counts and their conditional averages, independently checking the monotonicity argument with and without each flag;
- 35 actual side-chain concentration examples, including up to 107 vertices;
- all 4,950 pairs 2<=m<=100, 1<=j<m against the exact prefix-feasibility condition;
- explicit finite failures outside the universal parameter region and the Fibonacci-square equality family.

Normal and `python3 -O` outputs agree. The published JSON receipt is in `05-knowledge/results`; the script prints a fresh receipt without modifying files. The proof is the injection, subset coupling, conditional tail law, covariance inequality, and explicit concentration construction; finite checks are controls on these arguments.

The geometry still has a substantive restriction: every cross-relation points from a permitted chain prefix into an upset of the side poset. Relations in the opposite direction, restrictions involving the c-chain, or attachments beyond the permitted prefix need not produce decreasing averaged fibres. The theorem makes no claim that an arbitrary poset can be decomposed in this way.

Primary Gold-conjecture reference: Marcin Peczarski, *Comments on the Golden Partition Conjecture*, Contributions to Discrete Mathematics 12 (2017), 106–109:
https://cdm.ucalgary.ca/article/download/62349/46785/178291


Reproduce from the repository root with `python3 05-knowledge/experiments/gold_prefix_side_poset.py`. Helper scripts are colocated in the experiment directory.
