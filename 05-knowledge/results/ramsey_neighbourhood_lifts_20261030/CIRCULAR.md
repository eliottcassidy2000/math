# Uniform carousel rounding and the exact U5-free capacity

2026-10-30. **PROVED**, using the explicitly **CITED** structural theorem
of Gaku Liu. Packet-local labels C1-C3; no literature-priority claim.
This closes the two open extension targets in
[R2/R3, ramsey_circuit_interfaces_20261030/RAMSEY.md](../ramsey_circuit_interfaces_20261030/RAMSEY.md).

## C1. A chart for every odd circle

Let n=2r+1, r>=1, and let T_n have i->j when 0<j-i<=r modulo n.
Its maximal transitive supports are exactly the n intervals of r+1
consecutive vertices. Indeed, a transitive set has a source v; every
other member lies in the r consecutive outneighbors of v. That whole
interval is transitive. Thus maximality gives the asserted family.

For EVERY subset O and epsilon in{0,1}, there exists c:O->{0,1,2} such
that each of these intervals S satisfies, for j=0,1,2,

    |S intersect c^-1(j)| <= ceil((|S intersect O|-epsilon)/2)
                              + epsilon*[j=2].                    (C1)

This assertion covers every size, not just the checked n<=15.

**Even-budget chart.** If epsilon=0, properly color the cyclic order of
the m marked vertices with colors0,1,2: alternate0,1 when m is even;
when m>=3 is odd, alternate on the first m-1 vertices and give the last
vertex color2. For m<=1 any color works. An interval meets O in a
consecutive cyclic block. No color repeats consecutively, so its count
in a block of length s is at most ceil(s/2), including the whole circle.

**Odd-budget geometry.** Call a marked vertex forced if an interval S
meets O only there. It must have color2 in (C1). For m>=2, let g_i be
the positive clockwise gaps between successive marks, so sum g_i=n.
Mark i is forced exactly when

    g_(i-1)+g_i >= r+2.                                           (C2)

The available run excluding its two marked neighbors has length
g_(i-1)+g_i-1. It can contain an r+1 interval through i exactly under
(C2). For m>=4, two nonadjacent forced marks would use disjoint pairs
of gaps whose total is at least2r+4>n. Consequently all forced marks
are pairwise adjacent in a cycle of length at least4: there are at most
two, and if there are two they are consecutive. This is the key
all-size mechanism behind the earlier finite chart banks.

**Odd-budget chart.** Write H for color2 and L for a low color0 or1.
If m<=2 use H everywhere. If m=3 and there is an unforced mark, give
it color0 and the other two H. If all three are forced, use H,H,H.
In that case no interval contains all three: if it did, the middle
mark would have its two adjacent gaps summing to at most r, contradicting
(C2). These choices satisfy (C1).

For m>=4 rotate the list so the forced marks occupy positions0 and1
when there are two, or position0 when there is one. Use

    m odd:   H H L H L H L ... L,
    m even:  H H L L H L H L ... L,

and alternate the low labels0,1 in their order of occurrence. All forced
marks are H. Each individual low color has cyclic distance at least3
between successive appearances. In any block of s>=2 marks its count
is therefore at most ceil(s/3)<=floor(s/2). A block consisting of one
low mark cannot be an interval intersection, by the definition of forced.

For H, every block has H-count at most floor(s/2)+1. Here is a direct
check including wraparound. Assign H the increment+1 and L the
increment-1. In the even word the prefix sums range from0 to2 and the
total is0. In the odd word, after the initial0 all prefix sums are1
or2 and the total is1. An ordinary block has increment at most2; a
wrapping proper block also has increment at most2 (its removed interior
block has increment at least-1). The whole odd word has increment1.
Since the increment has the parity of s, this is exactly the required
H-count bound. This proves (C1).

The construction is linear in n once O is supplied. The finite checker
also derives forced marks independently from literal intervals, and
recovers maximal supports from the actual cyclic-triangle relation.
It chooses a starting mark in a supplied circular frame. This existence
construction is not asserted to be a canonical selector commuting with
every automorphism. If its individual layers are reused under unnamed
operations, the selected frame and induced color permutation must be
transported; the capacity inequality itself does not need those choices.

## C2. All integer attachments on every carousel

Let

    G(k)=sum_(0<=j<k) 2^popcount(j),  G(0)=0.
    G(2h)=3G(h),  G(2h+1)=2G(h)+G(h+1).

For any nonnegative integer weights a_i whose sum on each transitive
support of T_n is at most k,

    sum_i G(a_i) <= G(k).                                         (C3)

Apply C1 to the odd coordinates. Writing a_i=2b_i+e_i and
k=2h+epsilon produces the three layers

    a_i^(j)=b_i+[e_i=1 and c(i)=j],

of budgets h,h,h+epsilon. Coordinatewise the binary recurrence gives
G(a_i)=sum_j G(a_i^(j)). Induction on k proves (C3), with k=0,1
immediate since every pair is transitive. This is exactly the compiler
proved in inherited R2, now with a uniform rather than finite bank.
Zero extension covers induced subtemplates too.

Two consequences needed below are G(x)+G(y)<=G(x+y), and

    G(x)+G(y)+G(z)<=G(max(x+y,x+z,y+z)).                            (C4)

The latter is C2 for T3; set z=0 for the former. These are inequalities
for all integer arguments, not convexity assumptions about G.

## C3. Exact capacity for ALL U5-free tournaments

Let U5 be T5 with the arc0->1 reversed. For every integer k>=0,

    max{|T|: T is U5-free and tau(T)<=k} = G(k),                   (C5)

where tau is the largest transitive subtournament order, and the empty
tournament is admitted at k=0. Equivalently, the transitive-k Ramsey
threshold restricted to U5-free tournaments is G(k-1)+1.

**Cited input.** Liu, *Structure theorem for U5-free tournaments*,
[1208.0398v2, Theorem1.1](https://arxiv.org/pdf/1208.0398v2): a prime
U5-free tournament is an odd carousel or has a partition X,Y,Z such
that each of X union Y, Y union Z, Z union X is transitive. U5 is prime,
so the U5-free class is closed under substitution into U5-free quotients.
Liu already proves the sharp power-law exponent log_3(2) in Theorem6.1.
We refine the integer budget; we do not claim that exponent as new.

**Weighted prime step.** A carousel obeys C2. In the other case put
A=sum_X a_i, B=sum_Y a_i, C=sum_Z a_i. Each pair sum is at most k.
Superadditivity and (C4) yield

    sum_i G(a_i) <= G(A)+G(B)+G(C) <= G(k).

A transitive quotient obeys superadditivity. Apply these inequalities
up a prime/transitive substitution decomposition. Every child is
U5-free and smaller; its order is at most G of its transitive budget
by induction. The inherited exact substitution law is

    tau(Q[T_i])=max_(S transitive in Q) sum_(i in S) tau(T_i).

Thus the actual child budgets meet the required inequalities, proving
the upper bound. Notice that this argument does not require a successful
parity chart for every prime template in Liu's second case.

**Attainment.** Inherited R3 gives the literal Pascal-cell tournament
P_k on (j,t) with 0<=j<k and t a bit-submask of j. At the least
significant differing bit-pair use the cycle (0,0)->(1,0)->(1,1)->(0,0).
Its cyclic children have budgets (h,h,h) or (h+1,h,h), so its order is
G(k) and tau(P_k)=k. Iterated C3 substitution is U5-free because U5 is
prime. This proves equality for every k. The same literal construction
and its full symmetry recurrences survive the extension of the upper
bound's domain; no new isomorphism-class enumeration follows merely
from enlarging that domain.

**Algorithmic consequence.** Given k in binary, keep (G(m),G(m+1)).
Bit0 sends (a,b) to (3a,2a+b); bit1 sends it to (2a+b,3b). The exact
capacity takes O(log(k+1)) arithmetic operations. Finding a largest
transitive set in an arbitrary input graph is a different computational
problem; no unrestricted polynomial-time claim is made.

**Scope/hostiles.** Paley7 has tau=3 but order7>G(3)=5 and contains U5.
It refutes dropping the forbidden-subtournament hypothesis. The existing
twelve-vertex81/27 pair still has equal labelled cyclic triples and
different automorphism groups. Their common transitive-support observer
is adequate for these Ramsey queries but cannot classify full transport,
actual complements, or arbitrary re-marking. No A000568 term or universal
attachment-carrier theorem is asserted here.

**Generalization question.** Which other prime template classes admit
either uniform incidence rounding or a partition reducing their integer
capacity to a known small quotient? A failed chart is only a failure of
this sufficient proof method, not automatically a failed capacity law.
