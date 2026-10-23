# Independent addendum: silent square resets are Pell events

Status: **PASS.** This audits the coordinator's additional all-scale arrival
mechanism in `session32-root/PELL_ARRIVAL_CLOCKS.md` and its exact-integer
`check_pell_arrivals.py`/JSON controls. The conclusion concerns the square channels in an arrival window,
not Hamiltonian-path existence.

Write

    d(n)=floor(sqrt(2n−1))−floor(sqrt(n)),  n>=1.

For n>=2, the first floor can jump only when a square lies in
(2n−3,2n−1], an interval of length two. Consecutive positive squares differ
by at least three, and the interval starts at least at one. Thus this floor
jumps by either zero or one, never more. The second floor jumps by one
exactly when n is a perfect square. Consequently d(n) changes by one of
−1,0,1, and a negative change can occur only at a square n.

At n=k² with k>=2, the lower-floor jump is one. It is cancelled exactly
when an integer m satisfies

    2k²−3 < m² <= 2k²−1,

or, equivalently,

    m²−2k² = −1  or  −2.

Both interval endpoints and both integer possibilities matter. Equivalently,
the upper-square channel arrives at floor(m²/2)+1=ceil((m²+1)/2)=k².

The −1 branch has odd m and the negative Pell solutions

    m+k sqrt(2) = (1+sqrt(2))(3+2sqrt(2))^r,  r>=0.

Its k sequence begins 1,5,29,169,985,5741. The −2 branch has even m=2j,
so it is the ordinary positive Pell equation k²−2j²=1, with

    k+j sqrt(2) = (3+2sqrt(2))^r,  r>=0.

Its k sequence begins 1,3,17,99,577,3363. These standard Pell descriptions
are complete, not just numerical patterns. For example, multiplication by
3−2sqrt(2) maps every positive solution beyond the minimal one to a smaller
nonnegative solution: in the −1 case
(m,k) maps to (3m−4k,3k−2m), and in the +1 case
(k,j) maps to (3k−4j,3j−2k). Positivity and strict decrease follow from
the Pell equations, with the exceptional minimal solutions handled directly.
Iteration therefore reaches (m,k)=(1,1) or (k,j)=(1,0), respectively,
and reversing the operation yields precisely the displayed sequences.

The appearances of k=1 belong to n=1, for which the prior arrival-degree
formula d(0) is not defined. They are excluded from the difference claim.
In the −2 branch its m=0 is also not an entering positive square. The first
actual cancellations are therefore

    n=9,25,289,841,9801,28561,... .

At a silent square reset, the channel k² leaves the lower end of the
arrival-square window while m² enters its upper end. All common square
channels remain, but their actual neighbour labels decrease by one as n
increases. Thus cancellation preserves the number of available neighbours,
not their identities.

For n=25 specifically,

    arrival 24: square channels {25,36}, neighbours {1,12};
    arrival 25: square channels {36,49}, neighbours {11,24}.

The degree remains two. The new 49 channel supplies the endpoint-aware
11–25–24 ear once nontraceability of G24 is separately known. Pell
cancellation alone does not establish that obstruction or that ear's
Hamiltonian consequence.

The accompanying `audit_pell_resets.py` uses explicit integer checks for all
square roots k=2,...,10000 and all n=2,...,100000, verifies both parity
branches, reconstructs the channel exchange, and compares the complete
bounded cancellation set to the Pell generators. Normal and optimized
certificates are identical. The interval and descent arguments above, not
the finite bounds, establish the all-order statement.

## Density coefficient and final source check

The root note's asymptotic count is correct. Put lambda=3+2 sqrt(2).
In each Pell family, k is a positive constant times lambda^r plus a term
of order lambda^(−r), so n=k² is a positive constant times lambda^(2r)
to a bounded multiplicative correction. Each family therefore contributes
log(X)/(2 log(lambda))+O(1) events below X. Beyond k=1 the two families
are disjoint: an overlap would make 2k²−1 and 2k²−2 distinct positive
squares differing by one, which is impossible. Consequently

    S(X)=log(X)/log(lambda)+O(1).

In particular S(X)/sqrt(X) tends to zero, proving zero density even among
square birth times. The root note correctly distinguishes the real field
Q(sqrt(2)) in this norm equation from the imaginary field Q(sqrt(−2));
no transfer of unrelated CM or representation data is warranted.

The audit certificate pins the root prose, program and finite receipt and
compares every root Pell witness with the independently generated set.

## Complete arrival-level polynomials

The appended all-degree interval law in the root note also passes. Fix an
integer j>=0 and r=floor(sqrt(n))>=1. The other floor equals r+j precisely
when (r+j)^2 <= 2n−1 < (r+j+1)^2. Combining these inequalities with
r² <= n <= (r+1)²−1 yields exactly

    L=max(r², ceil(((r+j)²+1)/2)),
    U=min((r+1)²−1, floor((r+j+1)²/2)).

No floor/ceiling endpoint correction is needed. Distinct r bands are disjoint.
A contributing band has 2r² <= (r+j+1)², so
r <= (1+sqrt(2))(j+1). The exact integral bound is
(j+1)+isqrt(2(j+1)²). This proves finiteness of every level set and therefore

    D_j(z)=sum_(nonempty bands) z^L(1−z^(U−L+1))/(1−z).

The apparent singularity at z=1 is removable, with value U−L+1 per band.
The j=0 case is included: its intervals are [1,2] and [4,4], giving three
isolated births. The counts for j=0,...,9 are
3,12,24,34,47,59,68,84,92,105. They are not given by 12j beyond j=2.

For the producer's rows j=0,...,100, this audit checks every supplied interval
at both endpoints, checks maximality inside its fixed-r band, and certifies
each omitted r band by its endpoint degrees. On a fixed-r band, d(n) is
nondecreasing, so these tests prove the exact interval, not merely a sample.
This avoids repeating the producer's full vertex enumeration. The all-order
floor equivalence, not those finite controls, proves the general formula.
The root polynomial script and JSON receipt are included in the source pins.
