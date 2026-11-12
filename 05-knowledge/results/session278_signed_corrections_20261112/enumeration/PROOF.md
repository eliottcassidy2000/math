# Divisibility-defect gas: a signed correction calculus for A000568

**PROVED identities and inequalities; new packet, not a claim of a faster exact counter.**
Base: e07f05b9fee231e59ba7c474fd52c59bc63d3980. Ordinary unlabelled tournament
classes are counted, not Hamiltonian aggregates or labelled tournament mass.

## Inheritance and cross-problem route

Session277 `absolute/PROOF.md` supplies the positive odd-cycle Burnside formula
and the obstruction to rounding a proper positive-sector truncation. Its
`engine/PROOF.md` supplies interval/residue recovery. Both are inherited inputs.
A less direct source is `01-canon/theorems/THM-4346-lrc14-halfturn-current-brownian-kernel-and-cubic-moment-boundary.md`:
its Brownian `min(i,j)` kernel and warning about low-order moment blindness
suggest replacing a difficult arithmetic kernel by a nested kernel and retaining
the exact missing interactions. The present Burnside application is separate
from that LRC theorem. Bonferroni's principle is classical; the correction
objects, envelope, resummation and quantified obstruction below are the application.

## 1. Two Gram embeddings and a cone-restricted comparison

For an odd partition lambda=(a_1,...,a_r) of n, put

    z(lambda)= product_a a^(m_a) m_a!,
    E(lambda)= sum_i (a_i-1)/2 + sum_{i<j} gcd(a_i,a_j).

The inherited formula is a(n)=sum_lambda 2^E/z. Replace gcd by min:

    Ebar=sum_i (a_i-1)/2+sum_{i<j} min(a_i,a_j),
    W=2^Ebar/z,  d_ij=min(a_i,a_j)-gcd(a_i,a_j).

All d_ij are nonnegative integers. An edge ij exists precisely when neither
cycle length divides the other. Give it weight q_ij=1-2^(-d_ij). Then EXACTLY

    2^E/z = W product_edges (1-q_ij).                         (1)

Thus an arithmetic symmetry sector is an envelope times a finite repulsive
interaction partition function. Equal lengths and singleton cycles cause no
interaction. The envelope is exact for a partition iff its distinct lengths
form a divisibility chain.

Both kernels are Gram kernels:

    min(a,b)=sum_{t>=1} 1[t<=a]1[t<=b],
    gcd(a,b)=sum_{d>=1} phi(d)1[d|a]1[d|b].

The target replaces nested intervals by divisor incidence. The difference is
entrywise nonnegative, hence nonnegative on nonnegative occupancy vectors,
but NOT positive semidefinite: on lengths3,5 it is [[0,2],[2,0]]. Do not infer
an operator/covariance domination theorem from this valid multiplicity-cone
comparison. This is the exact sidecar required for the LRC-inspired analogy.

## 2. A collapsed upper-envelope recurrence

Process odd lengths a in increasing order. Let D_a(m) sum envelope weights
for partitions of m using lengths at most a. Start D_{-1}(0)=1, other entries0.
For each prior total m and j new cycles of length a, add

    D_{a-2}(m) 2^[jm+j(a-1)/2+a binom(j,2)]/(a^j j!)

into D_a(m+aj). Then U_n=D_n(n)>=a(n).

Proof: every old cycle has length at most a; its interaction with each new
cycle is its own length, so the total cross exponent is jm. This is why size
alone suffices for the envelope, whereas gcd requires arithmetic response
information. The recurrence is finite, rational and exact for U, not for a.
It uses O(n) coefficient storage and O(n^3) arithmetic operations (a loose
bound), with growing integer bit cost. No ordinary-class approximation is
silently relabelled exact.

## 3. Signed corrections, signs and sharp two-sided error

For a sector let e_j(q) be the elementary symmetric polynomial in its M
nonzero edge weights, e_0=1. Define

    C_j=sum_lambda W(lambda)e_j(q(lambda)),
    S_k=sum_{j=0}^k (-1)^j C_j.

Then for EVERY n,k:

    (-1)^k(S_k-a(n)) >= 0,                                   (2)
    sum_{lambda:M>k} W (k+1)/M e_{k+1}(q)
       <= |S_k-a(n)| <= C_{k+1}.                             (3)

These are bounds in the correct direction, not heuristic extrapolation.
They apply with the resummed baseline of section4 as well.

Proof: independently activate each edge with probability q_e and let X be the
number of active edges. E binom(X,j)=e_j. For integer X>=1,

    sum_{j=0}^k (-1)^j binom(X,j)=(-1)^k binom(X-1,k),

where the right side vanishes if X<=k. For X=0 both exact product and partial
sum equal1. The absolute error is therefore
E[binom(X-1,k)1_{X>k}]. On X>k this equals
((k+1)/X)binom(X,k+1). Since k+1<=X<=M, (3) follows.
The constants are sharp for arbitrary edge probabilities: upper equality at
M=k+1; lower equality when all M edges activate surely. For these tournament
weights q<1 the latter is a limiting universal constant, not an asserted
attained equality for a finite sector.

Partial sums of a fixed parity need NOT improve monotonically. Keep the maximum
of available lower bounds and the minimum of upper bounds, rather than assume
nested errors. A test initially made that monotonicity assumption and rejected
it before publication. Formulas(2)-(3), not monotonicity, drive the code.

At k>=max_lambda M the formula is exactly a(n). For a prescribed finite n this
is an all-order terminating correction algorithm. Computing the C_j here still
enumerates partitions; this is a structural calculus and certification method,
not a demonstrated asymptotic speed improvement over existing counters.

## 4. Resum the 3-cycle gas before expanding

Every edge from a 3-cycle to a cycle of length b contributes d=0 when3|b and
d=2 otherwise. With t 3-cycles and j_b b-cycles, all those edges contribute

    4^[-t sum_{b>=5,3 does not divide b} j_b].                 (4)

Multiply W by (4) EXACTLY and remove these edges from the expansion. Equivalently
use actual gcd whenever the smaller length is3. The residual graph contains
only cycles of lengths at least5. Apply(1)-(3) to that graph unchanged.
This is a genuine resummation of arbitrarily many signed terms into a single
positive factor. It is not merely keeping more small-support Burnside sectors.

The four terms S_3 are exact for every n<=23 in this model. If residual M>=4,
at least four cycle vertices are needed (three vertices have at most3 edges).
Four lengths>=5 with total<=23 must have at least three5's: otherwise two
lengths>=7 give total>=24. Three5's and a possible fourth length have at most3
incompatible edges. More vertices cost at least25. Contradiction. At n=24,
(5,5,7,7) has four incompatible edges, all q=15/16, so S_3 is not exact.
This is a sharp support threshold, not an empirical cutoff.

## 5. Fixed correction order still fails at large size

For any fixed k put r=k+1. The plain sector (1^f,3^r,5), n=f+3r+5,
has exactly r incompatible edges, all q=3/4. Its contribution to
|S_k-a(n)| is exactly

    (3/4)^r 2^[binom(f,2)+f(r+1)+r+2+3binom(r,2)+3r]
       /(f! 3^r r! 5).                                      (5)

All sector remainders have the same sign by(2), so(5) is a LOWER bound on the
whole error. It tends to infinity with f. For the 3-resummed model use instead
(1^f,5^r,7), with q=15/16 and lower bound

    (15/16)^r 2^[binom(f,2)+f(r+1)+2r+3+5binom(r,2)+5r]
       /(f! 5^r r! 7).                                      (6)

Therefore no fixed finite correction depth in either model supports eventual
nearest-integer recovery. Resummation can yield huge practical precision gains
without abolishing the absolute-error obstruction. Adaptive depth, a certified
residue, or further exact aggregate summation remains necessary.

## 6. Reach and obligations

The experiment includes direct permutation edge-orbit Burnside controls n<=7,
partition-vs-collapsed-envelope controls, both signed bounds, and sharp-threshold
checks. It does NOT classify Hamiltonian count support. Ordinary class mass and
an exact target fibre remain different objects. A plausible transfer to
Chan–Pak is to retain interaction graphs for physical edits, but that needs
its own realizability laws; no averaging-to-existence implication is made here.

Open algorithmic improvement: compute low C_j by a marked-size recurrence that
avoids enumerating full partitions, without dropping divisibility information.
Open structural improvement: choose/resum a larger divisibility-chain backbone
adaptively, certified by(3), and quantify its cost against saved residue bits.
