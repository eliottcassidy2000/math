# Niche recurrences: source laws, arithmetic boundaries, and composition

2026-10-29. Source pin `d3741fe6ed02f6072dbe6bfda2f3d04193d75ce0`.
This session recovers older mechanisms rather than restating recent Lucas or
A000568 formulas. All claims below have explicit operation menus, proofs and
exact computational controls. Classical tools are credited; no external
priority or major-open-problem resolution is claimed.

## 1. Five recovered routes and their new recurrences

### A. September floor reciprocity closes an implicit boundary

[Arithmetic boundary recurrences](arithmetic/ZERO_PRODUCT_RECURRENCES.md)
start from the old product-floor identity. Its correction counts tuples of
nonzero residues whose product is zero modulo n. For composite n, that
correction has exact minimal scalar recurrence order

    2^omega(n) + 2^(omega(n)-1) [Omega(n)-omega(n)].

The proof keeps the nonzero source dependence under CRT by first counting
all residues, then subtracting tuples containing a globally zero coordinate.
Prime-power depth becomes repeated spectral roots; distinct prime supports
become distinct roots. The floor sums themselves replace the leading root
n-1 by n(n-1)/2. This gives precise different formulas at11,12,13.

The arithmetic carrier for a *named* next multiplier still needs all tau(n)
divisor states. Its dimension first exceeds the uniform scalar recurrence
order at **36=6²**: nine states versus order eight. Two explicit positive
source distributions agree after every number of uniform random future
multiplications, but multiplication by2 separates them by probability1/4.
An exact eigenvector identity certifies agreement for all averaging times,
not just a computed prefix. Conversely, eight specified one-step multiplier
probes plus known total mass recover all nine state masses by divisor Möbius
inversion. This probe count is optimal among fixed linear banks on arbitrary
normalized masses, not a lower bound for adaptive single-state identification.

The reason36 is first is structural: a gap requires two different primes
both to occur with repeated powers. The formula is not a numerical bridge
from6 to every other six-object structure.

### B. Older rational integrability becomes recurrence-memory growth

[Pole recurrence memory](jets/POLE_RECURRENCE_MEMORY.md) recovers the
October19 metabelian two-flow classification and an older endpoint-Prony
mechanism. Translating a finite signed pole source repeatedly produces a
normalized moment sequence whose exact all-degree recurrence order satisfies

    p_n = q n + b

past an explicit source-span threshold. Here q is the number of translation
orbits carrying nonzero total source weight. In the inherited fourth-order
pole family, bounded recurrence memory is therefore equivalent to rational
integrability. The signed boundary formula gives b exactly.

A fixed moment degree r simultaneously satisfies a short *time* recurrence,
Delta_n^(r+2) m_r(n)=0. These are different axes. Every prescribed finite jet
can be blind for every iterate while the full recurrence rank grows linearly.
Conversely, knowing a support bound B makes the first2B moments sufficient,
and that uniform bound is sharp.

A finite cyclic variant keeps another exact clock: if coefficients and
spatial sites lie in Z/M, the source has least spatial period d and period
sum S_d, its cumulative translation clock is dM/gcd(M,S_d). Spatial period
alone misses the coefficient accumulation.

### C. Old palette gluing produces source-mass and entropy recurrences

[Orbit-mass recurrences](attachments/ORBIT_MASS_RECURRENCES.md) recover the
old named-output palette contraction, whose24-by24 fibres gave576. For iid
leaf colors in the ternary wreath action, let S_r(k) be the sum of the rth
powers of actual orbit probabilities. Then

    S_r(k+1)=3^(r-1) S_r(k)^3-[3^(r-1)-1] S_(3r)(k).

This single hierarchy includes supported class counts, normalization and
collision probabilities. A divisor–Möbius formula handles general cyclic
root actions. A single requested moment at depth k uses O(k²) scalar
arithmetic, with exponentially growing bit lengths still charged. The
centralizer recurrence retains one paired product palette, not two
independent orbit quotients.

The Shannon recurrence subtracts the expected log orbit size. This is an
expected information saving, not a fixed-width representation theorem.
For heterogeneous child sources, all nine pairwise overlaps can agree while
assembled collision probabilities differ:1534/16^6 versus1536/16^6.
A six-way shared-label term is the exact missing coordinate.

### D. Circulation reciprocity lifts to a closed boundary grammar

[Schur boundary recurrences](flows/SCHUR_BOUNDARY_RECURRENCES.md) recover
regular-tournament circulation cones. A scalar child circulation series is
not closed under cyclic substitution. Retain instead the matrix

    M_T(lambda,mu)=CT F_T(q;x) s_lambda(x) s_mu(x^-1).

Cauchy expansion gives

    R_(C3[T,T,T]) = sum_t q^(3t) tr(M_T[t]^3).

A Littlewood–Richardson recurrence generates the matrix entries themselves,
so this is a closed compositional rule rather than an unexplained sidecar.
It extends to any specified quotient tournament. Fixed flow-degree queries
require finitely many boundary partitions; no uniform finite matrix bank is
claimed.

The same method gives exact3/4/5-cycle counts in the recursive family:
at27 vertices these are819,8991,70470. Root rotation has its own trace
recurrence with q replaced by q³, preserving actual orbit clocks. This is
not the full automorphism quotient or a Hamilton-path count.

There is an actual composition hostile, not just a missing formal entry.
The five-vertex child C3 ordered-joined with T2 and its converse have the
same size and entire scalar flow series1/(1-q³). Attaching each in the same
cyclic context C3[child,K1,T2] gives matching coefficients through degree11,
then **9470 versus9462 at degree12**. A joint size-three boundary partition
is the first difference detected in this context. Independent enumeration
of cross-block matrices verifies the discrepancy without Schur/LR code.

The old reciprocity shift also has a sharp boundary: on C3 the raw series
has shift3, but its one-port Schur response has shift1. Charged boundary
responses cannot inherit the zero-divergence inversion law unchanged.

### E. Periodic valuation clocks compile every odd critical arity

[Residue-clock branching](valuation/RESIDUE_CLOCK_BRANCHING.md) extends
Session94's critical ternary law while recovering the older critical-run
extractor lesson: preserve source probabilities and unconsumed tails.
For iid accelerated valuations P(a=j)=2^-j and an integer d>=2, a periodic
residue test can give d offspring with probability1/d exactly when d is odd. Its least residue
clock is ord_d(2). Binary digits of (2^h-1)/d construct the unique selector.
The singleton cases are3,7,15; arity5 uses residues3,4 modulo4.

A bivariate generating recursion keeps both tree size and native division
clock, with a closed coefficient formula. Conditioning on tree size removes
a precisely computed between-type variance but retains a boundary correction.
The arithmetic source, residue labels, residual geometric waits and original
source order remain necessary for reconstruction.

Sharing a source across arities does not make their branching laws
independent. Arity3 and7 offspring have covariance-2/3 even though their
residue periods2 and3 are coprime. CRT address compatibility is not product
source measure. Auxiliary extinction has density one in the stated model,
but the inherited41 ->31 ->47 ->71 ->107 control still blocks any inference
of native Collatz descent.

## 2. A unifying statement: distinguish the recurrence axis from the operation menu

The common mechanism can be stated precisely without claiming that these
objects are isomorphic.

Let V be a finite-dimensional vector space, let L:V->K be a linear
observation, and let A be the unital algebra generated by the allowed linear future
operations. Two states are observationally equivalent exactly when their
difference lies in

    N_A = intersection_(B in A) ker(L B).

This kernel is invariant under every allowed operation: if v lies in N_A,
then L B C v=0 for every B,C in A. The quotient V/N_A therefore carries
all permitted future observations. Any exact linear encoding E must have
ker(E) contained in N_A, so V/N_A is a linear image of E(V) and its
dimension is a lower bound for that encoding. For the smaller menu
consisting of powers of a single T, the kernel is only

    N_T = intersection_(r>=0) ker(L T^r),

which can be strictly larger. The scalar output sequence obeys the minimal
polynomial of the induced observable T; its recurrence order does not by
itself recover the richer A-carrier. This is classical linear realization
logic, applied here with explicit new source examples and recurrences.

The modulo36 eigenmode is a concrete strict inclusion. The pole example
shows another boundary: allowing all time steps on a fixed small jet still
differs from allowing every jet degree. The Schur calculation identifies a
restricted attachment menu that closes exactly; arbitrary named arcs require
more coordinates. Orbit-mass powers instead require a nonlinear hierarchy
of replicated sources, visible in the r ->3r coupling and the six-way
hostile. The valuation tree uses an infinite stopping grammar and therefore
is not asserted to be an instance of the finite-dimensional theorem.

Thus a productive research procedure is:

1. Specify which index is being recurred: arity, depth, moment degree, time,
   flow weight, or physical clock.
2. Specify which source law and future operations the recurrence preserves.
3. Derive its closed carrier, including boundary terms and joint labels.
4. Find an observation outside that menu and test the cheapest hostile.

This produces usable formulas while exposing exactly when a seemingly free
bit, scalar mode, port, or clock has merely been hidden.

## 3. Reproduction and scope

Each lane contains a standard-library Python checker. [Independent audit](audit/AUDIT.md)
records separate implementations and hostile controls. Run all supplied
scripts under ordinary Python and `python -O`; checks remain active in both.
Producer and audit scripts are not a census of all tournaments or all
arithmetic trajectories. Final publication verification replays the complete
packet from the committed checkout and checks relative links and file bytes.

No new A000568 coefficient record, Lucas-prime infinitude, Collatz convergence,
LRC closure, or generic equivalence between major open problems is claimed.
The deliverables are exact recurrences, source-preserving constructions,
minimality results in stated models, and explicit failures of tempting
compressions.
