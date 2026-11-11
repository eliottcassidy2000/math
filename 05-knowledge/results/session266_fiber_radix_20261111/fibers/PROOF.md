# From multiplicity tables to an exact growth compiler

Status: inherited all-order SCC classification plus newly assembled cumulative,
automorphism-weighted and inverse-growth consequences. No literature-priority
claim. H always counts directed Hamiltonian **paths**, not Hamiltonian cycles.

## Inheritance and what the compiler actually completes

The closest proved mechanism is Session24's exact fixed-H strong-core series:
`05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md`.
The stronger finite inventory is opus S50:
`05-knowledge/results/h_spectrum_multiplicities_20261021.md` and THM-4744.
These already give eventual polynomials, self-converse quasipolynomials,
logarithmic cutoffs, the two-degree loss on 7 times powers of 3, and complete
strong inventories below 625. We do not claim those as new.

The hostile is H=39: padding its strong cores alone misses the 3*13 SCC words.
The second hostile is H=45: limiting core size to the census order misses
longer ordered products and produces the wrong polynomial degree.
The underused sidecar is automorphism order, which permits a genuinely different
labeled-frequency channel without confusing division by n! with isomorphism.

This packet turns the inherited finite-core description into an executable
all-order response for every odd h<=73, including cumulative counts. It adds
exact inversion of the observed growth table, a convexity law, and the weighted
channel. The cap 73 is an intentionally frozen reusable input, not a new census
or the current limit of knowledge. The existing below-625 input could be loaded
into exactly the same compiler. For arbitrary h the algorithm is finite, but
its strong-inventory calculation remains to be executed.

## 1. Four different frequencies

Let A_n(h) count ordinary isomorphism classes, S_n(h) self-converse ordinary
classes, M_n(h)=(A_n(h)+S_n(h))/2 classes after identifying converse, and L_n(h)
labeled tournaments on a specified n-element set. Put W_n(h)=L_n(h)/n!.
Thus W is a rational weighted count, not a count of ordinary classes.
An empty tournament at n=0,H=1 is used only as the convolution unit.
Cumulative user-facing counts run from n=1 to N.

An automorphism fixes a directed Hamiltonian path only if it fixes its ordered
vertex sequence pointwise. Hence its action on paths is free and |Aut T| divides
H(T). Distinct SCCs cannot be permuted by automorphisms, even when isomorphic:
their order is intrinsic. Consequently |Aut T| is the product of the component
automorphism orders, with no extra factorial for repeated components.

For each h>1 define strong inventories

    P_h(z)=sum_[C strong,H(C)=h] z^|C|,
    Q_h(z)=sum_[C strong,self-converse,H(C)=h] z^|C|,
    R_h(z)=sum_[C strong,H(C)=h] z^|C| / |Aut C|.

Camion's Hamilton cycle gives |C|<=h. The inherited Busch bound improves this
to |C|<=1+floor(3 log_5 h). These polynomials can therefore be computed by a
finite, canonically deduplicated tournament enumeration. This is an effective
classification, not a simple closed formula for all their coefficients.

## 2. Exact recurrence and cumulative coefficient formula

For either X=A with inventory P or X=W with inventory R, set X_0(1)=1,
X_n(1)=1 for n>=0, and zero at negative n. Let I_{m,d} be the relevant
strong-inventory coefficient. Then for h>1,

    X_n(h)=X_(n-1)(h)+sum_(d|h,d>1) sum_m I_(m,d) X_(n-m)(h/d).

The first SCC is either one singleton or one nontrivial strong component.
This is a triangular, exact recurrence in n and h. It applies to every n,
including the small orders before the eventual polynomial starts.

Equivalently, let B_(h,k,s) count ordered nontrivial strong core words with
product h, k letters and total size s. In the W channel weight each such word
by the product of inverse automorphism orders. In both channels,

    X_n(h)=sum_(k,s) B_(h,k,s) binom(n-s+k,k),
    sum_(n=0..N) X_n(h)=sum_(k,s) B_(h,k,s) binom(N-s+k+1,k+1).

Terms with N<s or n<s are ZERO, not polynomially extrapolated. For h>1 the
n=0 term vanishes, so the cumulative formula already answers sizes 1 through N.
For h=1 subtract the auxiliary empty object: A_<=N(1)=M_<=N(1)=N.
For labeled cumulative counts use sum n! W_n(h), not N! times the cumulative W.

If r(h) is the maximum core length and b_h=sum_s B_(h,r(h),s), then

    A_n(h) ~ b_h n^r(h)/r(h)!,
    A_<=N(h) ~ b_h N^(r(h)+1)/(r(h)+1)!.

For W replace b_h by the positive automorphism-weighted top-core mass.
Thus labeled counts are n! times an eventual polynomial; they are NOT themselves
eventual polynomials. Uniformly random labeled tournaments have this fixed-h
probability divided by 2^(n(n-1)/2), which tends to zero faster than exponentially.

The inherited safe polynomial onset for A_n is floor(3 log_5 h). A safe onset
for the cumulative polynomial is the same bound: in each cumulative binomial
it suffices that N>=s-k-1, and s-k<=3 log_5 h. The same statements hold for W.
No cancellation of the leading pole occurs because all core weights are positive.

## 3. A global convexity law and a recovery theorem

For every h>1 and every n>=0, with negative-index values zero,

    A_n(h)-2 A_(n-1)(h)+A_(n-2)(h) >= 0.

The same is true for W. This holds BEFORE polynomial onset. In generating
functions each core term contributes

    (1-z)^2 z^s/(1-z)^(k+1)=z^s/(1-z)^(k-1),

which has nonnegative coefficients since k>=1. Every fixed-h ordinary frequency
therefore grows with a nondecreasing increment. Its cumulative frequency has
nonnegative third backward difference. Converse-merged counts need not inherit
this convexity because the self-converse correction oscillates with parity:
for h=3, M_n(3)=floor((n-1)/2) at n>=3 has alternating increments 0 and 1.

Growth also recovers the lost strong inventory. Write F_h(z)=sum A_n(h)z^n.
The first-component recurrence implies

    P_h(z)=(1-z)^2 F_h(z)
           -(1-z) sum_(d|h,1<d<h) P_d(z) F_(h/d)(z).

Process h in increasing order. This is triangular divisor inversion: the
observed ordinary growth table determines every strong inventory uniquely.
The same formula recovers R from W. It distinguishes a structural inventory
from an unstructured histogram rather than silently conflating the two.
For a prime h, or any h with no factorization into two attainable values >1,
all correction terms vanish:

    [z^n]P_h = A_n(h)-2A_(n-1)(h)+A_(n-2)(h).

Thus a discrete acceleration at order n counts new strong classes exactly.
For composite h with a realizable factorization, acceleration contains a
positive background of composite cores; divisor inversion removes it.
The number of strong h-classes is the eventual slope only in the indecomposable
case. In particular h=63 is arithmetically composite but indecomposable in the
realizable monoid, whereas h=39 is not.

## 4. Converse channel and limits of a scalar histogram

The inherited exact center formula is

    J_h(z)=sum_(x^2|h) F_x(z^2) Qplus_(h/x^2)(z),
    Qplus_1=1+z, Qplus_d=Q_d for d>1.

Here F_x counts ordinary first halves, not converse-merged first halves.
Self-converse and merged cumulative counts follow by multiplying their series
by 1/(1-z), or simply summing their exact integer outputs.
Their eventual quasipolynomial period divides 2; cancellation can make it 1.
No global three-cycle or Collatz convergence claim is needed anywhere here.

An ordinary scalar histogram cannot recover the automorphism-weighted or
self-converse inventories without these sidecars. This is a useful obstruction
to over-interpreting a single H multiplicity. The compiler stores all three
strong channels separately and uses no unlicensed conversion between them.

## 5. Explicit all-order and finite outputs

The inherited complete strong H63 inventory is two converse-paired, asymmetric
classes at order 8; H7 and H21 are absent. Therefore for N>=8,

    A_N(63)=2(N-7), S_N(63)=0, M_N(63)=N-7,
    A_<=N(63)=(N-7)(N-6), M_<=N(63)=(N-7)(N-6)/2,
    L_N(63)=2(N-7) N!.

At N=30 this means 46 ordinary classes of order 30, but 552 across orders
1 through 30; 23 and 276 respectively after converse merging.
That distinction answers the two readings of “how many times up to this size.”
For H3, A_<=N(3)=(N-2)(N-1)/2 for N>=3 and L_n(3)=n!(n-2)/3.
For H1 there is one ordinary and one converse class at every positive order,
but n! labeled transitive tournaments.

`atoms.json` freezes inherited strong data through H73 with source hashes.
`build.py` computes arbitrary requested N and all odd h<=73. `report.json`
contains N=30, both ordinary/cumulative and converse/cumulative frequencies,
and exact labeled frequencies. No new large census was run. The 5,704 active
controls cover integer output, the global convexity, cumulative summation, and
inverse recovery. `python -O` gives the identical result; these constructor
checks do not purport to re-audit the inherited atom census.

## 6. What remains genuinely open

The finite reduction does not resolve the complete set of attainable odd h
(the Chan–Pak question), give an efficient closed formula for the strong atom
coefficients, or classify wiring within a metagraph H-fiber. It does completely
specify the multiplicity sequence once its finite strong inventory is supplied.
The clean next output is a persistent inventory service: add each newly computed
strong class once with its automorphism and converse data, then regenerate all
order and cumulative frequencies with no repeated tournament enumeration.

## Querying the frozen inventory

Run `python build.py --h 63 --n 100 --output h63_through100.json` for one
count, or omit `--h` for all odd counts through73. Uppercase `--H`/`--N`
aliases work too. Negative order and out-of-scope/even H are rejected instead
of silently returning an incomplete answer. Every cumulative field excludes
the auxiliary empty object. At H1 and N30 the ordinary and merged cumulative
answers are30, not31.

The three compact CSVs list ordinary, cumulative ordinary, and merged counts
at every size1 through12 and every odd H1 through73, explicitly including zero
rows for7 and21. The source bytes in `atoms.json` match the pinned repository
snapshot linked there. `receipt.json` records final replay and argument checks.
