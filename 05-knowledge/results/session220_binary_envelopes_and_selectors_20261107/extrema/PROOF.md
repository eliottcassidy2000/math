# Exact extremal envelopes for the three-core binary grammar

Status: PROVED algebraic extremal theorems; FINITE-EXACT verification of stated
bounded domains. The all-target three-core conjecture and Chan–Pak conclusion
remain OPEN. No literature-priority claim. Independent review pending.

## Inheritance and scope

The input is commit `acbad6c8636cf8cbb675b4c947fa19e33ac0df23`, specifically
`05-knowledge/results/session216_overnight_constructions_20261107/three-core/PROOF.md`
and its `CORE_STATES.txt`. We inherit its construction, strongness rule,
counting formula, and the definition of f. Its Section 6 is an OPEN sufficient
conjecture, not a theorem. Session217's
`exact_inverse_and_gate_limits_20261107/inverse/SUFFIX_INVERSE.md` (with the
`session217_` path prefix) gives complete square-root-time inversion, not
occupancy. Its `run-entropy/PROOF.md` shows that bounded-run recipes have
zero-density support. Session216's numeration proof rules out specified
universal radix rewrites and supplies the hostile balanced-band gaps 265,287.
The permanent restricted omissions 47,63 are further hostile controls.

Concept board: exact occupancy; binary rounding; submodular cut energies;
low-bit inverse cylinders; run entropy. The least-used sidecar here is the
order relation between a word, its reversal, and its quadratic f, rather than
only the scalar output. This yields exact extremal envelopes, not a new
physical operation or an induction on represented targets.

## 1. Binary rounding lemma

For n>=0 define f(0)=0,

    f(2u)=f(u),   f(2u+1)=f(u)+u+(u mod 2).

Write rev_k(v) for reversal of exactly k bits, and first_k(v) for their first
bit, with both zero for k=0. The selected-pair definition of f, separating
pairs internal to either word and crossing their boundary, proves

    f(2^k x+v) = f(x)+f(v)+x rev_k(v)+(x mod 2) first_k(v),
    0<=v<2^k.                                                    (1)

Suppose an internal zero has a selected bit after it, so that

    n=(2U)2^k+V,   U>=1, k>=1, 1<=V<2^k.

Round upward to n'=(2U+1)2^k. Then

    n'>n,       f(n')<=f(n),       rev_m(n')<rev_m(n),              (2)

for any common word length m, and the first bit stays unchanged. Also n' is
even. In fact (1) gives

    f(n) = f(U)+f(V)+2U rev_k(V) >= f(U)+2U,
    f(n') = f(U)+U+(U mod 2) <= f(U)+2U.

For reversal, the lowest changed position is a selected bit of V that is
cleared; it becomes the highest differing position after reversal, proving
strict decrease. These statements include equality in f: for n=5 and n'=6,
both have f=2. The two strict order comparisons still hold.

Repeated rounding terminates because n strictly increases within its fixed
length. The resulting nonzero word has one run of ones and trailing zeroes.
Leading zeroes are unchanged by internal rounding.

## 2. Exact prefix maximum, and constant-size evaluation

Let K,h>0, S>=0, and m>=1. Put

    g(P)=KP-h f(P)-S(P mod 2).

Over 0<=P<2^m, its exact maximum is

    max(0, max_{0<=j<m}
        [K(2^m-2^j)-h(2^(m-j)-2)-S 1_{j=0}]).                  (3)

Indeed each rounding step strictly increases g by (2). A remaining one-run
word not using the highest bit can be shifted left until it does; this
increases P, preserves f, and cannot increase the parity penalty. Therefore
every nonzero maximizer is P=2^m-2^j. Such a word has
f(P)=2^(m-j)-2, proving (3). If the domain is the positive m-bit words,
2^(m-1)<=P<2^m, simply omit the zero candidate. No assumption K>h is needed
for this maximum theorem. It applies also with h=0, but h>0 is the stated
case relevant here.

For j>=1, maximize by minimizing

    K 2^j+h 2^(m-j).

This strictly convex real function is minimized at
j*=(m+log_2(h/K))/2. Only the two neighboring integers, clipped to [1,m-1],
need testing, in addition to j=0 and P=0. This is an at-most-five-candidate
bound (in fact at most four distinct candidates); it is a mathematical
observation, not permission to use inexact floating-point threshold tests.
Integer comparisons of adjacent candidates provide an exact implementation.
For m=1, the j>=1 range is empty. For m=0, only P=0 exists.

The equality threshold between adjacent interior candidates j and j+1 is

    K 2^(2j+1)=h 2^m.                                          (3a)

Below this left-hand/right-hand ratio the j+1 candidate wins; above it j
wins. Equivalently the minimizing integer j is nearest to j*, with a tie
exactly when j* is a half-integer, subject to the indicated clipping. The
candidate j=0 must still be compared separately because it carries -S.

For fixed K,h>0 and S>=0, and every m with 1<=j*<=m-1, let

    d_m=distance(j*, integers).

Then the exact unrestricted maximum simplifies further to

    K 2^m - 2 sqrt(Kh) 2^(m/2) cosh((ln 2)d_m) + 2h.            (3b)

Indeed an interior candidate beats the j=0 candidate even before its
nonnegative penalty, and the j=m-1 candidate is positive. The deficit's
normalized coefficient is between 2 sqrt(Kh) and (3/sqrt(2))sqrt(Kh).
Here d_(m+2)=d_m: at fixed K/h the coefficient has period two in m, possibly
constant. Therefore replacing it by 2 sqrt(Kh) with O(1) error is generally
false. For K=h it alternates between 2K and 3K/sqrt(2), giving a hostile
control with an unbounded square-root difference on odd m. This is an exact
extremal asymptotic, not an approximation to typical occupied values.

When K>=h and 0<=S<=h, the minimum on the domain including zero is exactly
zero. For P>=1, the inherited elementary bound f(P)<=P-1 gives

    g(P)>=(K-h)P+h-S>=0.                                        (4)

The minimum statement is not asserted on a domain excluding zero.

### Low-bit cylinders

For t>=1 fixed low bits l and P=2^t u+l, formula (1) gives exactly

    g(P)=K'u-h f(u)-S'(u mod 2)+C',
    K'=2^t K-h rev_t(l),
    S'=h first_t(l),
    C'=Kl-h f(l)-S(l mod 2).                                   (5)

For t=0 use the unchanged parameters K'=K,S'=S,C'=0.
If K>h, then K'>0, so (3) applies to every bounded remaining u-word. This
replaces the inverse recursion's loose upper envelope by an exact one. The
remaining word length must be bounded explicitly; no infinite-cylinder
claim is intended. The scalar extrema do not specify which interior values
are occupied.

## 3. Exact maximum of each physical three-core rectangle

For a core state (h,q,A,G,T,S,E), the inherited count is

    F(P,R,b)=(qP+A)(2^(b+1)-rev_b(R))
       +h(2R+1-f(R))-PG-hf(P)-T-S(P mod 2)-E first_b(R).          (6)

Fix a,b. If a>0, impose 2^(a-1)<=P<2^a; if a=0, impose P=0 and a nonempty
core cut. In both cases take 0<=R<2^b. These are exactly the inherited legal
strong recipes in this fixed-length rectangle.

The maximum of (6) is attained among

    P in {2^a-2^j:0<=j<a}  if a>0, otherwise {0};
    R in {0} union {(2^ell-1)2^j:ell>=1,j>=0,ell+j<=b}.          (7)

Thus an exponential word search becomes at most
max(1,a) [1+b(b+1)/2] evaluations for each core state.

Proof. At fixed R the P objective is that of Section 2, with
K=q(2^(b+1)-rev_b(R))-G>0. Here q>=2h and G<=q-h, so K>=h.
The strict prefix-rounding increase therefore forces a one-run maximizer.
At fixed P set A0=qP+A>0. The R-dependent objective is

    2hR-hf(R)-A0 rev_b(R)-E first_b(R).

Every internal rounding from Section 1 strictly increases its first term,
does not decrease either negative term, and preserves the first-bit term.
It strictly increases the whole objective. Hence its nonzero maximizers
have only one run, possibly preceded by leading zeroes. This proves (7)
simultaneously: in a global maximizing pair neither word can have an internal
zero followed by a selected bit. No claim is made that a maximizer is unique.
The argument holds for any fixed physical core with the inherited coordinate
inequalities and strong recipe rule, not only this bank of 80 states.

## 4. Exact minima on arbitrary fixed-bit cylinders by a cut

Express (6) as a polynomial in the individual prefix and suffix bits x_i.
It has the form

    F=C+sum_i l_i x_i - sum_{i<j} c_ij x_i x_j,  c_ij>=0.        (8)

Indeed f supplies nonnegative pair coefficients, multiplied by -h; the only
cross-word quadratic term is -qP rev_b(R), also with nonpositive
coefficients. All other terms are linear or constant. Fixing any bits leaves
this same form. Thus even arbitrary cylinders, not only low-bit cylinders,
have polynomial-time exact minima.

For a self-contained reduction put d_i=2l_i-sum_{j!=i}c_ij. Then

    2F=2C+sum_i d_i x_i+sum_{i<j}c_ij |x_i-x_j|.                (9)

Make one node for each free bit; x_i=1 means the source side of the cut.
Join each pair by an undirected edge of capacity c_ij. Add source-to-i
capacity max(-d_i,0) and i-to-sink capacity max(d_i,0). Consequently

    min 2F = 2C+sum_i min(d_i,0)+minimum_cut_capacity.            (10)

The graph has O(a+b) vertices and O((a+b)^2) edges, with integer capacities
of polynomial bit length in the parameters. Standard augmenting paths give
a polynomial exact algorithm; the portable checker implements an elementary
breadth-first augmenting-path algorithm. A positive prescribed prefix length
is handled by fixing its leading bit to 1 before forming (8). This is an
extremum computation on existing recipes, not a new tournament construction.

## 5. Join-bank algebraic specialization

The join notation is different: let

    g(0)=1,  g(2u+1)=g(u),  g(2u)=g(u)+u+(u mod 2) (u>=1),
    e(v)=1-(v mod 2).

Then the exact identity

    g(v)+f(v)=v+e(v)

holds for all v>=0. It holds at zero; for either digit epsilon, the sum at
2u+epsilon is g(u)+f(u)+u+(u mod 2)=2u+1, because
e(u)+(u mod 2)=1. This equals 2u+epsilon+e(2u+epsilon), completing induction.
Consequently, defining F2(a,b,c;v)=a g(v)-b e(v)+c(1+v),

    F2(a,b,c;v)=(a+c)v-a f(v)-(a-b)(v mod 2)+(a-b+c).            (11)

For a>0,c>0,0<=b<=a, apply Section 2 using K=a+c,h=a,S=a-b,
C=a-b+c. This gives the exact maximum over 0<=v<2^m from (3), plus C,
and the minimum C at v=0. On positive m-bit v use (3) without zero, plus C.
The exact square-root-deficit formula (3b) and tie criterion (3a) transfer
unchanged. The lowercase a,b,c here are unrelated to the padding lengths.

For the physical interpretation, the inherited one-arc join identity is
THM-4875(D), file
`01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md`.
Its parameters are a=H(Y), b=the rooted endpoint Hamilton-path count, and
c=the rooted two-path-cover coordinate, so a,c>=1 and 0<=b<=a. Its v domain
is all nonnegative integers. The same notation is restated in the Setting
section of `05-knowledge/hypotheses/HYP-9478-one-apex-join-mass-laws-and-coverage-to-2-35.md`;
no unproved all-target hypothesis in that file is used.

## 6. Why sharp envelopes still do not close the conjecture

The exact prefix and suffix maximizers use at most two one-runs in total.
Session217's low-run theorem shows that recipes with any fixed run budget
have only polylogarithmically many targets up to X. Therefore solving these
extrema cannot be mistaken for an all-target construction. Exact minima and
maxima describe an interval hull, not its occupied odd points.

A literal hostile control: core index 0=(5,32,0,0,0,0,0), a=4,b=0, has
support {517,561,635,669,763,797,871,895}. Its sharp interval hull contains
519, which is absent from this rectangle. This is a rectangle-relative gap,
not a claim that 519 is absent from the full grammar.

The strongest achieved advance is a polynomial exact hull description plus
sharp inverse-cylinder maxima. No target-selection theorem, finite
inductive cone, carry automaton guaranteeing occupancy, or proof of all odd
N>=65 was obtained. The needed next step is an occupancy invariant that
survives inside these hulls; repeatedly sharpening endpoints alone cannot
provide it.

## 7. Reproduction

Run `python3 verify.py` and `python3 -O verify.py` from this directory. Normal
mode is load-bearing: its assertions perform the tests. Optimized mode is
only a smoke test and does NOT independently repeat assertions.

The retained `VERIFY.json` records 36,879 internal-rounding checks on all
words below 2^13, 3,573 exact maximum rectangles for all 80 states and all
splits a+b<=8, and 2,219 exact min-cut rectangles with a+b<=6, plus 2,219
fixed-interior-bit subcylinders. It also checks 324 arbitrary-parameter scalar
maximum problems and the join bridge on 8,192 inputs. A 400-bit two-free-bit
cylinder checks cut capacities exceeding 10^100; augmentations use the true
path bottleneck without a fixed numeric cap. Every count in
the rectangle maximum check is also compared to the independently expressed
returning coordinate compiler. The independent pair-sum evaluator checks f
against its recurrence for all words below 2^13. No large recipe file, finite
coverage extension, or unrestricted spectrum assertion is involved.

The targeted novelty check searched current canon and Session21x result
paths for single-run maxima, prefix extrema, and binary rounding; it found
no matching result. That bounded local search is not a literature claim.
