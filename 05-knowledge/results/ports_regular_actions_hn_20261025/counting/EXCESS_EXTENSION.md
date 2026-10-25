# Exact excess grading repairs the gauge, and why an ungraded scalar cannot

2026-10-25. **PROVED:** a bivariate extension of the ordinary defect gauge
to the inherited pair-free-core substitution; its exact individual-future
equivalence; and a scalar-normalization obstruction without that grading.
**FINITE-EXACT:** small independent substitution and literal graph controls.
No reduction in total core-engine states, runtime, or bit complexity is
claimed. The [reviewed ordinary recurrence](README.md) remains the main
counting result.

## 1. The two degrees have different jobs

The inherited identity is

    Z_C(p)=Z_T(p_d/(1+p_d)),
    (1+z^d)^(-k)=sum_(e>=0) beta(k,e) z^(d*e),
    beta(k,e)=(-1)^e binom(k+e-1,e),  k>0.

For k=0 define beta(0,e)=1 if e=0 and 0 otherwise. The source and its
normalization are in [the signed-core transfer, section 4](../../tournament_enumeration_20261021a/README.md)
and [the signed ordered-list transform](../../defect_boundary_20261020/boundary.md).

Introduce the formal bivariate series

    B(s,z)=Z_T(s^d/(1+z^d)),
    b_(M,E)=[s^M z^E] B(s,z).

Here M is the actual cycle degree of the source Burnside term and E is
the substitution excess. Both are auxiliary formal gradings. They are not
two vertex counts of a particular resulting core. On setting s=z=x,

    number of pair-free tournament classes on N vertices
       =sum_(M+E=N) b_(M,E).                          (X1)

The individual b_(M,E) are signed integers: (-1)^E*b_(M,E) counts
isomorphism classes of M-vertex tournaments whose vertices carry
nonnegative integer colors with total color sum E. An invariant coloring
is constant on each permutation cycle. Only odd cycles can fix a
tournament, and a length-d cycle contributes
sum_(j>=0)(-z)^(d*j)=1/(1+z^d). Weighted Burnside therefore gives exactly
B(s,z). These are colored source tournaments, not resulting cores with
two size statistics. The signed diagonal sum is retained exactly. This
corrects the draft claim that the b_(M,E) could be nonintegral; see the
[independent review](INDEPENDENT_REVIEW.md). This refinement does not
recover literal clone owners or ports from the aggregate coefficients.

## 2. Fixed actual degree removes cycle count again

Fix one coefficient b_(M,E). Process odd lengths in decreasing order.
For a prefix, let u be actual degree already used, v excess already used,
m its actual cycle count, and

    R=M-u, Q=E-v,
    f(d)=sum_(past actual c) gcd(c,d),
    D(d)=d*m-f(d),  D(1)=0.

Replace its prefix weight w by w_tilde=2^(R*m) w. Every decorated future
has actual cycle degree exactly R, independently of its excess Q. Thus

    sum_(future actual d) f(d)=R*m-sum_(future actual d) D(d).   (X2)

The same gauge therefore removes m. The new state is (R,Q,D) at the
current length. It is not (R+Q,D): the two degrees cannot be silently
combined. The excess choices made at each processed length are completed
when that length is processed, so no unresolved past allocation is hidden
in this state.

More explicitly let G_ell(R,Q,D) sum all future odd partitions lambda of
actual degree R with parts at most ell, together with excess assignments
e_d>=0 satisfying sum d*e_d=Q. A positive e_d is allowed only if the
actual multiplicity k_d is positive. Its summand is

    2^[e(lambda)-sum_(d in lambda)D(d)] / z_lambda
       * product_d beta(k_d,e_d),                    (X3)

using the same inherited exponent e(lambda) and centralizer z_lambda as
the ordinary recurrence. There is no prefix cycle count in this formula.

## 3. The exact bivariate recurrence

At odd length ell select k actual cycles and e excess units of that length.
Put R'=R-ell*k, Q'=Q-ell*e, and on affordable smaller odd d put

    D'(d)=D(d)+k*(d-gcd(ell,d)).

Then

    G_ell(R,Q,D)
      =sum_(0<=k<=floor(R/ell)) sum_(0<=e<=floor(Q/ell))
         beta(k,e) / (ell^k*k!)
         * 2^[k*(R-D(ell))-(ell*k*k+k)/2]
         * G_(ell-2)(R',Q',D').                      (X4)

Only e=0 contributes when k=0. If ell>R, skip that length without reading
an unavailable D(ell): no actual cycle of that length can support excess.
The exponent and D update are exactly those of the ordinary gauge; excess
has changed neither actual cycle count nor gcd interactions. The proof is
the ordinary batch-transition calculation using (X2), multiplied by the
known beta(k,e). This proves (X4) term by term, before signed cancellation.

Boundary values are

    G_ell(0,Q,D)=1 if Q=0, otherwise0;
    G_ell(R,Q,D)=0 for R>0 and ell<1;
    G_1(R,Q,D)=beta(R,Q)*2^binom(R,2)/R!,  R>0.      (X5)

The last expression follows because every remaining actual cycle and every
excess unit must have length one. For any sufficiently large odd L,

    b_(M,E)=G_L(M,E,0).                              (X6)

Equations (X1), (X4)--(X6) are an exact unrestricted pair-free coefficient
recurrence, with cycle count genuinely absent from each graded state.
They trade it for a second degree coordinate. This is an exact repair,
not evidence that the trade reduces the total number of states.

## 4. Exactly which defect coordinates remain observable

Fix R>=1, Q>=0 and a descending stage ell>=1. A coordinate for odd d
is observable by an individual decorated continuation exactly when

    d<=min(ell,R), and either d<R or R divides Q.     (X7)

To see necessity at d=R, its one actual cycle exhausts R and is the only
available length for excess; therefore Q must be a multiple of R. For
d<R, the actual partition consisting of one d-cycle and R-d singletons
exists, and its singleton block can carry all Q excess units. For d=1
use the all-singleton partition directly. Every such beta coefficient
is nonzero.

After removing the future's known intrinsic weight and beta coefficient,
the all-singleton continuation has normalized interaction1, since D(1)=0.
The stated probes have interaction2^(-D(d)). Thus two histories have
identical normalized multipliers for every individual decorated future
iff their D values agree on exactly the set (X7). The probes read each
coordinate in that set, apart from the permanently zero D(1); coordinates
outside it never occur. This identifies the exact equivalence relation
for the individual-future observable. It does not assert independent
coordinate freedom on a restricted set of reachable prefix profiles,
independent storage necessity, or a minimal bit encoding. It also makes
no claim about methods preserving only the final signed sum.

There is a small extra horizon refinement: if R is odd and Q is not a
multiple of R, the affordable coordinate D(R) is nevertheless invisible.
For example

    G_3(3,2,D)=8,
    G_3(3,3,D)=-40/3-(2/3)*2^(-D(3)).               (X8)

The first future is forced to have three actual singletons; the second
can also have one actual 3-cycle with one excess unit of length three.
This is a support condition, not a cancellation argument. The reference
checker keeps the larger sufficient coordinate mask for simplicity and
tests (X7) separately.

## 5. Ungraded individual futures recover the whole old response

There is also a sharp negative result. Fix only the total remaining
substitution degree S, and suppose length one remains available.
Let V_f be the vector of the nonzero individual decorated-future
multipliers as a function of the old prefix response f. Arbitrary scalar
prefix normalization can merge f and g only if V_f and V_g are
proportional.

**Theorem.** For S>=2, these vectors are proportional iff

    f(d)=g(d) for every odd d<=min(ell,S).            (X9)

**Proof.** Consider first the future with one actual singleton and S-1
singleton excess units, and then the one with two actual singletons and
S-2 excess units. Both have nonzero intrinsic signed coefficients. After
dividing those known coefficients, the ratio of the second response to
the first is2^f(1)=2^m. Proportional vectors therefore have the same m.

For each odd 1<d<S, use one actual d-cycle and one actual singleton, with
all S-d-1 excess units on the singleton. Dividing its normalized response
by the one-singleton probe gives2^f(d). This determines f(d). If d=S is
an available odd length, use its one cycle with no excess; its response
relative to the one-singleton probe is2^(f(S)-m), and m is already known.
The converse follows directly from equality of every response. QED.

Consequently, for S>=2 the old affordable f profile is already the
coarsest projective state for preserving every individual ungraded core
continuation. A scalar gauge cannot create the ordinary theorem's strict
mergers on this ungraded observable. At S=1 there is only one future and
its response can be absorbed into the scalar; at S=0 there is none to
distinguish. Omitting length one requires a different support analysis
and is outside this theorem.

For a concrete hostile use histories (27) and (9,9,9), which have the
same D(1),D(3) but cycle counts1 and3. The ratios of their interactions
on one actual singleton and on two actual singletons are respectively
1/4 and1/16. No single nonzero rescaling aligns both. Fixing actual
degree removes this second probe from a single graded fibre, which is
exactly why (X2) succeeds. A polynomial weight that records every actual
degree is another expression of this grading; it does not make the
information disappear for free.

## 6. Controls and the next useful target

Run `python check_excess.py` and `python -O check_excess.py` from this
directory. Both passed with explicit exceptions rather than assertions:

- All120 bivariate coefficients with M+E<=14 agree with independent
  unmerged odd partitions and direct products of the truncated factors
  (1+z^d)^(-k).
- Summing the grading agrees with literal enumeration and isomorphism
  canonicalization of all tournaments with no homogeneous pair for
  orders0..5. The counts are1,1,0,1,2,5. This also keeps the meaning of
  pair-free distinct from prime or indecomposable.
- All288 relevant pairs of small formal defect profiles in the declared
  R=1..8, Q=0..8, ell in{1,3,5} universe agree on individual supported
  future interactions iff they agree on (X7).
- The ungraded scalar hostile and both visibility cases in (X8) pass;
  invalid degree inputs are rejected under normal and optimized Python.

The precise next target is an excess-graded implementation of the
existing odd-factorial transfer, sharing memoized (R,Q,D) states across
the diagonal sum (X1), followed by a state and arithmetic-cost comparison.
Until that comparison, only the mathematical quotient within each graded
fibre is established. The obstruction (X9) explains why simply applying
the ordinary scalar normalization to the old ungraded core engine cannot
produce its strict response mergers.
