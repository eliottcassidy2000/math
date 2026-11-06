# Actual source slopes approximate a bit-reversal grid

Status: **RESERVED / UNREVIEWED PROOF CANDIDATE**. This is a constructive,
quantitative approximation theorem for actual source indices, with exact
rational controls. It does not prove exact target occupancy.

## 0. Source conventions and scope

Begin at the actual source (g,u,E)=(1,2,0). Its binary operations are

    O: (g,u,E) -> (g,2u,0),
    Z: (g,u,E) -> (g+u-E,2u-1,1).

Here E=u mod 2. These are the reviewed unit-source operations, and they
produce precisely the scalar recurrence

    g(1)=1, g(2h)=g(h), g(2h-1)=g(h)+h-(h mod 2).

The run-error derivation below gives a positive quantitative construction.
The optional joint discrepancy estimate includes an explicit rounding
correction. This packet is separate from
`../session201-digital-graph-obstruction/`.

The main object is an actual source, represented both by its chronological
O/Z word and by its ordinary index. The preserved information is the exact
index map; the discarded information in an approximate slope is a bounded
real error, which still matters for exact integer occupancy.

## 1. Every index in a dyadic interval has an explicit grid label

Let a word have length k>=0, with chronological bits epsilon_i in {0,1},
i=0,...,k-1, where 1 means Z and 0 means O. Put

    B=2^k,
    A=sum_i epsilon_i 2^i,
    R=rev_k(A)=sum_i epsilon_i 2^(k-1-i).

The fixed-length reversal includes leading zero bits; rev_0(0)=0. Then

    u=2B-R,                  t=(A+1/2)/B.                 (1)

Indeed u starts at2 and changes by u->2u-epsilon. Also t starts at1/2 and
changes by t->(t+epsilon)/2. These recurrences give (1) by induction.
Since bit reversal is a permutation of0,...,B-1, this is a bijection between
word labels A and all actual indices

    B+1 <= u <= 2B.

In particular all original indices remain >=2. At k=0 there is one empty
word: A=R=0, B=1, u=2 and t=1/2.

## 2. Exact run-error identity and one-sided approximation

Set e=g(u)-u t. Initially e=0. Operation O leaves e unchanged. For Z,

    e_new-e_old=(t_old+1)/2-E_old.                         (2)

This follows by substituting the two source updates into g_new-u_new t_new.

Consider a maximal consecutive Z-run of length r>=1. Its initial E is0.
Its initial t=s satisfies 0<s<=1/2: the run either starts at the seed, or
immediately follows O, which halves a number less than1. Inside the run,

    t_j=1-(1-s)/2^j.

The first error increment is(s+1)/2. Each subsequent increment is
-(1-s)/2^j, for j=2,...,r. Thus the WHOLE run increases e by exactly

    s+(1-s)/2^r.                                         (3)

This is positive and at most(1+s)/2<=3/4. A final Z-run ending at the end of
the word obeys the same formula; no terminal O is required. Summing over
all maximal Z-runs therefore gives

    e=sum_runs [s_run+(1-s_run)/2^(length_run)],
    0 <= e <= 3 Zruns/4
            <= Delta_k := 3 floor((k+1)/2)/4
            <= 3(k+1)/8.                               (4)

The run count bound is floor((k+1)/2), since distinct runs require separating
O letters. The lower bound is an all-word statement: individual later Z
steps may decrease e, but each complete or truncated run contributes a
positive amount. For k=0, Zruns=Delta_0=e=0.

Boundary checks are exact. A=0 gives the all-O word, u=2B, g(u)=1, e=0.
For k>=1, A=B-1 gives the all-Z word, u=B+1, g(u)=u, and

    e=(B+1)/(2B)=1/2+1/(2B),

as (3) with s=1/2 and r=k also gives.

## 3. Slope discrepancy and an explicit approximate constructor

Define the physical slope s_A=g(u_A)/u_A, with u_A=2B-rev_k(A), and let

    epsilon_k = Delta_k/(B+1),
    epsilon'_k = 3(k+1)/(8(B+1)).

By (1) and (4),

    t_A=(A+1/2)/B <= s_A <= t_A+epsilon_k
                                 <= t_A+epsilon'_k.      (5)

All slopes lie in[0,1], by1<=g(u)<=u. Some equal1; these endpoints are not
removed from the following definition.

For a finite list of B points in[0,1], use the half-open star discrepancy

    D*=sup_(0<=z<=1) | #{A:s_A<z}/B-z |.

The midpoint grid t_A has star discrepancy1/(2B). If F_s and F_t are their
strict-cutoff distribution functions, (5) gives

    F_t(z-epsilon_k) <= F_s(z) <= F_t(z),                 (6)

where F_t(w)=0 for w<0. The midpoint-grid bound and (6) imply

    D*(s_A) <= 1/(2B)+epsilon_k
            <= 1/(2B)+3(k+1)/(8(B+1)).                  (7)

More precisely, F_s(z)-z<=1/(2B), while z-F_s(z)<=1/(2B)+epsilon_k.
When z<epsilon_k, the latter follows directly from F_s>=0; otherwise apply
the midpoint-grid bound at z-epsilon_k. Strict cutoffs and slopes equal1
cause no exception to the argument.

For ANY requested real slope theta in[0,1], choose

    A=min(B-1,floor(B theta)),
    u=2B-rev_k(A).                                       (8)

Then |t_A-theta|<=1/(2B), so the actual index u satisfies

    theta-1/(2B) <= g(u)/u
       <= theta+1/(2B)+epsilon_k,
    |g(u)/u-theta| <= 1/(2B)+epsilon_k.                  (9)

This is an explicit actual-index constructor, not a search over free profile
coordinates. It covers theta=0 and theta=1, and also k=0 with u=2. In
particular the dyadic slope lists become uniformly distributed with the
explicit error in (7), and any prescribed slope can be approached by choosing
k larger. No exact prescribed rational slope is asserted.

## 4. The auxiliary bit-reversal grid: careful rectangle rounding

Consider the classical base-two bit-reversal grid

    H_k={(A/B,rev_k(A)/B):0<=A<B}.

This is the two-dimensional Hammersley point set in base two. The name is
standard; all estimates used here are proved below, not imported from a
sharp discrepancy theorem. A primary-source definition is recorded in
`SOURCES.md`.

### 4.1 Dyadic blocks reverse into arithmetic progressions

Take an aligned index block

    I={q 2^ell,...,(q+1)2^ell-1},  0<=ell<=k.

Writing A=q2^ell+j shows

    rev_k(A)=rev_(k-ell)(q)+2^(k-ell) rev_ell(j).

Hence the reversed values are, as a set, an arithmetic progression of
M=2^ell terms, with step B/M and offset r0=rev_(k-ell)(q), where
0<=r0<B/M. For any real cutoff y in[0,1], put eta=r0/(B/M) in[0,1).
The number C of points below y is

    C=ceil(M y-eta).

This ceiling already lies in[0,M]. Consequently

    -eta <= C-M y < 1-eta,
    |C-M y|<1.                                          (10)

This is the exact rounding estimate behind the discrepancy bound.

### 4.2 Grid-aligned first cutoff

Suppose the first cutoff is x=a/B with integer0<=a<=B. For0<a<B,
[0,a) splits into popcount(a)<=k aligned dyadic integer blocks. For a=B,
one full block suffices. For a=0 the count is zero. Summing (10) gives, for
k>=1 and arbitrary real y in[0,1],

    | #{A:A/B<x,rev_k(A)/B<y}/B-x y | <= k/B.              (11)

Thus only the FIRST cutoff needs to be grid-aligned for this estimate.
Both grid-aligned cutoffs are included as a special case.

### 4.3 Arbitrary real cutoffs

For arbitrary x in[0,1], let a=ceil(Bx). The first-coordinate selection is
exactly A<a, and0<=a/B-x<1/B. Applying (11) at a/B gives, for k>=1,

    | #{A:A/B<x,rev_k(A)/B<y}/B-x y | <= (k+1)/B.          (12)

The extra1/B pays for rounding the first cutoff. At k=0 there is one point,
(0,0), whose real-rectangle star discrepancy is1. This boundary is handled
separately, and also satisfies the numerical bound(k+1)/B=1.

**Correction to an informal stronger claim.** The unqualified bound k/B for
ALL REAL rectangles is false. At k=1 the grid is{(0,0),(1/2,1/2)}. Taking
both cutoffs just above1/2 gives discrepancy tending to3/4, greater than
k/B=1/2. The proof only supplied k/B before the first-cutoff rounding step.
This correction changes neither (1)--(9) nor the grid-aligned estimate (11).

No optimal joint discrepancy constant is claimed. Formula (12) is the
explicit elementary bound used here.

## 5. What this gives, and does not give, for the J form

Let u and v come from word lengths k and ell, with midpoint labels t_u,t_v
and actual errors e_u,e_v. Put p=u mod2 and q=v mod2. The exact identity is

    Psi_J(u,v)+pq-uv(1+t_u t_v)
       =u t_u e_v+v t_v e_u+e_u e_v.                     (13)

Since0<t_u,t_v<1 and0<=e_u<=Delta_k,0<=e_v<=Delta_ell,

    0 <= Psi_J(u,v)+pq-uv(1+t_u t_v)
      <= u Delta_ell+v Delta_k+Delta_k Delta_ell.          (14)

Thus the grid approximation transfers to actual bilinear values with a
quantified error. That error generally grows with the source sizes; it is
not a bound below one. Uniform slope approximation or rectangle occupancy
therefore does not imply that a prescribed odd integer is attained exactly.
The missing integer-value control is retained explicitly in (13).

The positive geometric description is compatible with the separate proof
that the exact binary graph of g is not regular. Approximate distribution,
exact graph recognition, and exact target support are different predicates.

## 6. Reproduction and limits

Run the standalone verifier:

    python3 session201-hammersley-source-approximation/verify.py
    python3 -O session201-hammersley-source-approximation/verify.py

It uses only Python's standard library and exact integer/Fraction arithmetic.
It checks every word through k=12, the run identities, all boundary cases,
the slope discrepancy and approximate constructor, dyadic-block arithmetic
progressions and cutoff errors, exact small joint discrepancies including
the k=1 hostile, and (13)--(14). The reports specify all finite control sizes.
These are controls of the all-order proofs, not a target-coverage census.

No new exact occupancy statement is claimed. The twelve-template all-target
lemma remains OPEN.
