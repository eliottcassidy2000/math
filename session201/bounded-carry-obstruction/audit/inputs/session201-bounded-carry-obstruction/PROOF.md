# Bounded source approximation error cannot cover a tail

Status: **RESERVED / UNREVIEWED PROOF CANDIDATE**. This is a separate,
all-order obstruction to bounding BOTH source approximation errors in a
finite-template construction. The twelve-template occupancy lemma remains
**OPEN**. No exact target counterexample is asserted.

## 1. The actual-source error and a lower bound per Z-run

A length-k binary source word starts at u=2,g=1,E=0,t=1/2. The updates are

    O: u'=2u,   g'=g,       E'=0, t'=t/2,
    Z: u'=2u-1, g'=g+u-E,   E'=1, t'=(t+1)/2.

For e=g-u t, O leaves e unchanged, while Z changes it by(t+1)/2-E.
Let the word have R maximal Z-runs. For run i, let o_i be the immediately
preceding O-run length, and r_i>=1 the Z-run length. The first o_i may be0;
all subsequent o_i are positive. Trailing O letters are not assigned to a
Z-run, so

    sum_i(o_i+r_i)<=k.                                   (1)

Just before the first O block, t=1/2. After any preceding Z-run, t>1/2.
Halving o_i times therefore shows that the start value t_i of run i obeys

    2^(-o_i-1)<=t_i<=1/2.                               (2)

The upper bound holds because a noninitial run follows at least one O,
and0<t<1 throughout; an initial run without O starts exactly at1/2.

Summing the Z-step error increments gives the exact run contribution

    eta_i=t_i+(1-t_i)2^(-r_i).                           (3)

For example, the first increment is(t_i+1)/2, while later increments are
-(1-t_i)2^(-j), j=2,...,r_i; the finite geometric sum gives (3).
Combining(2)--(3) with the arithmetic-geometric mean inequality yields

    eta_i>=2^(-o_i-1)+2^(-r_i-1)
          >=2^(-(o_i+r_i)/2).                           (4)

Thus every Z-run contributes positive error. If R=0 then e=0. If R>=1,
convexity of x->2^(-x/2), followed by(1), gives

    e=sum_i eta_i
      >=R*2^(-sum_i(o_i+r_i)/(2R))
      >=R*2^(-k/(2R)).                                  (5)

This is an exact all-word bound. Individual later Z steps can lower e, but
that does not invalidate the positive whole-run sum.

## 2. Fixed bounded error permits only sublinear run counts

Fix E>=0 independently of k. Consider words with e<=E.

If E=0, (4) forces R=0, so only the all-O word occurs at each length.
Assume E>0. For any fixed delta>0, R>=delta k would imply, by(5),

    E>=delta k*2^(-1/(2delta)).                          (6)

This is impossible whenever

    k>(E/delta)*2^(1/(2delta)).                         (7)

Consequently the maximum possible R among bounded-error words satisfies
R=o(k). The estimate is uniform over all such words at each length, not
merely an assertion about a chosen sequence of words.

A useful equivalent necessary inequality, when R>E, is

    2R log_2(R/E)<=k.

No rounding or numerical root approximation is needed for exact finite
checks: (5) is equivalent to

    e^(2R) * 2^k >= R^(2R),

which is a rational inequality for every actual source word.

## 3. The bounded-error source set is subpolynomial

The exact number of length-k binary words with R Z-runs is

    binom(k+1,2R).                                      (8)

To see this, select the2R distinct run boundaries from the k+1 boundaries
between or outside letters; consecutive selected pairs specify the starts
and ends of the Z-runs. Distinct boundaries between adjacent runs leave at
least one O letter.

If r_k=o(k) bounds all permitted run counts, (8) gives

    #bounded-error words of length k
      <=sum_(R<=r_k) binom(k+1,2R)=2^(o(k)).             (9)

For an elementary estimate, eventually2r_k<=(k+1)/2, so the sum is at most
(2r_k+1) times the largest binomial coefficient in that range. The bound
binom(n,j)<=(euler_number*n/j)^j then makes its base-two logarithm at most

    log_2(2r_k+1)+2r_k log_2(euler_number*(k+1)/(2r_k)),

which is o(k); the r_k=0 case has count1. Here euler_number is the constant
exp(1), not the source error e.

For each actual u>=2, there is exactly one binary source word, whose length
k is characterized by2^k<u<=2^(k+1). Its approximation label is

    B=2^k, A=rev_k(2B-u), t=(A+1/2)/B,
    e(u)=g(u)-u t.

Therefore summing(9) over k<=floor(log_2 N) proves

    C_E(N):=#{2<=u<=N : e(u)<=E}=N^(o(1)).               (10)

Precisely, for every eta>0 there is N0(E,eta) such that C_E(N)<=N^eta for
N>=N0. To derive this from(9), first bound every sufficiently large length-k
count by2^(eta*k/2), sum the resulting geometric series, and absorb its
constant and finitely many initial lengths into N^eta. For E=0,
C_0(N)=floor(log_2 N) for N>=2, since the indices are exactly powers of two.

## 4. Obstruction to a bounded-error finite-template constructor

The reviewed central-arithmetic product bound for each of the twelve forms
is

    Psi_W(u,v)>=uv/2,   u,v>=2.

This is the already audited source-cone inequality from the central
arithmetic proof, Section4; it is not inferred from the finite checks here.
A value <=N therefore requires u,v<=N. By(10), the number of distinct target
values <=N with a witness satisfying BOTH

    e(u)<=E and e(v)<=E

is at most

    12 C_E(N)^2 = N^(o(1)).                             (11)

It cannot contain a tail of the odd integers, whose counting function has
linear growth. More generally the same conclusion holds for any fixed
finite menu with an analogous positive product bound; two different fixed
error bounds may be replaced by their maximum.

This is a barrier to uniformly bounded correction error in BOTH source
coordinates. It does not exclude a method that bounds only one error, lets
the error allowance grow with N, or uses a different kind of correction.
It does not say that every target or every witness must have large errors.
Neither(10) nor(11) is a statement about the unrestricted target support.

## 5. An explicit unbounded-error family

As a hostile boundary check, take the actual word(ZO)^m. Before the j-th
pair, j=0,...,m-1,

    t_j=1/3+1/(6*4^j).

Each isolated Z adds(t_j+1)/2 to e, and O leaves e unchanged. Hence

    e((ZO)^m)=2m/3+(1-4^(-m))/9.                         (12)

This grows linearly with word length2m, and is compatible with the general
O(k) upper bound in the positive Hammersley approximation theorem. The
bounded-error obstruction is not based on an assumption that all errors
are bounded in the first place.

## 6. Portable exact controls and dependencies

    python3 session201-bounded-carry-obstruction/verify.py
    python3 -O session201-bounded-carry-obstruction/verify.py

The verifier is standalone, using only standard-library integer and rational
arithmetic. It checks the lower bounds for every word through length14,
including the rational-power form of Jensen's inequality. It verifies(8),
small bounded-error counts for E=0,1/2,1,2, and the exact family(12) through
m=80. These are word-level controls, not a target-coverage census.

The dyadic index bijection and upper approximation bounds are proved in the
separate `../session201-hammersley-source-approximation/PROOF.md`. The
product-bound dependency is the released
`../session201-compact-central-release/session201/compact-central-arithmetic/dependencies/CORE_PROOF.md`,
Section4. No code import or external runtime dependency is required.

Uniform twelve-template occupancy remains OPEN.
