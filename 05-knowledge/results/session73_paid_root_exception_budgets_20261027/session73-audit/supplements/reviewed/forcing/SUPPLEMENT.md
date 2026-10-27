# Scalar payment below 1/42; the spatial observation still pays inverse masses

2026-10-27. Additive Session73 continuation; the independently audited core
`REPORT.md` and `CORE_MANIFEST.json` are unchanged. **PROVED** carrier identities,
**FINITE-EXACT** rational scalar certificate, and **OPEN** spatial-current bound.
Dependencies are the source snapshots already pinned in `SOURCE_PINS.json`.

## S1. A global scalar budget from a finite authenticated cutoff

The core proves, for the native shortcut-phase-aware passive forcing,

```
E = sum_(n>1 odd) W(n)/(6n+1) <= 13/57.
```

Keep the **same** native weights. For an odd integer cutoff N, let

```
M_N=sum_(3<=n<=N,odd) W(n),
E_N=sum_(3<=n<=N,odd) W(n)/(6n+1).
```

When those finitely many weights have authenticated strict ROOT receipts,
nonnegativity and the inherited `sum_(n>1)W(n)<=13/3` give

```
E_N <= E <= E_N + (13/3-M_N)/(6(N+2)+1).                 (S1)
```

There is no convergence assumption about larger sources. If an omitted source
is unrooted, its defined W is zero and the estimate remains valid. The simple
weight `1/(6n+1)` supplies the tail decay absent from unweighted atom positivity.

At N=5001 the checker authenticates all 2,500 nonroot odd sources independently
of any conjecture. The largest observed odd-step count is 87. Exact rational
arithmetic and outward rounding give

```
47106529/2000000000 <= E <= 236498201/10000000000 < 1/42.
```

Thus the original native source probability pays the entire first-hit
phase-aware forcing with **E<1/42**. The lower decimal is approximately
0.0235532645 and the upper is 0.0236498201; the proof uses the displayed
fractions, not decimal rounding.

The producer checker uses cached actual counters and the integer-binomial
identity `W=2/[(tau+K+1) binom(tau+K,K)]`. The parent's independent replay uses
source-by-source literal words and factorials; both reproduce the same interval.
The parent owns that independent script and its audit provenance. Neither the
cutoff nor the finite verification of its 2,500 sources is a universal ROOT proof.

The same principle makes E computable to arbitrary rational precision even
without deciding every source's support: finite source values have inherited
uniform pointwise approximations, while the cost `1/(6n+1)` controls the tail.
Computability of this one aggregate still does not decide whether a designated
W(n) is positive.

## S2. A marked signed occupation current

Let a source s be sampled once with its original `lambda(s)`. Before its first
ROOT hit, put `v_t=U^t(s)` and retain the actual accumulated shortcut parity

```
epsilon_(s,t)=(-1)^(a_1+...+a_t).
```

Its literal integer clock is separately `t+a_1+...+a_t`. The native equation
`f(U n)=(-1)^a(n) f(n)` implies

```
epsilon_(s,t) f(s)=f(v_t).                               (S2)
```

Define the marked target occupation and signed transported current

```
mu_t(v)=sum_(s:v_t=v,t<tau(s)) lambda(s),
J_t(v)=sum_(s:v_t=v,t<tau(s)) lambda(s)epsilon_(s,t)f(s).
```

Every term retains the same source, shared target port and clock. Equation (S2)
gives `J_t(v)=f(v)mu_t(v)`. The inherited occupation identity therefore yields

```
J(v)=sum_t J_t(v)=W(v)f(v),
sum_(v>1:W(v)>0) |J(v)|^2/W(v)=sum_(v>1)W(v)<=13/3.       (S3)
```

All sums are absolutely justified by finite total occupation. For a finite
collection of authenticated lambda-sources, the same construction gives
`J_F(v)=f(v)mu_F(v)`, `0<=mu_F(v)<=W(v)`, and
`sum |J_F(v)|^2/W(v)<=sum mu_F(v)<=13/3`.

This is a legitimate source-preserving phase-transported weighted current. It
is **not** the raw sum `sum lambda(s)f(s)` over arrivals. The exact finite test
already finds a discrepancy at shared target 5. In particular, source 3 has
sign -1 while its target 5 has sign +1: the intervening native odd clock must
be applied. Dropping the phase is not an innocuous simplification.

## S3. Why this still does not bound the original spatial ROOT current

The actual inherited target is

```
R_f(N)=sum_(N<n<=4N) (3[n=2 mod3]-1) f(n)
      =sum_(m odd) K_N(m) f(m),
K_N(m)=w(m)[N>=ceil(m/4)]
       +(-3 chi(m)-w(m))[N>=ceil(m/2)],
```

with `chi=0,1,-1` on residue classes `0,1,2 mod3`, and `w=3[.=2 mod3]-1`.
This is a spatial counting-weighted signed observer, not a lambda trajectory
observer. The root contribution is always `K_N(1)f(1)=-3`.

On a finite certified nonroot support V, (S3) does give the exact readout

```
R_f(N)+3 = sum_(v in V) [K_N(v)/W(v)] J(v)                (S4)
```

provided V contains the nonzero terms of that row. It costs the inverse
source masses. Precisely, the linear row observer on the Hilbert space with
norm `sum |J(v)|^2/W(v)` has squared operator norm

```
C_N^2 = sum_(v in V) K_N(v)^2/W(v).                      (S5)
```

This follows by Cauchy-Schwarz and is attained by a scalar multiple of
`J(v)=K_N(v)`; it is a statement about the whole declared linear input space,
not an assertion that the native current attains this worst case.

For a window `X<=N<2X`, the squared operator norm is at least
`max_v sum_N K_N(v)^2/W(v)` and at most the sum of those column quantities.
The lower bound tests one retained source; the upper is the Hilbert-Schmidt
bound. The exact controls give these lower bounds versus the actual centered
current energy:

| X | Actual `sum(R_f+3)^2` | Lower bound on the linear observation norm squared |
|---:|---:|---:|
| 4 | 124 | 11,274,451,650 |
| 16 | 592 | 3,066,650,848,800 |
| 64 | 4,808 | 906,443,258,043,456 |
| 128 | 16,832 | 2,116,758,473,188,128,000 |

These are finite exact operator lower bounds, not an all-X growth theorem.
They demonstrate that the small temporal payment and the source-weighted
current norm cannot simply be called a useful bound for the specified spatial
observer. Restoring the original predicate requires its actual sign-covariance
structure or a much more selective observation estimate.

At a source with unknown W=0 the positive-weight Hilbert observation (S4) is
not available: defining an inverse mass there would be circular. Mathematically
its native eventual f is also zero by definition, but proving that a prescribed
integer has positive f-support remains the original unsolved obligation.

A phase gauge likewise has to conjugate the spatial kernel. Replacing f(v) by
1 removes the native clocks from `f^T G_X f`; keeping `G_X` unchanged changes
the target. The genuine folded-current uniform estimate from the October17
packet is therefore **still OPEN**. Neither (S1) nor (S3) supplies it.

## S4. Reproduction and exact scope

```
python3 check_supplement.py > supplement.normal.json
python3 -O check_supplement.py > supplement.optimized.json
cmp supplement.normal.json supplement.optimized.json
```

There are 3,169 always-active exact gates. The universe comprises all 2,500 odd
nonroot sources through 5001; 84 lambda-source trajectories through source501,
with 1,694 marked nonroot visits and 350 merged target vertices; all direct
versus folded current rows N=1,...,256; and six exact observation windows
X=4,8,16,32,64,128. No censused exponent is extrapolated. The scalar source
counters and spatial response vector have separate SHA-256 receipts.

The full all-order proof remains in the core. This supplement closes the
explicit scalar certificate and a typed weighted-to-spatial interface; the
independent audit must verify these additions before their publication.
