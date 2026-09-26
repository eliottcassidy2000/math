# A contracting Bellman operator computes the logarithmic optimum

Status: **PROOF CANDIDATE; root independently accepted the main derivation.
Exact finite controls retained; final audit pending.**

The global two-step pairing constant `B_*` lies in the certified interval

```
750086585295922675/2369190669160808448 <= B_*
    <= 750158422907862637/2369190669160808448,
0.31660034587322217... <= B_* <= 0.3166306674564004... .
```

Its width is less than `0.000030322`. The interval is exact; no displayed
midpoint or decimal extrapolation is claimed to equal `B_*`.

## 1. Inheritance and scope

[THM-4491](../../01-canon/theorems/THM-4491-pairing-two-step-tree-extra-density.md)
identifies the legal pairing tree and minimum lower natural density alpha.
[THM-4492](../../01-canon/theorems/THM-4492-pairing-two-cutoff-density-separation.md)
shows alpha cannot be attained as a natural density.
[THM-4496](../../01-canon/theorems/THM-4496-pairing-kernel-convolution-and-log-density.md)
identifies `B_*` as the finite-kernel supremum, the limit of the harmonic
finite optimum divided by `log X`, and an attained minimum logarithmic
density. Its preceding numerical lower bound was `0.3067797897...`.

The closest mechanism is its harmonic subtree root-gap recursion. The
hostile inherited policy `epsilon_i=v3(i) mod2` satisfies selected private
constraints but is not globally legal: `6 -> 9 -> 14`. The corrected near
miss is replacing joint constraints by marginal prices. The least-used
sidecar is the signed difference of the two fixed-root subtree costs.

Board: harmonic root differences; Haar measure on binary residues;
three-way ancestral resets; legal finite policies; primal/dual residual;
natural-density infimum versus attainment. Anchor: a quantitative upper
bound complementing the kernel lower certificates. Niche: an exact finite
policy polytope. Wildcard: a contraction that holds in mean but fails in
the supremum norm. The contraction overtakes another coefficient search.

All claims concern the modified pairing model. A pair bit gives

```
0: 2i-1 -> 3i-1,  2i -> i;
1: 2i-1 -> i-1,   2i -> 3i.
```

Global two-step descent is equivalent to

```
epsilon_(2k)+epsilon_(3k)=1,
epsilon_(3k+1)<=epsilon_(2k+1)<=epsilon_(3k+2).
```

The vertices are pair indices, and the edges are these exact legal
constraints. No tournament orientation or transfer to ordinary Collatz
is made.

## 2. Explicit legal policies have natural densities

Fix `M=2^s` and a bit `t_r` for each residue modulo `M`. Whenever a child
bit is free under its parent's clause, set it to `t_(child mod M)`.
Equivalently, setting the harmless bit `epsilon_1=0`, recurse for `n>=2`:

```
n=0 mod3: epsilon_n=1-epsilon_(2n/3),
n=1 mod3: epsilon_n=t_(n mod M) epsilon_((2n+1)/3),
n=2 mod3: epsilon_n=t_(n mod M)
                       +(1-t_(n mod M)) epsilon_((2n-1)/3).
```

Every parent index is smaller for `n>=2`. This defines one global legal
pairing. The same free-choice bit is used at both types of free child;
this coupled policy class should not be confused with two arbitrary
independent decisions.

Read an index's successive ancestor residues modulo three. Branch zero
flips the bit; branch one resets to zero if `t=0` and otherwise copies;
branch two resets to one if `t=1` and otherwise copies. At every state
exactly one of the three branches resets. For each fixed initial residue
modulo `M`, any specified length-`k` branch word defines exactly one
residue modulo `3^k`, by successive affine inversion; multiplication by
two is invertible modulo three. The Chinese remainder theorem supplies
one cylinder modulo `M*3^k`. Therefore unresolved mass after `k` levels
is exactly `(2/3)^k`, independent of the chosen policy.

Resolved cylinders fix the final bit, including the parity of earlier
flips. Their densities converge, and the unresolved mass tends to zero.
Thus conditional natural densities `x_r` and overall natural density
`delta=(1/M)sum_r x_r` exist. Put

```
P_0(r)=2r/3 mod M,
P_1(r)=(2r+1)/3 mod M,
P_2(r)=(2r-1)/3 mod M,
```

where division by three means multiplication by its inverse modulo `M`.
The same reset/CRT argument works after any further dyadic conditioning
`n=r mod 2^a`, with `a>=s`: all future ternary branch words are still
equidistributed, and the policy sees only the residue modulo `M`, so the
conditional mean remains `x_(r mod M)`. This matters because conditioning
on a first branch sends a residue modulo `M` to a single parent residue
modulo `2M`. Applying this refined statement gives the unique linear system

```
3x_r+x_(P_0(r))
  -(1-t_r)x_(P_2(r))-t_r x_(P_1(r)) = 1+t_r.       (2.1)
```

Uniqueness follows because the corresponding affine update is a
supremum-norm contraction with factor `2/3`. This contraction belongs to
the fixed-policy mean system, not to the Bellman operator below.

A quantitative cylinder estimate is also available. Each terminal or
unresolved cylinder has a counting error at most one; there are fewer
than `M*2^k` terminal cylinders and `M*2^k` unresolved cylinders. Hence

```
|A_F(X)-delta X| <= X(2/3)^k+2M*2^k.
```

For `X>=M`, choosing `k=floor(log_3(X/M))` gives the explicit bound
`5 M^(1-c) X^c`, where `c=log_3 2<1`. The ancestor address process here
is not the forward Collatz orbit and must not be interpreted as its
stopping-time distribution.

The all-zero free policy has density `1/3`. Modulo four, choosing `t_2=1`
and all other bits zero gives state means `(8,5,17,5)/27` and density
`35/108`. Modulo 32, the one residues

```
2,10,18,26,30
```

give the exact natural density `8209/25920=0.31670524691358026...`.

## 3. An exact finite policy polytope

Consider the polytope on `0<=x_r<=1` given by

```
3x_r+x_(P_0(r))-x_(P_2(r)) >= 1,
3x_r+x_(P_0(r))-x_(P_1(r)) <= 2.                  (3.1)
```

Every vertex is precisely the solution of one pure coupled policy (2.1).
Indeed `x_r=0` would force `x_(P_0(r))=1`, and `x_r=1` would force
`x_(P_0(r))=0`. Iterating `P_0` reaches zero, whose `P_0` self-loop
cannot alternate between these endpoints. Thus every feasible coordinate
is strictly inside `(0,1)`. Both bounds for a single row cannot be tight,
because this would require `x_(P_2)-x_(P_1)=1`. A vertex needs `M`
independent active constraints, so exactly one bound per row is tight.
These are (2.1), choosing the lower bound for `t_r=0` and the upper for
`t_r=1`. Conversely every such system is uniquely soluble and its
probabilistic solution is feasible, so it gives a vertex.

Thus minimizing the mean coordinate is an exact linear program for this
finite policy class. It is not a universal lower bound on all global
pairings: general assignments can retain residue correlations outside
the declared policy state.

The retained modulus-32 policy has a rational dual certificate. Let
`A_t x=1+t` be (2.1), and solve `A_t^T y=1`. Its dual signs are
`y_r>=0` on lower rows and `y_r<=0` on upper rows. Consequently every
feasible `x` satisfies `sum x>=sum_r y_r(1+t_r)`, with equality at the
retained policy. The script computes both rational vectors exactly and
verifies every sign and equality. No floating LP result is used as a
certificate. It also exhausts all sixteen modulus-four policies.

## 4. The harmonic Bellman operator contracts in Haar L1

Let `r=2/3`, and use probability Haar measure on the binary integers
`Z_2`; for a function periodic modulo `2^s`, its integral is simply the
average over those residues. Define

```
(F d)(2k)   = 1-r d(3k),
(F d)(2k+1) = 1+r min(d(3k+1),0)+r max(d(3k+2),0).    (4.1)
```

This is the scaled harmonic root-gap recursion: child indices are
`(3/2)i+O(1)`, so child harmonic weights have ratio `r` in the leading
term. The additive `O(1)` shifts remain in the actual residue arguments.

**Contraction.** For any `d,e` in `L1(Z_2)`,

```
||F d-F e||_1 <= (2/3)||d-e||_1,
integral F d = 1.                                  (4.2)
```

The even half contributes at most `(r/2)||d-e||_1`. On the odd half,
each map `k->3k+c` preserves Haar measure, and

```
|min(a,0)-min(b,0)|+|max(a,0)-max(b,0)|=|a-b|.
```

It therefore contributes at most the same amount. For the mean, the
negative even integral cancels the sum of the odd minimum and maximum
integrals. This proves both statements.

Starting with `d_0=0`, let `d_h=F^h(0)`. Then `d_1=1`; for `h>=1`,
`d_h` is periodic modulo `2^(h-1)` and has mean one. The contraction
makes the iterates Cauchy in `L1`, with a unique fixed point `d_*` and

```
||d_(h+1)-d_h||_1 <= r^h,
||d_*-d_h||_1 <= 3r^h.
```

The pointwise bound `||d_h||_infinity<=3` follows separately from
`||F d||_infinity<=1+r||d||_infinity`; the fixed point has the same
essential bound. This is not a supremum-norm contraction. The exact
hostile example `d=(2,-1), e=(1,-2)` modulo two has input difference one
everywhere and output difference modulo four

```
(-2/3,4/3,-2/3,0).
```

Its output supremum is `4/3`, while its Haar `L1` norm is exactly `2/3`.

## 5. A universal lower certificate and a matching legal policy upper bound

Take any bounded periodic real potential `d_i` and write

```
R_d=F d-d,
c_0(d)=(integral d+integral min(d,0))/3.
```

Then

```
c_0(d)+integral min(R_d,0) <= B_*
    <= delta(sign policy d)
    <= c_0(d)+integral max(R_d,0).                 (5.1)
```

The sign policy is the legal policy of Section 2 with
`t_i=1` exactly when `d_i<0`, choosing zero at a tie. In particular its
natural density exists.

**Proof of the lower bound and the policy upper bound.** For each legal
parent and its actual children form

```
L_i=epsilon_i/i +sum_(child j) epsilon_j d_j/j -epsilon_i d_i/i.
```

Replace `1/j` by `r/i`; the total error is `O(||d||_infinity/i^2)`.
At fixed parent bit `b`, minimize over its allowed child bits. The
scaled value at `b=0` is

```
Q_0(i)=r d(3i/2)                         if i is even,
Q_0(i)=r min(d((3i+1)/2),0)              if i is odd.
```

The value at `b=1` minus that at zero is exactly `R_d(i)`.
Consequently every legal assignment gives
`L_i >= [Q_0(i)+min(R_d(i),0)]/i+O(i^(-2))`.
The sign policy minimizes every free child potential, for either parent
bit, so it gives equality before the parent bit is selected and hence
the corresponding upper bound with `max(R_d(i),0)`.

Summing over `2<=i<=X` cancels all internal child potentials. The only
remaining terms are a fixed root term and the children with
`X<j<=(3X+1)/2`; their total magnitude is bounded independently of `X`.
The `i^(-2)` errors are summable as well. A periodic harmonic average is
its residue mean times `log X+O(1)`. The mean of `Q_0` is precisely
`c_0(d)`, by the same even/odd Haar changes of variable as above.
Division by `log X` proves the universal logarithmic lower bound and
the sign-policy upper bound. THM-4496 identifies the minimum logarithmic
density as `B_*`; the sign policy's natural and logarithmic densities
agree, proving (5.1).

For `d=d_h`, `h>=1`, the residual has mean zero. Put

```
e_h=||d_(h+1)-d_h||_1,
b_h=(1+integral min(d_h,0))/3.
```

Equation (5.1) becomes the exact computable enclosure

```
b_h-e_h/2 <= B_* <= delta(sign policy d_h) <= b_h+e_h/2,
e_h <= (2/3)^h.                                    (5.2)
```

Sending `h` to infinity identifies the constant:

```
B_*=(1+integral min(d_*,0))/3.                      (5.3)
```

The sign-policy natural densities approach `B_*` from above, with error
at most `e_h`. Every existing natural density, and every upper natural
density, is at least the lower logarithmic density and hence at least
`B_*`. Therefore

```
inf_(global P2 with natural density) natural_density = B_*,
inf_(global P2) upper_natural_density = B_*.
```

The argument through this section alone leaves attainment of either
infimum OPEN; Section 8 supplies a separate proposed construction. This
conclusion must not be changed to a minimum by taking a pointwise sign of an almost-everywhere
fixed point: Haar equivalence classes do not determine values at the
countable set of positive integer addresses. Logarithmic attainment is
already proved in THM-4496, but does not resolve this natural-density
boundary.

## 6. Exact computation and controls

For `h>=1`, store `d_h` as integer numerators with denominator `3^(h-1)`
and period `2^(h-1)`. Formula (4.1) gives the next array using only
integer additions, multiplications by two, minima and maxima. The
residual is computed against the correctly repeated old array on the
doubled period. Thus both endpoints in (5.2) are exact fractions.

The retained depth is `h=24`: the policy period is `2^23`, and the
residual uses all `2^24` residue states. Chunked signed 64-bit NumPy
arithmetic is explicitly bounded before each step; each chunk sum is
converted to a Python integer before accumulation. The resulting exact
interval is the one at the start of this note. The calculation checks
mean one, the contraction envelope and each successive observed `L1`
contraction. An independent fixed-root subtree optimization with
discounted level costs checks every residue through `h=8`, rather than
reusing the root-difference formula.

The same script proves the modulus-32 finite-policy optimum by rational
primal/dual equality, exhausts modulus four, checks all legal clauses
through pair index 20,000, and executes actual two-step descent for
sources 3 through 40,000. The valuation-policy and supremum-norm
counterexamples prevent the two closest invalid generalizations.

Reproduce with
`python3 04-computation/experiments/crossroads_family_20260926_flow.py`,
or with `-O`. The [retained output](crossroads_family_20260926_flow.out)
is raw program stdout. This is a finite computation plus a uniform
convergence proof; the numerical center is not treated as an exact value,
and no original-Collatz claim follows.

## 7. Finite policy LPs are complete in the limit

Let `V_s` be the minimum mean coordinate in the exact policy polytope
(3.1) for modulus `2^s`. Lifting a policy to twice its period proves
`V_(s+1)<=V_s`. Every policy has a natural density, and the sign policy
of `d_(s+1)` has period `2^s`. Consequently (5.2) gives

```
B_* <= V_s <= B_*+e_(s+1) <= B_*+(2/3)^(s+1).       (7.1)
```

Thus the decreasing finite LP optima converge to the global infimum,
with an explicit error. This does not reinterpret a finite LP value as
a universal lower bound: it is an achievable upper bound. Completeness
comes from the separate universal Bellman certificate.

## 8. Slow phase sweep gives natural-density attainment

Status of this section: **PROOF CANDIDATE under independent audit.**
The preceding proved infimum statement is not used as an attainment
claim. This section constructs a single integer-address policy rather
than choosing representatives of an almost-everywhere function.

For each integer `h>=1`, let `P_h` be the sign policy of `d_h`, of period
`M_h=2^(h-1)`. Write `b_n^(h)` for its global legal bits, with common root
bit `b_1^(h)=0`, and `delta_h` for its natural density. The preceding
sections give

```
0 <= delta_h-B_* <= (2/3)^h,
|sum_(n<=X) b_n^(h)-delta_h X|
    <= 5 M_h^(1-c) X^c  for X>=M_h,   c=log_3 2<1.   (8.1)
```

Put `R=3/2`. Define a continuous increasing function `S` on the
nonnegative real line by `S(t)=1` for `0<=t<=1` and by linear
interpolation between `S(h^3)=h` for integers `h>=1`. On its successive
linear pieces the derivative is `1/(3h^2+3h+1)`. In particular,

```
S(t)=t^(1/3)+O(1),
0<S'(t)<=1/7  for t>1,
sup_(u>=T/3) S'(u)=O(T^(-2/3)).                    (8.2)
```

For integer `n>=2`, set

```
t(n)=log_R n,  phi(n)={t(n)},
ell(n)=floor(S(t(n))+1-phi(n)).                    (8.3)
```

Use the free-choice bit of policy `P_(ell(n))` at address `n`. Namely
define `a_n=t_n^(ell(n))`, where `t_n^(h)` is its periodic free-choice
bit, and use the recursion in Section 2 with this varying `a_n`.
Set `epsilon_1=0`. All parents are smaller, so these bits are uniquely
defined and satisfy every legal P2 clause. No density is assumed here.

We prove

```
sum_(n<=X) epsilon_n = B_* X+O(X/(log X)^(1/3)).     (8.4)
```

The estimate is asymptotic with an absolute constant; no small-X or
efficient numerical convergence claim is intended.

### 8A. A patched stationary reference already has density B_*

Define `z_n=b_n^(ell(n))`. This reference sequence need not itself be
legal; it will only be a comparison. Fix large real `X`, put
`T=log_R X`, and restrict first to `sqrt(X)<=n<=X`. The omitted initial
segment has at most `sqrt(X)` terms. Every selected level belongs to an
interval of integers

```
J_X=[(T/2)^(1/3)-O(1), T^(1/3)+O(1)] intersect Z.
```

In particular `|J_X|=O(T^(1/3))`, its minimum tends to infinity, and
`M_max=max_(h in J_X) M_h=exp(O(T^(1/3)))`.

On each unit logarithmic cell `q<=t<q+1`, the function
`G(t)=S(t)+1-{t}` is continuous and strictly decreasing with slope
between `-1` and `-6/7`. The knots of `S` are integers and hence do not
break the interior of such a cell. Its range has length less than one,
so its floor is constant on at most two intervals in the cell. Thus
the selector partitions `[sqrt(X),X]` into `O(T)` real intervals on
which the reference uses one stationary policy. Counting at the two
endpoints and applying (8.1) gives total error at most

```
O(T M_max^(1-c) X^c)=X^c exp(O(T^(1/3)))=o(X/T^(1/3)).
```

All stationary densities used differ from `B_*` by at most
`exp(-Omega(T^(1/3)))`. Endpoint rounding contributes `O(T)`. Therefore
the reference count is `B_*X+o(X/T^(1/3))`.

### 8B. The selector is stable on almost every short ancestor chain

Let `p(n)` be the parent map in Section 2. Until the root is reached,
`p(n)=n/R+u_n` with `|u_n|<=1/3`. Consequently for the `j`-th ancestor

```
|p^j(n)-R^(-j)n|<=1.                              (8.5)
```

Choose `k=floor(T^(1/3))`. For `n>=sqrt(X)` and `0<=j<=k`, all these
ancestors exceed the root for sufficiently large `X`. Taking logarithms
in (8.5) gives, uniformly,

```
t(p^j(n))=t(n)-j+eta_j,
|eta_j|<=eta=O(R^k/sqrt(X)).                       (8.6)
```

The logarithms throughout lie above `T/3`. By (8.2), their values of
`S` differ by at most `O(k T^(-2/3))+O(eta)`. Take
`Delta=C k T^(-2/3)+C eta` with one sufficiently large fixed constant.
Exclude source indices satisfying either

```
dist(t(n),Z)<=eta,
dist(G(t(n)),Z)<=Delta.                            (8.7)
```

Outside the first band no fractional-part wrap occurs in (8.6), so
`{t(p^j(n))}={t(n)}+eta_j`. Outside the second band the resulting change
of `G` cannot cross an integer. Thus

```
ell(p^j(n))=ell(n)  for all 0<=j<=k.                (8.8)
```

On a unit logarithmic cell the first bad band has length `O(eta)`.
Since `G` decreases with slope magnitude at least `6/7`, the second
has length `O(Delta)` and is a union of a bounded number of intervals.
After the change of variable `n=R^t`, the number of bad integer indices
in that cell is `O(Delta R^q)+O(1)`. The same estimate applies to the
last truncated cell. Summing the geometric series over all cells gives

```
number of bad sources <= O(Delta X+T)=O(X/T^(1/3)).  (8.9)
```

Closed bad bands handle equality and endpoints, without changing this
bound. This is the reason for sweeping along logarithmic phase: the
parent reduces `t` by one up to a vanishing error.

### 8C. Resets erase the remaining dependence on earlier choices

For a source outside (8.7), (8.8) says that its actual variable policy
agrees with the stationary policy `P_(ell(n))` for the first `k`
ancestral transitions. If that policy resets during those transitions,
the resulting source bit is independent of all earlier choices. Hence
`epsilon_n=z_n` unless the stationary policy has not reset by depth `k`.

For a fixed level `h`, the CRT count in Section 2 bounds the number of
unresolved sources up to `X` by

```
X(2/3)^k+M_h 2^k.
```

Taking a union over all levels in `J_X` bounds the remaining differences
by

```
|J_X| X(2/3)^k+2^k sum_(h in J_X) M_h
    =o(X/T^(1/3)).                                (8.10)
```

Together with (8.9), the omitted initial segment, and the count of the
reference sequence, this proves (8.4). In particular a single global
P2 pairing attains **natural density `B_*`**, and also attains the
minimum upper natural density. This concerns `B_*`, not the smaller
minimum lower natural density `alpha` from THM-4491.

### 8D. Exact address selection and the failure boundary

No floating-point choice is needed to define the construction. Given
integer `n>=2`, determine `q=floor(log_R n)` by rational power
comparisons. Let `h^3<=q<(h+1)^3`, `D=(h+1)^3-h^3`, and `A=Dq-h^3`.
Then (8.3) has the exact equivalent form

```
ell(n)=h+1  if n^(D-1)*2^A <= 3^A,
ell(n)=h    otherwise.                            (8.11)
```

Indeed the first case is precisely `{t}<=S(t)-h`. The equality convention
agrees with the floor in (8.3). Every required finite policy and every
integer-address bit is therefore explicitly computable.

An abrupt change of policy on height bands does not have the stability
property (8.8); a proportional set of nearby ancestors can cross such
a boundary. Taking a Haar fixed-point representative does not solve
that issue either. The present proof preserves both the integer address
and the logarithmic phase, and uses the uniform reset estimate to remove
the remaining ancestry. These are extra coordinates absent from either
shortcut. Nothing in this construction transfers P2 legality to the
ordinary, unmodified Collatz map.
