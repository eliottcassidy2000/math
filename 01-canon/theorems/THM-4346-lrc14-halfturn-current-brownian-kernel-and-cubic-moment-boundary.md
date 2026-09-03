---
id: THM-4346
title: "LRC(14) half-turn current Brownian kernel, divisor-square filtration, and cubic-moment boundary"
status: >
  PROVED + VERIFIED-EXACT; LRC(14) OPEN. The signed lower/upper-sheet
  occupation current of two odd speeds has an exact gcd-reduced mod-14
  covariance. Its seven-residue numerator has rank three and, after the
  order 1,5,3, is twice the three-step Brownian covariance min(i,j), with
  opposite residues carrying opposite orientation and residue 7 null. For
  an arbitrary odd family, its full energy is the square norm of an exact
  Dirichlet convolution and splits into literal reverse-martingale
  differences indexed by the 7-adic valuations of the speeds.
  Separately, two explicit exchangeable twelve-label sheet laws obey the
  exact one-label LRC marginals and exclusion law and agree on every labelled
  occupation tensor through total degree three, yet have zero-minimum masses
  5/28 and 0. Thus polynomial cubic occupation data alone cannot force a
  free sheet without an arithmetic/geometric sidecar; degree four is the
  first separating moment for this pair. The laws are abstract and the exact
  current kernel proves that neither is realizable by runner intervals.
source: root + cubic_depth + current_hierarchy / LRC14 continuation session, 2026-09-02
depends_on: []
related:
  - THM-638-signed-pair-mass-law-rational-thresholds
  - THM-594-pair-overlap-law-mirsky-newman-floor
  - THM-4345-lrc14-halfperiodic-anchor-strip-euclidean-remainder-and-current-envelope
  - THM-2236-pointwise-nested-binomial-minorants-and-cubic-vertex-fan
  - HYP-3223-lrc14-green-current-lorentzian-exchange-angles
  - THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal
  - THM-4338-lrc14-all-pair-rank-four-cubic-majorant-and-uniform-one-appender
artifacts:
  - 04-computation/lrc14_halfturn_current_kernel_cubic_blindness_probe_20260902.py
  - 05-knowledge/results/lrc14_halfturn_current_kernel_cubic_blindness_probe_20260902.out
  - 04-computation/lrc14_halfturn_cubic_depth_exact_audit_cubic_depth_20260902.py
  - 05-knowledge/results/lrc14_halfturn_cubic_depth_exact_audit_cubic_depth_20260902.out
current_script_sha256: eff85b6fe14fbe225bd116c47aa1927b2c4fa03d7e7de0a6f214c6e0e8003f98
current_output_sha256: 73abb0129bc39f4f0d8048ccf13787ccd3a4adc33d2931b7a08da8359aba7e64
cubic_audit_script_sha256: c35afe7fdf2ff40c5a8ccbb84e1e34370472e98a1b536d15f01171966d8572e1
cubic_audit_output_sha256: c4acd92d448b8667fe5004cc145edd51408bb297e88cda714cf1817ebac71d9e
hash_basis: raw LF bytes
audit: >
  PASS. A literal rational wall sweep on the normalized half-base checks all
  2,500 ordered odd pairs u,v<100 against the closed covariance formula. An
  exact formal-sine audit checks the Dirichlet-convolution coefficients through
  frequency 315 on all 256 subsets of the odd speeds through 15, as well as
  7-adic shell energies, Fourier projectors, private fundamental modes, and
  the exact twelve-term 9-power hostile. An independent exact finite-law
  calculation checks all factorial moments through degree three, the
  fourth-difference stencil, the free masses, and the Brownian residue
  factorization. The analytic proofs below do not depend on the finite ranges.
---

# THM-4346 -- half-turn current, Brownian residues, and cubic blindness

**PROVED + VERIFIED-EXACT. LRC(14) REMAINS OPEN.**

## 1. The signed half-turn current

Put

```text
A(x)=1_{||x||<1/14},
sigma_w(t)=A(wt)-A(w(t+1/2))                         (1)
```

for a positive odd integer `w`. Endpoints have measure zero and may be
assigned arbitrarily. For an integer `n`, define

```text
d_14(n)=min(n mod 14, 14-(n mod 14)).                (2)
```

For positive odd `u,v`, let

```text
g=gcd(u,v),   U=u/g,   V=v/g.                        (3)
```

Then `U,V` are coprime and odd, and

> **Exact current kernel.**
>
> ```text
> integral_0^1 sigma_u(t)sigma_v(t)dt
>   =[d_14(U+V)-d_14(U-V)]/(7UV).                    (4)
> ```

Because each `sigma_w` changes sign under `t -> t+1/2`, the product in
`(4)` is half-periodic. Thus `(4)` is also the integral on `[0,1/2)` with
normalized measure `2dt`.

The formula includes the diagonal: setting `u=v` gives `2/7`, the measure
on which either one of the two sheets is dangerous for that runner.

### Proof

The nonzero Fourier coefficients of `A` are

```text
hat A(k)=sin(pi k/7)/(pi k).                          (5)
```

The half-turn difference kills even frequencies and doubles odd ones.
In the product integral the only surviving frequency pairs are

```text
k=Vr,   l=-Ur,   r odd.                              (6)
```

Writing

```text
F(theta)=sum_(r>=1, r odd) cos(r theta)/r^2,
```

the standard cosine series, obtained by subtracting the even terms from
the full quadratic cosine series, is

```text
F(theta)=pi^2/8-(pi/4) dist(theta,2pi Z).             (7)
```

Product-to-sum applied to `(5)--(7)` gives

```text
4/(pi^2 UV) [F(pi(U-V)/7)-F(pi(U+V)/7)],             (8)
```

which is exactly `(4)` after reducing the two angular distances modulo
`2pi`. This proves the formula for all positive odd `u,v`.

## 2. The residue kernel is a Brownian covariance

In the conventional residue order

```text
R=(1,3,5,7,9,11,13),                                 (9)
```

the numerator `M(r,s)=d_14(r+s)-d_14(r-s)` is

```text
     1  3  5  7  9  11 13
1    2  2  2  0 -2  -2 -2
3    2  6  4  0 -4  -6 -2
5    2  4  4  0 -4  -4 -2
7    0  0  0  0  0   0  0
9   -2 -4 -4  0  4   4  2
11  -2 -6 -4  0  4   6  2
13  -2 -2 -2  0  2   2  2.                         (10)
```

Define an orientation `epsilon` and a Brownian time `ell` by

```text
r          1   5   3   13   9   11   7
epsilon    +   +   +    -   -    -   0
ell        1   2   3    1   2    3   0.             (11)
```

Then the whole table has the exact factorization

```text
M(r,s)=2 epsilon(r)epsilon(s) min(ell(r),ell(s)).     (12)
```

Thus its rank is three. In the positive order `(1,5,3)` its block is

```text
2 [min(i,j)]_(1<=i,j<=3)
 =[[2,2,2],[2,4,4],[2,4,6]],                         (13)
```

twice the covariance matrix of Brownian motion sampled at times `1,2,3`.
The inverse block is the path precision matrix

```text
[[1,-1/2,0],[-1/2,1,-1/2],[0,-1/2,1/2]].            (14)
```

Equivalently, `(12)` is the Gram factorization by the three nested signed
features

```text
epsilon(r) 1_{ell(r)>=j},   j=1,2,3.                 (15)
```

This is a literal covariance normal form, not an analogy.

For example, let `W` be `m` pairwise-coprime positive odd speeds and set

```text
a(t)=sum_(w in W) A(wt),
b(t)=sum_(w in W) A(w(t+1/2)),
C(t)=a(t)-b(t).                                      (16)
```

Writing `epsilon_w=epsilon(w mod 14)` and `ell_w=ell(w mod 14)`, equations
`(4)` and `(12)` give the exact stochastic variance identity

```text
E C^2=(2/7) [m-sum_w ell_w/w^2
                 +sum_(j=1)^3 (sum_(w:ell_w>=j) epsilon_w/w)^2]. (17)
```

For a general odd family, `(4)` remains exact, but each pair must use its
own gcd-reduced residues; the global three-square compression `(17)` need
not survive.

## 3. The arbitrary-gcd Dirichlet square and a septimal martingale

The loss of the three-square Brownian compression does not destroy
positivity.  It moves positivity from three residue coordinates to the full
divisor lattice.  Define, for a positive integer `k`,

```text
beta(k)=2 sin(pi k/7)/(pi k)   if k is odd,
       =0                      if k is even,          (17a)
```

and let `1_W` be the indicator of the finite odd-speed family `W`, regarded as
an arithmetic function.  The Fourier coefficient of `C` from `(16)` at a
positive frequency `n` is

```text
C_hat(n)=sum_(w in W, w|n) beta(n/w)
        =(1_W *_D beta)(n),                           (17b)
```

where `*_D` is Dirichlet convolution.  All even coefficients and the constant
coefficient vanish.  Parseval therefore gives the exact arbitrary-gcd square

```text
integral_0^1 C(t)^2 dt
 =2 sum_(n>=1, n odd) |(1_W *_D beta)(n)|^2
 =8/pi^2 sum_(n>=1, n odd) 1/n^2
       [sum_(w in W, w|n) w sin(pi n/(7w))]^2.        (17c)
```

This is the global positive completion of the pair kernel `(4)`: expanding
the squares and summing common frequencies recovers `(4)` for every ordered
pair, with no coprimality assumption on `W`.

There is also a literal stochastic-process structure.  Since
`sin(pi k/7)=0` exactly when `7|k`, every Fourier mode contributed by
`sigma_w` has

```text
nu_7(n)=nu_7(w).                                     (17d)
```

Put

```text
W_e={w in W: nu_7(w)=e},       C_e=sum_(w in W_e) sigma_w. (17e)
```

The `C_e` have disjoint Fourier supports, hence

```text
C=sum_(e>=0) C_e,
integral C^2=sum_(e>=0) integral C_e^2.              (17f)
```

For completeness, define the averaging operators

```text
(Q_e f)(t)=7^(-e) sum_(j=0)^(7^e-1) f(t+j/7^e),     (17g)
```

with `Q_0` the identity.  On Fourier mode `n`, `Q_e` is multiplication by
`1_{7^e|n}`.  Equivalently, it is conditional expectation onto the decreasing
sigma-fields generated by `7^e t mod 1`.  Since every current has mean zero,

```text
Q_e C=sum_(j>=e) C_j,
(Q_e-Q_(e+1))C=C_e.                                 (17h)
```

Thus `(17f)` is exactly the orthogonality of a finite reverse-martingale
difference sequence, not an electrical analogy.  The Bohr/Dirichlet-series
view says the same thing: `beta` has no coefficient carrying the prime `7`,
so the exponent of the `7`-coordinate grades the product orthogonally.

The decomposition gives a small unconditional floor.  In each `W_e`, let
`M_e` be the elements minimal under divisibility, and put
`m_7(W)=sum_e |M_e|`.  At frequency `n=w` for `w in M_e`, no other speed in
the same shell contributes.  A divisor from a lower shell has quotient
divisible by `7` and hence contributes zero.  Therefore

```text
C_hat(w)=beta(1)=2 sin(pi/7)/pi,
integral_0^1 C^2
 >=m_7(W) * 8 sin(pi/7)^2/pi^2.                     (17i)
```

This is the current-specific, shell-refined version of the divisor-minimal
Fourier mechanism in the surviving Parts A and C--E of THM-594.  It is an
exact private-frequency contribution, not a claim that `(17i)` is the sharp
global energy infimum.

The new positivity also has a sharp stopping signal.  For

```text
W_9={1,9,9^2,...,9^11},                              (17j)
```

all speeds lie in the single shell `e=0`.  If two exponents differ by `d>0`,
then `9^d mod 14` cycles through `9,11,1`, so `(4)` gives pair covariance
`2s_d/(7*9^d)`, where `s_d` cycles through `-1,-1,+1`.  Consequently

```text
integral_0^1 C_(W_9)^2
 =585613119200/219667417263
 ~=2.66590797350 <4.                                (17k)
```

More generally the energy per speed of the initial `m`-term 9-power chain
tends to

```text
2/7+(4/7) sum_(d>=1) s_d/9^d
 =2/7+(4/7)(-89/728)=275/1274.                      (17l)
```

Thus the septimal martingale can be trivial on a natural hostile family, and
the full quadratic current energy itself can lie below the threshold `4` in
the variance rebate `(35)`.  On the ordinary half-base its total energy is
half of `(17k)`, already below `12/7`; restricting to any even-anchor safe
core only decreases it.  Hence no anchor can make the *variance-only* current
correction positive on this family.  The exact nonlinear current profile can
still contribute.

Finally, restricting `(17c)` to an anchor strip multiplies by the strip
indicator and couples distinct Fourier modes; the diagonal square sum is no
longer the answer.  THM-4345 identifies the missing datum exactly as a
primitive residual-wall integral.  Therefore `(17c)--(17h)` supply a clean
full-circle filtration, but do not erase the anchor-address sidecar.

## 4. An exact cubic-moment blindness pair

Consider twelve abstract labels and two sheets. Conditional on a depth pair
`(a,b)`, choose uniformly an ordered pair of disjoint label subsets of sizes
`a,b`. Define the following sheet-symmetric laws:

```text
Law P: (0,4) 5/56, (4,0) 5/56, (2,2) 15/28, (1,1) 2/7;
Law Q: (1,3) 5/14, (3,1) 5/14,                 (1,1) 2/7. (18)
```

Both are probability laws. Every label has marginal `1/7` on each sheet,
and no label can occur on both sheets. These are precisely the universal
one-label marginal and opposite-sheet exclusion laws of twelve odd runner
intervals on the normalized half-base.

For every `r,s>=0` with `r+s<=3`, the factorial moments agree:

```text
E[binom(a,r)binom(b,s)]

(r,s)  (0,0) (1,0) (0,1) (2,0) (0,2) (1,1)
value      1   12/7  12/7  15/14 15/14  17/7

(r,s)  (3,0) (0,3) (2,1) (1,2)
value    5/14  5/14  15/14 15/14.                   (19)
```

Uniform conditional labelling turns these factorial moments, after division
by the corresponding falling-factorial denominator, into every labelled
same-sheet or mixed-sheet intersection probability through rank three.
Repeated same-sheet labels reduce to a lower rank, while a repeated label on
opposite sheets has probability zero. Hence the complete labelled tensor
through rank three is identical under `P` and `Q`.

Nevertheless,

```text
P(min(a,b)=0)=5/28,       Q(min(a,b)=0)=0,            (20)
E_P[1_{a=0}+1_{b=0}]=5/28,
E_Q[1_{a=0}+1_{b=0}]=0.                              (21)
```

The signed difference `P-Q` on the five points `(a,4-a)`, `a=0,...,4`, is

```text
(5/56) [1,-4,6,-4,1].                                (22)
```

This is the fourth finite-difference stencil. It annihilates every polynomial
of total degree at most three and first fires in degree four. In particular,

```text
Delta E[binom(a,2)binom(b,2)]=15/28.                 (23)
```

The two laws even have the same signed-current second moment,

```text
E_P(a-b)^2=E_Q(a-b)^2=20/7,                          (24)
```

whereas their fourth current moments are `320/7` and `80/7`. Thus the first
current-moment separator for this pair is also quartic.

The exact kernel also diagnoses the deliberately discarded arithmetic. Under
uniform labelling, every distinct pair of labels has signed-current covariance

```text
2*(15/7)/(12*11)-2*(17/7)/(12*11)=-1/231.            (25)
```

No pair of positive odd runner speeds has this covariance. Indeed `(4)` would
force

```text
M/(7UV)=-1/231, hence 33M=-UV,                       (26)
```

but the numerator `M` in `(10)` is even while `UV` is odd. Thus both laws are
provably non-realizable already after restoring the mod-14 arithmetic kernel.
This does not distinguish `P` from `Q`; it identifies the exact sidecar that
the abstract moment model omitted.

### Proof

The marginal and moment claims follow by direct substitution in `(18)`.
Conceptually, the common `(1,1)` atom cancels, and every degree-at-most-three
polynomial restricted to `a+b=4` becomes a one-variable cubic; `(22)` kills
it. The fourth difference of `a^4` is `24`, so degree four separates. The
uniform-labelling statement follows from

```text
P(I subset left, J subset right | a,b)
  =(a)_|I| (b)_|J|/(12)_(|I|+|J|)                   (27)
```

for disjoint fixed label sets `I,J`.

## 5. Where the nonlinear `min` certificate gets its signal

Let

```text
r(x)=1-x+binom(x,2)-binom(x,3)
    =-(x-1)(x-2)(x-3)/6.                            (28)
```

On the nonnegative integers, `r` is nonincreasing, with the plateau
`r(1)=r(2)=r(3)=0`. Therefore, for `q=min(a,b)`,

```text
r(q)=max(r(a),r(b)).                                 (29)
```

Under every sheet-symmetric law,

```text
E r(q)=E r(a)+(1/2)E|r(a)-r(b)|.                    (30)
```

If `S=a+b` and `C=a-b`, direct factorization gives

```text
|r(a)-r(b)|
 =|C|[a^2+ab+b^2-6(a+b)+11]/6
 =|C|[C^2+3(S-4)^2-4]/24.                           (31)
```

The bracket in the first line is nonnegative for unequal nonnegative integer
depths; it vanishes exactly on the unequal pairs inside the plateau
`{1,2,3}`. Consequently the familiar cubic lower certificate becomes

```text
E[1-q+binom(q,2)-binom(q,3)]
 =E[1-a+binom(a,2)-binom(a,3)]
  +(1/48)E{|C|[C^2+3(S-4)^2-4]}.                    (32)
```

The first term is the one-sheet third Bonferroni truncation. The second is a
nonlinear cubic-current rebate. For both laws in `(18)` the first term is
zero. Law `P` obtains `5/28` entirely from the rebate, while law `Q` remains
on its zero plateau. This explains exactly how a pointwise cubic in
`min(a,b)` can distinguish two laws whose polynomial moment tensors through
degree three agree: the `min`, equivalently the absolute current in `(30)`,
is the missing nonlinear coupling coordinate.

On the twelve-label state simplex `a+b<=12`, exact integer comparison also
gives

```text
|r(a)-r(b)| >= max((a-b)^2-4,0)/6.                   (33)
```

Let `E` be any half-turn-invariant region and integrate with its ambient
(not necessarily normalized) measure. Put

```text
B_3(E)=integral_E r(a)=integral_E r(b),
V(E)=integral_E C^2.                                  (34)
```

Equations `(30)` and `(33)`, followed by
`max(C^2-4,0)>=C^2-4`, yield the rigorous variance certificate

```text
integral_E r(min(a,b))
 >= B_3(E)+[V(E)-4|E|]/12.                            (35)
```

For the anchor-safe half-base `E=G_{2h} intersect [0,1/2)`, one has
`|E|=3/7`, so the correction is `[V(E)-12/7]/12`.  Formula `(4)` computes
the full-half variance exactly.  Writing `D_{2h}` for the complementary
anchor-danger strip,

```text
V(E)=V([0,1/2))-V(D_{2h} intersect [0,1/2)).          (36)
```

Thus the stochastic bridge has one sharply named missing sidecar: an upper
bound for the signed-current energy absorbed by the anchor strip.  Full-
circle Brownian covariance alone does not supply that bound.

## 6. Scope and consequence

The pair `(18)` is an information-theoretic obstruction, not an LRC(14)
counterexample. It proves only the following sharp limitation:

> The universal one-label sheet laws plus the complete polynomial labelled
> intersection tensor through rank three do not determine, and cannot by
> themselves force, positive zero-minimum mass.

An arithmetic interval constraint, an address/renewal condition, a nonlinear
current statistic, or a genuine rank-four coordinate can still rule out Law
`Q`. This does not conflict with THM-4338: that theorem uses a cubic coverage
majorant *after restricting to rank-at-most-four cells*, so rank four is
already present in its state space. Equations `(12)` and `(32)` identify two
compact sidecars worth carrying into the entry problem: three signed residue
currents for quadratic imbalance, and the absolute-current rebate for the
first nonlinear separation.

## 7. Reproduction

```bash
python3 -B 04-computation/lrc14_halfturn_current_kernel_cubic_blindness_probe_20260902.py
```

The frozen output is
`05-knowledge/results/lrc14_halfturn_current_kernel_cubic_blindness_probe_20260902.out`.
