# Independent audit: TP2 trace bounds, rational bridges, and the q boundary

**Status:** PROVED audit of the companion argument; section 5 PROVED and
independently audited by the root and flow lanes; FINITE-EXACT controls;
VERIFIED numerical sine comparisons. Date: 2026-09-26. The reviewed source
is [crossroads223_20260926_bridge.md](crossroads223_20260926_bridge.md).
No mathematical gap was found. No Collatz convergence claim is made.

## 1. The trace lemma, including zeros and periodic components

Let P be a finite square nonnegative TP2 matrix. On a selected principal
index set, sorting the columns of any permutation into increasing order
replaces each crossing pair of factors by a weakly larger aligned pair.
All other factors are nonnegative. There is no division, so zero entries
are harmless. Applying this to a simple directed cycle shows that its
edge product is at most the product of its diagonal entries.

A closed walk splits into simple cycles, preserving the multiset of
visited vertices with multiplicity. Its edge product is therefore bounded
by the product of the corresponding diagonal entries. The trace expansion
sums over all ordered n-tuples `(i_0,...,i_(n-1))`, with the return edge
back to i_0. Summing their diagonal products gives exactly `(trace P)^n`.
There is no division by n and no rotation factor:

    trace(P^n) <= (trace P)^n.

The standard finite Perron decomposition gives
`limsup trace(P^n)^(1/n)=r(P)`, hence `r(P)<=trace P`. There is also a
stronger zero-pattern observation specific to this hypothesis: every
vertex on a positive cycle must have a positive diagonal entry. Thus a
nontrivial cyclic irreducible component is aperiodic. Nilpotent and zero
components cause no difficulty. The unrestricted nonnegative matrix
`[[0,1],[1,0]]` is the decisive hostile: trace zero, spectral radius one,
and a negative 2-by-2 minor.

The independent probe checks every binary 3-by-3 matrix. Exactly 147 are
TP2, and all 1,176 comparisons `trace(P^n)<=trace(P)^n`, for n=1,...,8,
pass using integers. This is a finite control, not the proof.

## 2. Rational bridge and indexing audit

For `e=ceil(cL)` and `p=e/L`, the time-j states are `(Z-pj) intersect
(0,M)`. Initial and final states are both `1,...,M-1`, because pL is an
integer. Intermediate state counts can be M or M-1; it would be incorrect
to silently replace the time-dependent rectangular kernels by one fixed
matrix. The probe explicitly forms each grid using integer numerators
`Lx` and multiplies the rectangular matrices.

Each one-step matrix is bidiagonal after the grids are ordered, hence
TP2; products retain TP2. Independently, the note's path-swapping proof
is valid: two crossing paths have integer height difference, changing by
at most one at each time, so they meet at an integer time. Swapping tails
at the first meeting is injective, preserves survival, and preserves the
product of their weights. Zeros and changing intermediate grids remain
within this argument.

For `f(x)=exp(beta_p*x) sin(pi*x/M)`, children removed at the lower or
upper boundary have nonpositive f. Killing therefore gives the required
*lower* inequality, `B_j f>=mu_p f`; this direction is the opposite of
the padded sine supersolution used for the reflected upper bound.
After L steps the phase returns, and the positive initial vector satisfies
`P f>=mu_p^L f`. The Perron lower bound and TP2 trace lemma apply to this
finite square P.

Every trace path has exactly e ones, so all have the same Bernoulli-p
weight. A word has at most M-1 starting heights. Rotating at a minimum
of its rationally centered partial sums keeps those sums in `[0,M)`.
Adding `j(p-c)` makes the c-centered sums strictly positive for j>0,
while their peak stays below M+1. Here `0<L(p-c)<1` is essential.
Forgetting the starting height and rotation loses at most `(M-1)L`
preimages, including imprimitive words and repeated minima. The strict
integer Bad_L condition and the preterminal peak convention both survive.

The probe checks 108 universes: q=3,5,7; L=3,...,11; M=2,...,5.
Its exact matrix traces equal independent enumeration of 17,313 closed
paths. Every produced rotation satisfies `q^(e_j)>2^j` and the strict
peak bound, and every fiber obeys the stated multiplicity. These general-q
finite controls test the bridge map, not a general-q private-price theorem.

## 3. Constants and the polynomial comparison

For q=3 and L>=3, `c<p<=4/5`; the small lengths and the bound
`c<19/30` cover the entire claimed range. The entropy derivative bound
`h'(p)>=-log4` gives exactly the factor 1/4. Since `p>=c>1/2`, every
bracket `1-p^(2j+1)-(1-p)^(2j+1)` decreases with p. Thus `mu_p>=mu_c`
at the same width, without a perturbative error.

Setting `M=m+4` consequently yields

    N_m(L) <= 4e(m+3)L A_(m+5)(L).

The inherited private comparison accepts `K=5` and the L-dependent
constant `K_0=4eL(L+3)`. Its finite theorem is being used within its stated
range; this is not the uniform one-unit Robin comparison.

The sine remainder is bounded correctly by
`t^4/[180(1-t^2/pi^2)] <= t^4/135` for `t<=pi/2`.
With `M=ceil(theta L^(1/3))`, the lower constant `108e` is conservative:
the losses are entropy 4, rounding 9, multiplicity 3, and the exponential
remainder e. The private upper bound plus `pi_L>=2rho_peak(L)` gives the
matching peak upper bound. The proved arbitrary-edit sandwich has only
a polynomial factor. Thus all three claimed q=3 logarithms indeed have
error O(log L); no Mogul'skii citation is needed for this new conclusion.

## 4. Exact parameter boundary and the missing general-q sidecar

For any fixed odd q>=3, write `c=log_q2`, `lambda=(1-c)/c`, and
`D_q=exp[-log2+h(c)]`. The same closed bridge, with sufficiently large L,
does give an elementary lower bound of the form

    log rho_peak,q(L)
      >= L log D_q - kappa_closed,q L^(1/3) - O_q(log L),
    kappa_closed,q=(3/2)(pi^2 c(1-c))^(1/3)(log q)^(2/3).

For c<1/2 the monotonic comparison `mu_p>=mu_c` is unavailable, but this
last statement needs only `p-c=O(1/L)`: the change in the leading sine
coefficient contributes O(1/M^2), and the entropy change is bounded.
The sine remainder remains O(L/M^4). Optimizing M of order L^(1/3)
therefore gives the displayed lower bound.

This is the sharp endpoint regime only when `lambda<=1`, equivalently
`c>=1/2`. In the odd integer multiplier family, q=3 is the case in this
range. For q>=5, `lambda>1`, and the inherited sharp constant uses

    z_q=q/lambda=q*c/(1-c),
    kappa_q=(3/2)(pi^2 c(1-c))^(1/3)(log z_q)^(2/3),

which is strictly smaller than kappa_closed,q. The optimal terminal tilt
then rewards paths ending near the upper wall. The closed bridge ends
within one unit of the bottom and discards that reward.

A trace lower bound gives no bottom-to-top entry bound. Indeed TP2 gives
`P_ii P_jj >= P_ij P_ji`, an upper bound on the off-diagonal product.
Diagonal TP2 matrices show why the required lower bound cannot follow
from TP2 alone. A quantitative killed-kernel lower bound from near one
wall to near the other is an additional sidecar. Section 5 subsequently
supplies it by a separate construction.

One natural repair also fails at the level of the available estimate.
Round the endpoint up by R of order L^(1/3), taking `p-c~R/L`, and use a
closed p-centered bridge of width W. Its c-centered peak is bounded only
by W+R. The entropy gain is `lambda^R` to bounded error, so the resulting
cost exponent is

    W log q + R log(q/lambda) + A_q L/W^2.

Since `q/lambda>1`, this estimate is minimized at R=0. Merely changing
the bridge slope cannot extract the missing top reward. This is a limit
of this bound, not a theorem excluding other endpoint constructions.

There is a second boundary for the private reflected argument: its top
tilted mass is `2(1-c)`. It is sub-Markov only for c>=1/2 and strictly
absorbing only for c>1/2. At c=1/2 the boundary becomes reflecting; for
c<1/2 it has mass greater than one. The q=3 reflected estimate and its
private pairing consequence therefore cannot simply be reused for q>=5.

## Reproduction

    python 04-computation/experiments/crossroads223_20260926_bridge_audit.py

The script uses only the standard library and active checks under `-O`.
Integer controls are FINITE-EXACT. Its 108 sine/trace comparisons use
floating transcendental functions and are labelled VERIFIED. The proof
above and the reviewed note establish the inequalities for all parameters
in the stated domain. No literature priority claim is made.

## 5. Subsequent positive signal: polynomial-cost connectors

**PROVED + INDEPENDENTLY AUDITED by the root and flow lanes.** This construction
addresses the off-diagonal sidecar left open in section 4. It uses a
short rational bridge with a moderately changed slope; it does not append
a forced linear-length word. The previous statements that a trace bound
alone and simple endpoint rerounding are insufficient remain valid.

Fix any c in (0,1), put d=1-c, and use Bernoulli-c increments `bit-c`.
All constants in this section may depend on c. In the Collatz application
c=log_q2 for a fixed odd integer q>=3.

### 5.1 One-sided rational bridges

For integers K>=1, 0<a<K, h>=2, set p=a/K. The trace argument already
audited above gives

    Q_p(0<=S_j^p<h for all j, S_K^p=0)
      >= mu_(p,h)^K/[K(h-1)].                         (B1)

The same lower bound holds with `-h<S_j^p<=0`, by reversing each word.
Indeed the reversed centered prefix is the negative of a terminal
centered prefix of the original bridge. The common-count change of
measure from Q_p to Q_c multiplies either event by

    exp[-K KL(p||c)].

For every 0<=p<=1,
`KL(p||c)<=(p-c)^2/(cd)`, by the elementary bound of relative entropy by
the Bernoulli chi-square divergence. Also the sine series gives, uniformly
for 0<p<1 and h>=2,

    log mu_(p,h) >= -pi^2/(8h^2)-pi^4/(135h^4).       (B2)

No sign assumption on p-1/2 is needed in (B1)--(B2).

### 5.2 A guarded connector lemma

Work inside the real interval (0,M). Set

    h=floor(M/8), K=floor(M^2/log M).

Take M sufficiently large that h>=2 and M/K<min(c,d). Let s,t lie in
`[1,M-1]` and satisfy the phase condition

    t-s+2cK is an integer.                            (B3)

Then the Q_c probability of a 2K-step path from s to t, staying in (0,M),
is at least `M^(-C_c)` for a constant C_c. The same holds for s=0, where
zero is permitted only at time zero.

**Proof.** Choose an intermediate point

    u=s+a-cK in [M/2-1/2,M/2+1/2]

by choosing an integer a nearest to `cK+M/2-s`. The second half has an
integer number of ones by (B3). For a K-step path from s toward u,
write `p'=c+(u-s)/K`. If u>=s, use the nonnegative centered bridge (B1);
the actual path is

    s + S_j^(p') + j(u-s)/K,

which lies between s and u+h. If u<s, use the nonpositive version; its
path lies between u-h and s. Both ranges lie in (0,M). For s=0 the first
case applies and u-s>0, so every positive-time state is strictly positive.

To go from u to t, first construct a path from M-t toward the bulk point
M-u, then reverse time and reflect space. This operation preserves each
increment `bit-c` and reverses the word, preserving its Bernoulli weight.
It gives a path from u to t inside the same interval. The phase condition
guarantees the required integer number of ones.

Each half has `|p'-c|<=M/K`. By (B1)--(B2), its probability is at least

    exp[-pi^2 K/(8h^2)-pi^4 K/(135h^4)-M^2/(cd K)]
       /[K(h-1)].                                    (B4)

Here `K/h^2=O(1/log M)`, `M^2/K=O(log M)`, and `K(h-1)<=M^3`.
Thus (B4) is a reciprocal polynomial in M. Multiplying the independent
halves proves the connector lemma. QED.

### 5.3 Insert a long trace bridge between two connectors

Choose an integer M of order L^(1/3), the K above, and set

    N=L-4K, W=M-4, e=ceil(cN), p=e/N, delta=e-cN.

For the Collatz values c, `0<delta<1`; N is positive and large for all
sufficiently large L. Put `alpha={-2cK}`. For i=1,...,W-1, the first
connector ends at

    x_i=alpha+1+i.

These points have the correct phase after 2K steps and lie in the guarded
interval. Run an N-step trace path in the p-centered strip (0,W), starting
and ending at integer state i. When translated into the c-centered walk,
its state at local time j is

    alpha+1+X_j^p+j(p-c).

It therefore remains inside (1,M-1), and ends at
`y_i=x_i+delta`. Its Q_c weight is its Q_p weight times
`exp[-N KL(p||c)] >= exp[-1/(cd N)]`.

For lambda=d/c>1, choose the final endpoint

    z=floor(cL+M-1)-cL in (M-2,M-1).

The last 2K-step connector from y_i to z has the correct phase, because
`z+cL` and `y_i+c(2K+N)` are integers. Its endpoints are guarded.
For lambda<=1 instead choose `z=ceil(cL+1)-cL in (1,2)`.

For a generic rational c one may have delta=0, alpha=0, or an endpoint
on the guard values 1 or M-1. These are allowed: they are strictly inside
the killed interval (0,M), and the same proof works with closed guard
ranges. For the actual odd-q parameters, c is irrational and the displayed
endpoint ranges are strict.

Concatenating the two connectors and the central trace path gives disjoint
full-word sets for different i: the first 2K letters determine x_i, hence
i. Conditional on i, central words and final connector words multiply
under Q_c. The uniform connector bound can therefore be multiplied by
the *sum* of the central diagonal weights. The TP2 trace lower bound gives

    Q_c(S_j in (0,M), 1<=j<=L; S_L=z)
      >= M^(-C_c') exp[-N KL(p||c)] mu_(p,W)^N.         (B5)

Since `|p(1-p)-cd|<=|p-c|<1/N`, (B2) in its sharper leading-coefficient
form implies

    log mu_(p,W)^N
      >= -A_c N/W^2-pi^2/(2W^2)-pi^4 N/(135W^4),
    A_c=pi^2 cd/2.

For M of order L^(1/3), replacing W=M-4 by M costs O_c(1), while the
fourth-order term is bounded. Consequently (B5) is at least

    M^(-C_c'') exp[-A_c L/M^2-O_c(1)].                 (B6)

The connectors cost polynomial probability although their lengths are
of order M^2/log M; this is the point missed by a forced O(M) suffix.

### 5.4 Sharp general-q lower bound

In the Collatz tilt,

    rho_peak,q(L)=D_q^L E_Qc[lambda^(S_L) q^(-peak); Bad_L].

On (B6), the peak is below M. If lambda>1, the endpoint gives
`lambda^(S_L)>=lambda^(M-2)`; if lambda<=1 it gives a positive constant
depending only on c. Thus, with `z_q=q/max(1,lambda)>1`,

    rho_peak,q(L)
      >= D_q^L M^(-C_q) exp[-M log z_q-A_c L/M^2-O_q(1)].

Choosing `M=ceil((2A_c/log z_q)^(1/3) L^(1/3))` yields

    log rho_peak,q(L)
      >= L log D_q-kappa_q L^(1/3)-O_q(log L),
    kappa_q=(3/2)(pi^2 cd)^(1/3)(log z_q)^(2/3).       (B7)

This proves the sharp general-q lower bound with a logarithmic remainder.
Its additional obligations are exactly the connector
geometry, phases, and disjoint concatenation checked above; it does not
assert a general-q reflected/private-pairing upper bound. A matching
general-q peak upper bound is supplied separately in the next paragraph.

### 5.5 The matching hard-wall upper bound for every c

This uses an unreflected kernel, so it does not inherit the restriction
`2(1-c)<1` of the private reflected construction. For R>=1, the event
`Bad_L` with preterminal peak below R stays in `(0,R+1)` through time L,
because the largest increment is d<1. Put U=R+1, H=U+2=R+3, t=pi/H, and

    f(x)=exp(beta_c*x) sin(t(x+1)),
    beta_c=log[d sin(ct)/(c sin(dt))].

The unrestricted kernel sends f to `mu_c f`. Every killed child has
nonnegative f, because it lies in `(-1,0]` or `[U,U+1)`, within the
positive sine region shifted by one. Therefore `Bf<=mu_c f`. For either
sign of beta_c,

    |beta_c| <= t^2/3,
    |beta_c| U <= pi^2 U/[3(U+2)^2] < 1.

Both distances of x+1 from the sine endpoints are at least one on
`0<=x<=U`; hence `f(x)>=e^(-1)f(0)`. Iteration, followed by the sine
coefficient inequality, gives the uniform estimate

    Q_c(Bad_L, peak<R)
      <= e mu_c^L <= e exp[-A_c L/(R+3)^2].            (B8)

Partition the preterminal peak into bins `[m,m+1)`, m=0,...,L-1. On
Bad_L, `S_L<=peak+d`. Consequently

    lambda^(S_L) q^(-peak)
      <= max(1,lambda^d) z_q^(-m).

Apply (B8) with R=m+1 and set x=m+4. Each term obeys

    z_q^(-m) exp[-A_c L/(m+4)^2]
      <= z_q^4 exp[-kappa_q L^(1/3)].

There are at most L terms, giving the explicit upper bound

    rho_peak,q(L)
      <= e max(1,lambda^d) z_q^4 L D_q^L
           * exp[-kappa_q L^(1/3)].                   (B9)

Together, (B7) and (B9) prove, by the elementary bridge and sine
arguments in this section,

    log rho_peak,q(L)=L log D_q-kappa_q L^(1/3)+O_q(log L)

for every fixed odd q>=3. The same conclusion for the arbitrary-edit
price follows from the inherited THM-4480 polynomial sandwich. The root
and flow lanes independently audited the connector proof; the root also
independently derived this hard-wall upper estimate. No general-q claim
about globally consistent or private pairings follows.

### 5.6 Cheap controls for the connector construction

The independent probe now also checks 15 exact one-sided bridge versus
trace-count inequalities, with M=16,32,64 and five rational slopes at each
width. It constructs 48 phase-compatible connector examples for
q=3,5,7,31, M=32,64, and three central starting heights. The central ceil
inequalities are checked using integer powers; geometric connector states
are numerical probes and labelled VERIFIED. These examples check the
time/space reversal, endpoint phase, and constant-width guards. They do
not replace the uniform probability estimate (B4). For all 48 connectors,
an exact DP counts the one-sided half-word families; their Q_c probability
is numerically compared with (B4), with minimum logarithmic margin greater
than 41. The probabilities are labelled VERIFIED because c is evaluated
numerically. Normal and `-O` output agree.
