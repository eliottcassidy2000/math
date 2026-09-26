# One orbit's dyadic crossings: a finite-carry potential and a Sturmian carrier

**Status: PROVED scoped identities and conditional consequences, with
FINITE-EXACT controls.** The consequence for an infinite nonperiodic orbit
uses the proved reciprocal-summability theorem named below. It does not
prove that such an orbit exists or exclude one. No literature-priority
claim. This is a different observable from the owner's three-colour word.

## 1. Inheritance, board, and connection contract

The closest proved mechanism is the ordered-carry logarithmic identity in
[the discrepancy note, Section 1](collatz_guards_20260921_discrepancy.md),
combined with [THM-4476, thin divergent orbits](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md).
The latter proves that a nonperiodic positive integer orbit has summable
reciprocals. The canonical hostile is a completed excursion that grows:
`27 -> 41 -> 62 -> 31`. Both endpoints have dyadic height 4, but 31>27.
The corrected near miss is endpoint-free HYP-9161, refuted by the actual
hover/drop construction in [the integer audit, Section 5](crossroads_poset_20260926_integer.md).
The least-used sidecar here is the exact winding of logarithmic mantissa,
retained jointly with the valuation depth of every actual odd transition.

The five live concepts are (i) dyadic flux, (ii) laminar excursions,
(iii) the finite additive carry, (iv) irrational-rotation languages, and
(v) the low-bit valuation clock. The anchor is an orbit-coupled potential;
the niche is an exact finite-language restriction; the wildcard is a
Fourier potential for the proposed log-periodic trigonometric functions.

| source -> target | map and preserved predicate | destroyed information | sidecar / decisive test |
|---|---|---|---|
| one integer orbit -> dyadic path | `h(n)=floor(log2 n)`; exact up/down crossings | mantissa and parity | retain odd source and `v2(3n+1)`; test growing return 27->31 |
| dyadic path -> excursion forest | stack-match each down-step to latest unmatched up-step | arithmetic values and chronology of disjoint branches | retain all edge labels; generic balanced paths are insufficient |
| odd orbit -> rotating phase | `phi={log2 x}` and accumulated carry | total even depth | retain `a=v2(3x+1)`; no independence assumption |
| finite phase word -> Sturmian language | exact rational order of phase cuts | initial integer realization and unbounded valuations | finite-height condition below; source-height cutoff remains present |
| poset -> uniform extensions | would randomize incomparable events | the deterministic law of this actual orbit | no such law established; the 27/31 covariance counterexample still applies |

In particular, a laminar forest can be ordered by containment, but that
does not make the observed chronology a uniformly sampled linear
extension. Its existence is a legitimate poset-adjacent construction;
an application of the 1/3--2/3 balancing theorem would need an additional
measure-preserving map or an inequality valid for its actual weights.

## 2. Exact crossing flux and the finite-carry potential

Use the shortcut map `T(n)=n/2` on even n and `(3n+1)/2` on odd n. For its
successive odd values write

```
x_(j+1)=(3*x_j+1)/2^a_j,       a_j=v2(3*x_j+1)>=1,
h_j=floor(log2 x_j),           phi_j=log2 x_j-h_j,
alpha=log2(3/2),               c_j=log2(1+1/(3*x_j)),
u_j=h((3*x_j+1)/2)-h(x_j) in {0,1},
e_j=a_j-1.
```

Here u_j counts the up-crossing at the odd half-step, while e_j counts
the following even halvings, each of which crosses exactly one dyadic
level downwards. An odd half-step increases a positive integer by a
factor between 3/2 and 2, so u_j is 0 or 1. The exact up-crossing test is
`3*x_j+1 >= 2^(h_j+2)`.

The first central identity is

```
u_j = phi_j + alpha + c_j - phi_(j+1),             (2.1)
h_(j+1)-h_j = u_j-e_j.                            (2.2)
```

Even halvings preserve mantissa, which proves (2.1); counting their
heights proves (2.2). Thus, for J odd steps, with `U_J=sum u_j`,
`E_J=sum e_j`, and `C_J=sum c_j`,

```
U_J-alpha*J = C_J+phi_0-phi_J,                     (2.3)
h_J-h_0 = U_J-E_J,
log2(x_J)-C_J = log2(x_0)+J*log2(3)-sum a_j.       (2.4)
```

This is a potential on one fixed orbit, not an average over source
residues. Its error is the actual positive additive carry. It is exact
even for the terminal odd fixed point x=1: then u=1, a=2 and
`c=log2(4/3)=1-alpha`, so no spurious negative drift is obtained.

**Conditional theorem.** Suppose the actual positive integer orbit is
not eventually periodic. Its states are distinct and tend to infinity.
THM-4476 gives `sum_j 1/x_j<infinity`, hence

```
0<=C_J<=C_infinity <= (1/(3*ln2))*sum_j 1/x_j < infinity.
```

Consequently up-crossings have bounded prefix discrepancy:

```
|U_J-alpha*J| < C_infinity+1.                      (2.5)
```

For an interval beginning at i, the same error is bounded by the tail
carry `sum_(j>=i)c_j+1`, which tends to 1. The removed information is
precisely E_J, the total depth of the even descents. Equation (2.2) shows
why a theorem controlling U alone does not control net escape.

## 3. A finite-height Sturmian-language theorem

This section does not require an infinite orbit or summability. For any
fixed L, let `S_L(alpha)` be the length-L words

```
w_j=floor(theta+(j+1)*alpha)-floor(theta+j*alpha),
0<=j<L, 0<=theta<1.
```

There are exactly L+1 such words. Irrationality of alpha follows because
`(3/2)^q=2^p` would contradict unique prime factorization.

Here is an exact computable sufficient height bound for actual crossing
words to belong to this language. In multiplicative phase coordinates
`t=2^theta`, the L cuts are

```
b_j = 2^ceil(log2(3^j))/3^j,  1<=j<=L,            (3.1)
```

together with the endpoints 1 and 2. Sort these rational numbers as
`1=b_(0)<...<b_(L+1)=2`, and put

```
rho_L=min_i b_(i+1)/b_(i)>1.                      (3.2)
```

**Finite-height theorem.** If L successive odd sources all satisfy
`x_j>=m` and

```
(1+1/(3*m))^L < rho_L,                            (3.3)
```

their actual up-crossing word belongs to `S_L(alpha)`.

Proof. Let `P_s=sum_(j<s)c_j`; then the prefix crossing counts are
`floor(phi_0+s*alpha+P_s)`. The corresponding cut b_s is moved to
`b_s/2^P_s`, where `1<=2^P_s<rho_L`. Every cut remains strictly between
its old predecessor and its old location. Hence none crosses another
cut or either fixed endpoint; their circular order remains unchanged.
The floor values at initial phase zero remain unchanged as well. As
initial phase increases, each cut changes exactly one prefix floor by
one, in the same order as before. Thus the set of prefix-floor vectors,
and therefore their successive-difference words, is unchanged. Boundary
values use the corresponding right-hand interval and introduce no extra
word. This proves the assertion.

All quantities in (3.1)--(3.3) admit exact rational comparisons. For L=1
through 16, the smallest integer m certified by this sufficient test is

```
2, 6, 9, 12, 32, 39, 45, 52, 58, 64, 71, 296, 320, 345, 369, 394.
```

These are sufficient cut-stability thresholds, not claimed optimal
thresholds. The jump at L=12 records a closer return of the rotation;
it is a concrete appearance of approximation of log2(3/2), not the
golden-ratio rotation underlying the Fibonacci word.

Equivalently, `log2(rho_L)=min_(1<=d<=L)||d*alpha||`, the minimum
circular separation of the L+1 phase cuts including zero. In particular,
a non-Sturmian observed length-L crossing block **certifies that at least
one of its odd sources is smaller than the certified m**. This is an
exact finite certificate for a visit to a bounded arithmetic core.

**Consequence for one hypothetical divergent integer orbit.** For every
fixed L, all sufficiently late length-L crossing blocks belong to
`S_L(alpha)`, since the odd states tend to infinity. This conclusion in
fact uses only escape to infinity, not the stronger reciprocal bound.
It does not assert that one fixed index starts an exactly Sturmian
infinite tail; the cutoff may depend on L. Interchanging those two
quantifiers would lose an essential condition.

Two cheap local consequences are useful before any asymptotic claim:

* `00` is impossible on every ordinary positive orbit. After u_j=0,
  `phi_(j+1)=phi_j+alpha+c_j>=alpha`; since `2*alpha>1`
  (equivalently 9>8), the next u is 1.
* `111` is impossible if its three odd sources are all at least 7.
  Three up-crossings would require
  `3 < 1+3*alpha+sum c_j`. But
  `3*alpha+3*log2(22/21)=log2((11/7)^3)<2`, because
  `11^3=1331<4*7^3=1372`. This contradicts that requirement.
  The small core matters: `3,5,1` gives `111`.

The owner's supplied blue-gap word has both 00 and 11 after encoding
gaps 2 and 4 as 0 and 1, and is therefore not a binary Sturmian factor.
That earlier exact obstruction is retained in
[the energy audit, Section 6](collatz_blueprint_20260921_energy.md).
The present crossing word has a specified arithmetic observable, a
different slope, and a theorem explaining its restrictions; it does not
repair the colouring by silently changing the listed data.

## 4. One-orbit Benford law and a bounded Fourier potential

Under the nonperiodic-orbit hypothesis of Section 2, (2.1) gives

```
phi_J = {phi_0+J*alpha+C_J},  C_J -> C_infinity.    (4.1)
```

Thus the phases asymptotically shadow one fixed irrational rotation
with intercept `phi_0+C_infinity`. For every `1<=a<b<=2`, the fraction
of odd orbit values with binary mantissa in [a,b) tends to
`log2(b/a)`. This is a base-2 Benford law along this very orbit, not a
statement obtained by putting a Haar law on starting integers. It is
conditional on the existence of the nonperiodic orbit.

There is also a quantitative bounded-sum statement for every fixed
nonzero integer Fourier mode m. Set `z_j=exp(2*pi*i*m*phi_j)` and
`q=exp(2*pi*i*m*alpha)`, with q!=1. Then

```
z_(j+1)=q*z_j*exp(2*pi*i*m*c_j),
(q-1)*sum_(j<J)z_j = z_J-z_0
                    -q*sum_(j<J)z_j*(exp(2*pi*i*m*c_j)-1).
```

Therefore

```
|sum_(j<J) exp(2*pi*i*m*log2(x_j))|
 <= (2+2*pi*|m|*C_infinity)/|exp(2*pi*i*m*alpha)-1|. (4.2)
```

The right side is independent of J. Dividing by J proves vanishing
Fourier averages; approximating interval indicators by continuous
trigonometric polynomials proves the stated mantissa distribution.
This supplies an elementary derivation without importing a probabilistic
mixing assumption. A mode-by-mode extension to an infinite Fourier
series needs convergence against the displayed small denominators.

Equivalently, for `f(theta)=exp(2*pi*i*m*theta)`, the bounded phase
potential `psi(theta)=f(theta)/(q-1)` satisfies
`psi(theta+alpha)-psi(theta)=f(theta)`. Actual orbit steps add only the
summable perturbation controlled in (4.2). Thus a finite trigonometric
log-periodic correction has bounded accumulated centered contribution
on any hypothetical divergent orbit. It cannot by itself supply an
unbounded negative drift. The missing term must couple phase to actual
valuation depth or an equally strong arithmetic coordinate.

### 4A. The same orbit has a second clock in base 3

The odd-event clock is not an arbitrary choice. There is a companion
clock on every shortcut time t. Let `n_t=T^t(n_0)`, let O_t count its odd
sources before t, and let

```
D_t=sum_(i<t, n_i odd) log3(1+1/(3*n_i)),
beta=log3(2), theta_t={log3(n_t)}.
```

Then exactly

```
log3(n_t)=log3(n_0)-t*beta+O_t+D_t,
theta_t={theta_0-t*beta+D_t}.                       (4.3)
```

Under the same nonperiodic-orbit hypothesis, D_t converges. Hence the
full shortcut orbit is Benford in base 3, and its fixed nonzero Fourier
mode sums obey the analogue of (4.2), with rotation -beta and carry D.
The proof is the identical telescoping calculation. The two clocks
observe the same arithmetic process at different sections: the odd clock
has return times a_j in the full shortcut clock. Those return times are
precisely the valuation sidecar that must not be replaced by a random
or independent sampling rule. Two valid phase laws do not automatically
prove an independence law for phase and parity.

This extends the inherited ordered-carry calculation, rather than
claiming an unrelated new dynamical conjugacy. The root session is
checking primary literature for related log-coordinate constructions.

## 5. Excursion posets, exact flux, and the failed contraction

For the full shortcut orbit, dyadic height moves by -1,0,+1. Match each
down-step to the latest unmatched up-step, if one exists. A matched pair
delimits an excursion above its starting level. Such intervals are
disjoint or nested, so inclusion gives a rooted forest/laminar poset.
The unmatched steps record endpoints and net change of height. This is
an intrinsic poset on actual orbit events; it needs no cosmetic
tournament or invented tie-breaking relation.

On the complete 70-step shortcut orbit of 27 to 1, there are 25 matched
excursions, no unmatched up-steps, and four unmatched down-steps, agreeing
with the change from height 4 to height 0. Eight matched excursions
increase their return value; examples are 27->31 in 3 steps and 47->61
in 52 steps. Thus every-excursion contraction is false even after the
excursion is fully closed. Relabelling an increased return as a new
branch does not remove its cost.

A containment poset's incomparability means two excursions are disjoint,
not that one may reorder their arithmetic maps while keeping the same
source. The affine intercept depends on chronological order. Consequently
the 1/3--2/3 balance property or the new large-width results cannot be
applied to the observed ordering without a proved sampling law. The
[previous exact poset/height audit](crossroads_poset_20260926_bridge.md)
already shows that even a source-height restriction can reverse an XYZ
covariance sign. This remains the cheapest hostile for proposed transfers.

### 5A. Every finite permitted crossing word coexists with pure growth

**Finite saturation theorem.** Fix L and any word in `S_L(alpha)`. There
are arbitrarily large positive integers realizing this crossing word
with `a_j=1` for every one of its L odd steps. Thus E_L=0 throughout that
block and its source grows strictly at each step.

Proof. Each rotation word has a phase interval of positive length.
Choose a compact subinterval of its interior. Integers
`n=-1 mod 2^(L+1)` realize L successive valuations equal to one, with

```
x_j=(3/2)^j*(n+1)-1, 0<=j<=L.                    (5.1)
```

Their normalized mantissas in a large dyadic band form an arithmetic
grid of mesh `2^(L+1-h)`, which tends to zero as the height h increases.
Choose one with initial phase in the selected interior subinterval.
Its total carry over these L steps tends to zero as n tends to infinity,
so its crossing word eventually equals the selected rotation word.
This proves the claim for arbitrarily large sources. The modulus
`2^(L+1)`, rather than `2^L`, ensures that the L-th odd successor is still
odd and that the final specified valuation is exactly one.

This is a useful saturation obstruction: **no finite crossing word in
the allowed language forces even one extra even halving.** A proof that
couples phase to valuation therefore has to use information beyond a
fixed finite crossing pattern, such as the evolution and regeneration
of an arithmetic resource across successive blocks.

### 5B. A natural integer resource is consumed, then cheaply regenerated

The resource `R(n)=v2(n+1)` does retain exactly the length of a growing
run. If a=1, then R(n)>=2 and

```
R((3*n+1)/2)=R(n)-1.
```

The integer

```
Q(n)=3^R(n) * ((n+1)/2^R(n))                       (5.2)
```

is invariant along such a run. Equivalently,
`log2(n+1)+alpha*R(n)` is constant. The run of R(n)-1 consecutive a=1
steps ends at the exact odd integer `(2/3)*Q(n)-1`, whose R-value is one.
This gives an arithmetic quotient of the growing run, not merely a
statistical classification.

However, shallow resets replenish arbitrarily much resource. For every
odd k>=3, put

```
s_k=(2^(k+2)-5)/3,      y_k=2^k-1.
```

Both are positive odd integers, and

```
3*s_k+1=4*y_k, a(s_k)=2,
R(s_k)=1, R(y_k)=k.                               (5.3)
```

Only one extra even halving buys k-1 units of new run resource. For
`k=1 or 5 mod6`, s_k is also coprime to three, so the obstruction survives
the necessary congruence condition for an internal odd orbit value.
The example k=5 is exactly **41 -> 31**, inside the orbit of 27.

Consequently no constant c makes
`V_c(n)=log2(n)+c*v2(n+1)` nonincreasing on every sufficiently large odd
Collatz transition. Along large a=1 transitions its increment tends to
`alpha-c`, forcing c>=alpha if the proposed monotonicity held. But on
(5.3) the increment is
`log2(y_k/s_k)+c*(k-1)`, which tends to positive infinity for such c,
since `y_k/s_k -> 3/4`. This explicitly refutes the simplest potential
combining height and the odd-run resource. It does not refute a potential
with further arithmetic coordinates or a nonlocal excursion charge.

There is a stronger obstruction, independent of linearity. Follow (5.3)
by the k-1 growing steps that consume the new resource. This actual
integer macroblock is

```
s_k -> 2^k-1 -> ... -> z_k=2*3^(k-1)-1,
R(s_k)=R(z_k)=1,
z_k/s_k ~ (1/2)*(3/2)^k -> infinity.              (5.4)
```

Every value in the block tends to infinity as k does. Hence **no
potential of the form `V(n)=log2(n)+f(v2(n+1))+g(n)`, with arbitrary
real-valued f and bounded g, is nonincreasing on every sufficiently
large odd transition.** The f terms cancel between these endpoints,
the g difference stays bounded, and the height increase tends to
infinity. This rules out adding any bounded finite-residue or continuous
log-periodic phase correction to this one-resource ansatz. It does not
rule out other unbounded arithmetic coordinates.

For the internal-image subfamily `k=1 or5 mod6`, the entire odd macroblock
consists of units modulo 30. Indeed its values after the reset are
`3^j*2^(k-j)-1` for `0<=j<=k-1`. They are odd; modulo3 they are 1 at j=0
and -1 thereafter; modulo5 the product before subtracting 1 is
`(-1)^j*2^k`, which is 2 or3 since k is odd. The source s_k is also
coprime to2,3,5. Thus passing to the small wheel modulo `2*3*5` does not
delete these resets or their subsequent growing run. For example,
41->31 is the residue transition 11->1 modulo30. This is an exact
compatibility with the user's mod30 seed, not a claim that residue11
is intrinsically exceptional for all prime comparisons.

### 5C. No fixed finite list of prime valuations repairs this potential

The root session found a strengthening; this lane independently audited
the algebra and added exact controls. Let P be **any finite set of
primes**, and define

```
r_P(n)=(v_p(n+1))_(p in P),
L=lcm(6, ord_p(2), ord_p(3): p in P, p>=5).
```

Take arbitrarily large `k=1 mod L` in the previous construction and append
one more odd transition. Since k is odd,

```
3*z_k+1=2*(3^k-1),  v2(3*z_k+1)=2,
z_k -> w_k=(3^k-1)/2.                            (5.5)
```

Then the two endpoint vectors are exactly the same:

```
r_P(s_k)=r_P(w_k)=(v_p(2))_(p in P).             (5.6)
```

For p=2, both `s_k+1` and `w_k+1` have valuation one: the first is
`2*(2^(k+1)-1)/3`, and the second is `(3^k+1)/2`, with k odd.
For p=3, the latter is a unit; also k=1 mod6 gives
`2^(k+2)=8 mod9`, hence `s_k+1=2 mod3`. For p>=5 in P, the multiplicative
order conditions give

```
s_k+1=(2^(k+2)-2)/3 = 2 mod p,
w_k+1=(3^k+1)/2     = 2 mod p.
```

These congruences give valuation zero, not just congruence of two
potentially nonzero valuations. No assumption on independent prime
events is used.

Meanwhile

```
w_k/s_k ~ (3/8)*(3/2)^k -> infinity,              (5.7)
```

and every odd source along the macroblock tends to infinity with k.
Therefore, for every **c>0**, every real-valued function f of the entire
finite resource vector, and every bounded function g, the potential

```
V(n)=c*log(n)+f(r_P(n))+g(n)                       (5.8)
```

cannot be nonincreasing on all sufficiently large odd Collatz
transitions. The f contributions at the endpoints cancel exactly;
the g difference is bounded; the logarithmic height increase diverges.
The qualification c>0 is necessary: c=0 permits a constant potential.

In particular P={2,3,11} gives L=30. Adding the prime 5 gives
P={2,3,5,11} and L=60. Thus the obstruction directly tests both the
user's selected prime triple and the smallest-prime wheel. Every odd
state of the block, including w_k, remains a unit modulo30: the earlier
states were covered in Section 5B, and for odd k, `(3^k-1)/2` is odd and
coprime to3 and5. Adding an arbitrary bounded function of residues,
colour, or logarithmic phase cannot repair (5.8).

This result is deliberately scoped. It concerns valuations of **n+1**
at a fixed finite set of primes, with a bounded remainder and positive
logarithmic height term. It does not rule out other shifts, a growing
prime set, unbounded phase/height sidecars, or a potential defined on
whole excursions rather than individual odd states. Those are explicit
ways to formulate a stronger next construction instead of repeating
the refuted finite-prime ansatz.

## 6. The remaining inequality, stated in the retained coordinates

The new potential controls U_J exactly up to finite carry; the remaining
quantity is E_J. On a hypothetical divergent orbit,

```
E_J-alpha*J = h_0-h_J + C_J+phi_0-phi_J -> -infinity. (6.1)
```

Equivalently, its accumulated valuation excess
`sum a_j-J*log2(3)` tends to minus infinity. This last discrepancy
consequence is inherited from THM-4476; the crossing construction locates
its two constituents and provides a low-complexity language for one of
them. It must not be presented as a new no-divergence theorem.

Three precise next targets survive the controls:

1. **Arithmetic extension of the crossing automaton.** For each L, the
   allowable high-source crossing words are the L+1 rotation words.
   Decorate these with exact source residues and the valuations a_j.
   Seek an endpoint-stable inequality forcing a cumulative excess of
   e over u before a source drops below its starting value. A finite
   residue model must retain its height condition and rule out compatible
   infinite shadows; an unconstrained automaton cycle is not an integer
   counterexample.
2. **A weighted excursion inequality.** Charge a growing return to its
   actual enclosing excursion or a later return, retaining the exact
   affine intercept. Test any proposed charge on all eight growing
   returns in the 27 orbit and on arbitrary hover/drop cylinders. A
   summable or well-founded bank would be substantive; a bank defined by
   unknown future stopping time would merely restate the target.
3. **A carry-sensitive poset inequality.** Study a measure on legal
   rearrangements that preserves a fixed source by explicit compensating
   carry coordinates. Prove the measure and its inequality before using
   balance. The uniform law on unconditioned word permutations already
   fails the fixed-source requirement.

The strongest concrete gain is (3.3): the apparent crossing freedom of a
large integer orbit is restricted to an exactly computable Sturmian
language, while its arithmetic halving depth remains explicitly visible.
Neither arbitrary residue statistics nor source-independent phase
oscillations can stand in for the remaining coupling inequality.

## 7. Reproduction and controls

Run `python 04-computation/experiments/crossroads_crossing_20260926_potential.py`.
The retained stdout is
[crossroads_crossing_20260926_potential.out](crossroads_crossing_20260926_potential.out).
The universe is all odd sources 1<=n<=100000, language lengths 1..16,
with only the explicit minimum-height filter (3.3) inherited for language
tests. There are 50000 exact edge/flux controls and 665811 eligible exact
language controls. The rotation cuts, safe heights and language
membership are computed with integers and rational numbers. A separate
multiplicative identity verifies the logarithmic potential on eight
starts for 80 odd steps each. Printed floating errors are diagnostics,
never pass/fail tests.

The script also constructs all 152 allowed rotation words of lengths
1..16 with actual all-a=1 sources, checks 100 shallow-reset examples,
and checks the integer growing-run invariant for every odd n<=100000.
For the finite-prime extension it checks three increasing k-values for
each of six prime sets, including {2,3,11}, {2,3,5,11}, and every prime
at most19. The computed periods for the two owner-linked sets are30
and60, and the full prime set through19 has period720.

Positive controls include 27, the later first-descent records 703,10087,
35655,626331, and a long-growth Mersenne start. Hostiles are the ordinary
odd fixed point 1, the positive 3n-1 cycle with odd states 5,7, and the
positive 5n+1 cycle with odd states 13,33,83. Their reciprocal sums
diverge, so the finite-carry premise cannot be silently applied. The
forest check is independent of the odd-state phase computation. Every
assertion that gates success uses an explicit exception, so Python -O
does not disable verification.
