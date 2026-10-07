# Arithmetic transfers from the OpenAI mathematics collection

**Status:** the cited OpenAI theorems are accepted inputs for this research
session, as requested by the owner; this is not a theorem audit. **PROVED
relative to those inputs:** the fixed-word correlation and bounded-prime-filter
corollaries below. **PROVED elementary:** the map identities and transfer
boundaries. **FINITE-EXACT:** the companion controls. The proposed source-specific
positivity transfers remain **OPEN**.

## 1. Inheritance and the live questions

The recent route is [cycles, tubes and debts](collatz_cycles_tubes_debt_walk_openai_20261006.md),
especially its two coupled exponent streams and HYP-9217. Its previous OpenAI
reading concerned limit cycles, Ostmann and witnessed choice. The earlier
[seven/twenty-one reading](seven_twentyone_mersenne_openai_math_20261006.md)
already related the pi paper to the **different** number `log_2(3)` and discussed
a joint Dickman statistic. This note adds an exact legal-cylinder application
of ordinary Liouville correlations, not another informal independence claim.

For positivity, the closest mechanisms are
[arithmetic-shell denominators](collatz_shell_denominator_budget_20261005.md),
[the finite moment localizer](collatz_moment_localizer_feasibility_20261005.md), and
[threshold-to-receipt compilation](collatz_weight_threshold_receipts_20261005.md).
The canonical hostile is an everywhere-positive smoothed field with an absent
atom. The least-used sidecars are **arithmetic dependence under target absence**
and **which averaging law a correlation theorem actually controls**.

The concept board is: legal affine cylinders; coupled exponent streams;
arithmetic moment determinants; genus versus class; finite prime filters;
stationary probability versus rooted source weights.

## 2. The supplied list is Euler's idoneal list

The 65 integers supplied by the owner are exactly the idoneal numbers in
[Watkins's Table 1](https://magma.maths.usyd.edu.au/~watkins/papers/IDONcorr.pdf).
Their invariant belongs to the imaginary quadratic **order of discriminant
`-4n`**: there is one proper form class per genus. Equivalently its class group
has exponent at most two. The genus quotient `Cl/Cl^2` therefore loses no
class information in these cases.

This is not the list of fundamental discriminants of class number one.
For example, `n=105` has eight classes; `n=11` is absent even though the
field of fundamental discriminant `-11` has class number one. The relevant
order for `n=11` instead has discriminant `-44`. Watkins explicitly separates
the idoneal and fundamental-discriminant lists.

The script enumerates all primitive reduced positive forms of discriminant
`-4n` for **every `1<=n<=1848`**, with the usual boundary convention. A reduced
form `(a,b,c)` is its own inverse class precisely on the ambiguity boundaries
`b=0`, `|b|=a`, or `a=c`. Exactly the supplied 65 values have only ambiguous
classes. For `n=11` the forms are `(1,0,11),(3,-2,4),(3,2,4)`.

The transferable pattern is genuine: find a quotient whose kernel is proved
trivial on a specified domain. It does not supply such a theorem for a Collatz
residue quotient. The restored coordinate there is still the source and its
guarded word. The finite enumeration does not classify all larger idoneal
numbers, and no such classification is inferred merely from a zero-free region.

## 3. A direct ordinary-correlation theorem on every fixed legal word

[OpenAI family 007, Theorem `q-affine`](https://raw.githubusercontent.com/openai/math/main/preprints/Ordinary-two-point-correlations-of-multiplicative-functions-September-24-2026/build/introduction.tex)
states a power-of-log saving for ordinary Liouville correlations along any
fixed pair of nonproportional positive-slope integer affine forms. The exponent
is absolute; the form-dependent constant need not be effective. Its corrected
Elliott extension assumes a uniformly nonpretentious original multiplicative
factor. Here use the Liouville specialization, denoted `lambda_L` to distinguish
it from our leaf law.

Let a fixed nonempty valuation word `w` have length `l`, valuation sum `A`, and

```
F_w(n)=(P n+B)/Q,       P=3^l, Q=2^A, B>0.
```

Its exact odd source cylinder has a unique odd residue `r` modulo `2Q`, determined
by `P r+B=Q (mod 2Q)`. Set `e=(P r+B)/Q`. For integer `t>=1`,

```
n_t=r+2Q t,       U^l(n_t)=e+2P t,
(2Q)e-(2P)r=2B>0.                                      (1)
```

All prescribed valuations are exact: the final oddness congruence inductively
forces each earlier valuation. Also `n_t>Q`, so no proper prefix can hit ROOT:
such a hit would give `3^j n_t+B_j=2^{A_j}<=Q`, impossible.

**Corollary.** For each fixed `w`,

```
sum_(1<=t<=X) lambda_L(n_t) lambda_L(U^l(n_t))
    = O_w(X/(log X)^c).
```

Consequently equal and unequal total-prime-factor parities each have limiting
relative frequency `1/2` within this cylinder. Equation (1) checks the entire
affine nonproportionality premise. Example: word `12` gives
`n_t=11+16t`, endpoint `13+18t`, determinant `10`.

This decorates exact Collatz edges with an arithmetic statistic. It does not
describe the parity of the **2-adic valuation**, already prescribed by `w`,
nor the long-lag joint law in HYP-9217. Those observables are different functions.
Moreover the word and its coefficients must stay fixed as `X` grows. A growing
word portfolio owes a uniform-in-word estimate or an independently bounded
discarded tail.

The actual rooted weight cannot replace Liouville: exact receipts give
`W(3)=1/6`, `W(5)=1/3`, but `W(15)=1/252 != W(3)W(5)`. A ROOT indicator is not
known to meet the multiplicativity/nonpretentiousness assumptions either;
if Collatz holds, that indicator is the constant one, which is pretentious.

## 4. Catalan, odd zeta and the precise arithmetic gate

[OpenAI family 005](https://raw.githubusercontent.com/openai/math/main/preprints/Catalans-constant-is-irrational-September-24-2026/build/sections/introduction.tex)
proves irrationality of Catalan's constant `G=beta(2)`. Its determinant strategy
cancels an unwanted `zeta(2)` coordinate, obtains rational determinants under
`G in Q`, establishes nonvanishing, and compares finite-prime denominator cost
against real-place decay. This is distinct from Catalan's exponential equation
`9=8+1` used in our word-collision work. It is also not an irrationality theorem
for the individual real value `zeta(5)`.

The repo's [Apéry framework](../reference/apery-style-irrationality-framework.md)
and [odd-zeta inheritance](collatz_depth_layers_zeta5_preprints_20260927.md)
already separate nonvanishing, decay and denominator clearing. For our actual
localized moments `H_k=sum_j p_j h_m(j)^k`, infinitely many explicit rooted leaves
give infinite off-target support. Their finite Hankel Gram determinants are
therefore strictly positive even under the assumption `p_m=0`.

That identifies a focused new transfer hypothesis:

> Construct a determinant `Delta_N` of arithmetic combinations of the `H_k`
> such that `p_m=0` forces a controlled rational value; prove it nonzero, bound
> its denominator after all cancellations, and make its real decay beat that
> bound.

The missing premise is **target absence implies arithmetic collapse**. Our
current moments have no known finite period span that supplies this implication.
Integer polynomial coefficients alone do not make their expectations rational.
The known denominator `prod_(d<=N)(2^d-1)^2` and the small signed shell readout
are not a substitute. The conditional small-nonzero-integer gate in the shell
denominator note is the exact interface to pursue.

## 5. Pi's approximation theorem and zero-free regions

[OpenAI family 017](https://raw.githubusercontent.com/openai/math/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026/build/main.tex)
gives `mu(pi)=2`, with an ineffective eventual threshold for each exponent
greater than two. The Collatz near-cycle cost uses `alpha=log_2(3)`, not `pi`.
The earlier repo proposal to transfer the interpolation method to this `alpha`
remains a research hypothesis.

There is an exact payoff **if** such an analogue is established. For
`A_min(l)=ceil(l alpha)`, suppose eventually
`|alpha-a/l|>=l^(-2-epsilon)` for every integer `a`. Then

```
e_l=A_min(l)-l alpha >= l^(-1-epsilon),
q_max(l)=ceil(l/e_l) <= 1+l^(2+epsilon).                 (2)
```

An effective version would bound the expense of a word of known length. It
would still not bound that length from the supplied source. Equation (2) is
a transfer implication, not a claim that the pi theorem applies to `alpha`.

[OpenAI family 003](https://raw.githubusercontent.com/openai/math/main/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/build/paper.tex)
gives the common zero-free half-plane `Re(s)>7/8` for Dirichlet L-functions and
finite-order Hecke L-functions over `Q(sqrt(-3))`, with principal poles allowed.
This can feed arithmetic estimates built from those exact Euler products.
Our source-weight series does not inherit an Euler product: the weight
multiplicativity control above already shows why a direct identification fails.

For idoneal classification the relevant proposed route is a quantitative
class-number lower bound, compared with the small exponent-two class-group
size. It still needs its constants and a finite exclusion calculation. For
Collatz, a useful route would specify a Dirichlet/Hecke expansion of the
required residue statistic and control the coefficients uniformly in the
growing modulus. Neither pointwise atom positivity nor a debt-walk recurrence
follows from the zero-free region alone.

## 6. Serre positivity gives a precise, demanding interface

[OpenAI family 193, main theorem](https://raw.githubusercontent.com/openai/math/main/preprints/Positivity-of-Serres-Intersection-Multiplicity-September-23-2026/build/sections/introduction.tex)
asserts `chi^R(M,N)>0` for nonzero finitely generated modules over a regular
local ring, when the tensor product has finite length and their dimensions
are complementary. Here `chi` is the alternating sum of Tor lengths. Its
strict sign belongs to that algebraic intersection invariant.

There is an exact toy realization of our word maps. Every nonempty word has
`P!=Q` and a rational fixed point `a=B/(Q-P)`. In local coordinates at `(a,a)`,
use `R=Q[[u,v]]`, `M=R/(Qv-Pu)`, `N=R/(v-u)`. Their intersection is transverse;
the tensor product is `Q`, higher Tor vanishes, and `chi=1`. This works for
word `1` with fixed point `-1`, word `12` with fixed point `-5`, and word `14`
with fixed point `5/23`, despite their different real directions and absence
of a positive-integer cycle at those anchors.

Thus a positivity transfer needs a much stronger construction: modules
`M_m,N_m` whose nonzero finite proper intersection is independently proved,
together with an identity or lower comparison from `chi` to **the designated
atom** `p_m`. Encoding only the rational affine graph gives the constant one
and forgets integer source legality, ROOT boundary and the infinite orbit.
This says exactly which extra representation would be useful.

## 7. Three additional catalogue selections

**Priority 1: family 021, Jacobsthal.**
[The main theorem](https://raw.githubusercontent.com/openai/math/main/preprints/A-quadratic-bound-for-Jacobsthals-function-September-25-2026/build/sections/introduction.tex)
uniformly bounds the longest gap avoiding `k` prescribed prime divisors by
`C k^2/(log log(3k))^2`; the constant is not supplied numerically.

For a legal source family `n_t=r+a t`, impose `p` not dividing `n_t` for a
fixed finite prime set with `p` not dividing `a`. Each prime forbids precisely
`t=-r a^(-1) (mod p)`. CRT combines them into one shift `b`, so the condition
is exactly `gcd(t-b,prod p)=1`. Jacobsthal therefore gives an asymptotically
short search for a filtered **new family member**, while retaining the whole
word guard. If `p` divides `a`, the filter is either automatic or impossible;
that case must be removed first. Multiple independent forbidden classes at
one prime require a different theorem. This improves parameter selection,
not coverage of an already supplied integer.

**Priority 2: family 145, multiple mixing.**
[The main theorem](https://raw.githubusercontent.com/openai/math/main/preprints/Rokhlins-multiple-mixing-problem-for-one-transformation-September-23-2026/build/sections/01-introduction.tex)
upgrades mixing of an invertible probability-preserving transformation to
all finite orders. It could upgrade a suitably constructed stationary debt
joining after two-mixing is proved. It does not construct that joining or
its mixing property. Marginal geometric streams alone leave the joint law open.

A useful elementary boundary is stronger than a generic warning: on positive
odds with ROOT made absorbing, no invariant probability can have positive
mass at 5, because `U(1)=U(5)=1` forces
`mu({1})>=mu({1})+mu({5})`. Our positive rooted weights are subinvariant,
not invariant. A stationary transfer must use a different state space and
retain how arithmetic sources are selected from it. The two-sided iid
valuation shift is already mixing of all orders, so the genuinely new target
is an arithmetic-conditioned, non-product joining.

**Priority 3: family 022, weak inhomogeneous Duffin--Schaeffer.**
[The stated theorem](https://raw.githubusercontent.com/openai/math/main/preprints/The-weak-inhomogeneous-Duffin-Schaeffer-conjecture-September-25-2026/build/main.tex)
turns divergence of `sum phi(q) psi(q)/q` into infinitely many
`||q x-gamma||<psi(q)` for almost every `x`, for any fixed shift `gamma`.
This gives a metric comparison model for near-resonance visits. It does not
place the specific `x=log_2(3)` outside the exceptional set, nor supply a
one-sided prescribed carry or a source-specific stopping time. It ranks below
Jacobsthal for a literal current compiler improvement.

## 8. Reproduction and the next useful obligations

[Script](../../04-computation/experiments/collatz_openai_arithmetic_bridges_20261007.py)
and [output](collatz_openai_arithmetic_bridges_20261007.out):

```
python -B 04-computation/experiments/collatz_openai_arithmetic_bridges_20261007.py
python -B -O 04-computation/experiments/collatz_openai_arithmetic_bridges_20261007.py
```

The explicit universe is all 1848 candidate idoneal integers, all 120 words
of lengths 1--4 over valuations `{1,2,3}` with five actual sources per word,
2310 fixed-filter parameter controls, and the stated exact weight/intersection
and type controls. It contains no numerical test of an accepted theorem's proof.

The immediate positive tasks are now precise: use the fixed-word correlation
as an arithmetic decoration of the debt cylinder census; implement the finite
prime-filter parameter interface when a controller asks for one; and search
for a target-absence arithmetic relation among moments. For long-lag exponent
correlations, a possible route is an arithmetic stationary joining with
controlled conditioning and recurrence. The companion
[paired-exponent proof](collatz_paired_exponent_information_20261007.md)
instead works directly with nonstationary joint cylinders: it proves finite
observation forgetting and quantitative remote decorrelation. Recurrence
under debt conditioning remains the missing conclusion in either approach.
