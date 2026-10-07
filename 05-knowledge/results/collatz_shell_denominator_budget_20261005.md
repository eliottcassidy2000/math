# The denominator bill behind the arithmetic-shell selectors

2026-10-05. **PROVED:** the exact least coefficient denominator, its prime
clock, invariance under integral basis changes, an unavoidable normalization
bill even with additional degree, and a damped integer-polynomial repair.
**PROVED:** the actual localized-moment generating function has infinitely
many uncancelled poles and is not D-finite. **CONDITIONAL:** the explicit
readout-lattice hypothesis in section 5 would force the selected atom positive.
That hypothesis has not been established. **FINITE-EXACT:** the declared
arithmetic controls. No odd-zeta irrationality or universal Collatz theorem
is claimed by this package.

Script: `04-computation/experiments/collatz_shell_denominator_budget_20261005.py`.
Output: `05-knowledge/results/collatz_shell_denominator_budget_20261005.out`.

## 1. Inheritance: three separate proof gates

The warm-start reference
`05-knowledge/reference/apery-style-irrationality-framework.md` isolates
integrality, nonvanishing, and decay **after** denominator clearing. The
external p-adic-zeta source audit
`05-knowledge/reference/p-adic-zeta-irrationality-source-audit-20260825.md`
explicitly separates a verified numerical margin from the additional
geometric and arithmetic hypotheses that would make it a theorem.

For the real odd-zeta mechanism, Brown's section 1.2 explicitly separates
nonzero linear forms in a fixed collection of values, analytic decay, and
coefficient-denominator control. Only a favorable combined bound proves an
irrationality statement; small forms alone do not:
[Brown, *Irrationality proofs for zeta values, moduli spaces and dinner parties*](https://www.ihes.fr/~brown/IrratModuliMotivesv8.pdf).

For a precise primary example, Lai, Sprang and Zudilin construct nonzero
rational approximations with a quantitative height-versus-2-adic-error
inequality for **the 2-adic value** zeta_2(5). Their statement is not about
the real value zeta(5):
[author preprint, arXiv:2505.05005](https://arxiv.org/abs/2505.05005).
We use only this distinction of proof obligations; none of their zeta
arithmetic is transported to the Collatz moment law.

The closest internal algebra is:

* THM-4091,
  `01-canon/theorems/THM-4091-integral-coordinate-change-lcm-depth-boundary.md`:
  coordinate changes must carry the denominator lattice, not merely a
  formal coefficient formula. The theorem supplies no irrationality.
* THM-4056,
  `01-canon/theorems/THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock.md`:
  exact prime-power addresses retain valuation depth and lose the named
  target unless it is separately retained.
* `collatz_posets_dags_zeta5_20260927.md`, section 3: orbit linear forms have
  a valuation-feedback rule; their coefficients cannot be chosen as freely
  as in a separate approximation construction. Only this typed comparison
  is inherited here, not every historical proposed consequence.
* `collatz_refinement_energy_dual_20261005.md`: the normalized shell
  minorants have direct-moment coefficient norm below 512. This controls
  sensitivity to moment errors, not arithmetic bit height or integrality
  of a measured linear form.

The board is **normalization / coefficient lattice / integer shell support /
readout lattice / nonvanishing / analytic carrier**. The decisive hostile
below has integer polynomial coefficients and a nonzero readout smaller
than one, because the readout is a rational noninteger.

## 2. The least denominator is exact

For `d>=1` let

```text
x_d = 4*2^d/(2^d+1)^2,
Q_N(h) = product_(d=1)^N (h-x_d)/(1-x_d),
D_N = product_(d=1)^N (2^d-1)^2,
A_N(h) = product_(d=1)^N [(2^d+1)^2 h-2^(d+2)].             (1)
```

Then `Q_N=A_N/D_N` and `A_N(1)=D_N`. The denominator identities hold for
every positive `N`; the atom-minorant sign requires odd `N`.

The exact constant coefficient is

```text
A_N(0)=(-1)^N 2^(N(N+5)/2).                               (2)
```

Since `D_N` is odd, the constant coefficient of `Q_N` is already reduced.
Thus **the least common denominator of all coefficients is exactly D_N**.
There is no global cancellation hidden in the product. Also `A_N` is
primitive: its constant coefficient is a power of two, while its leading
coefficient is odd, so their gcd, and hence the content, is one.

The q-product bound proved in the preceding package gives

```text
(81/1024) 2^(N(N+1)) <= D_N < 2^(N(N+1)).                  (3)
```

Indeed `D_N=2^(N(N+1))*product_(d<=N)(1-2^-d)^2`, and that product without
the square is at least `9/32`. Consequently the denominator has quadratic
bit height. The integer numerator does too: its leading coefficient is
larger than `2^(N(N+1))`, whereas the inherited coefficient-norm bound
gives `sum |[h^k]A_N|<512D_N` for odd `N`.

This is compatible with a small norm for the **rational** coefficients.
Magnitude, sensitivity, denominator, and representation size are distinct.

### Extra degree cannot remove the normalization bill

Suppose a rational polynomial `P` satisfies

```text
P(1)=1,             P(x_d)=0 for d=1,...,N.                 (4)
```

If `E` is a positive integer with `EP` integral coefficientwise, then

```text
D_N divides E.                                             (5)
```

Proof: each primitive linear factor in `A_N` divides `EP` in `Z[h]`.
For completeness, if `sh-r` has coprime integer `s,r` and divides an
integer polynomial rationally, evaluation at `r/s` shows its leading
coefficient is divisible by `s`. Subtract the corresponding integral
multiple of `h^(degree-1)(sh-r)` and continue by degree induction. Distinct
roots let us repeat this argument for every factor. Thus `EP=A_NS` for
some `S in Z[h]`, and evaluation at 1 gives `E=D_NS(1)`.

Therefore added degree, a different rational representation, or extra
factors can improve approximation, but cannot lower this least necessary
normalization denominator while retaining all the specified zeros.

An integer unimodular change of polynomial basis preserves the coefficient
lattice exactly: `Ec` is integral iff `EUc` is integral when both `U` and
`U^(-1)` are integral. A rational basis may move the denominator into the
observable. Writing the answer as “one times Q_N” does not make its
expectation an integer. Scaling by `D_N` makes the target value `D_N`, so
dividing by that same value is required to recover an atom floor.

## 3. Cyclotomic and prime-power accounting

Regrouping `2^d-1` by cyclotomic factors gives the exact identity

```text
D_N = product_(r=1)^N Phi_r(2)^(2 floor(N/r)).              (6)
```

These numerical factors are not pairwise coprime: `Phi_2(2)=Phi_6(2)=3`.
One must retain overlapping prime valuations rather than count the displayed
factors as independent primes.

Let `ell` be an odd prime, `e=ord_ell(2)`,
`a=v_ell(2^e-1)`, and `K=floor(N/e)`. Then

```text
v_ell(D_N)=2[K a+v_ell(K!)].                               (7)
```

Only multiples `d=et` contribute. The elementary odd-prime binomial
identity gives `v_ell((2^e)^t-1)=a+v_ell(t)`: powers coprime to `ell`
preserve the initial valuation, and raising to the `ell`th power adds one.
Summing over `t=1,...,K` proves (7). This is a precise valuation clock,
not a claim that the Collatz dynamics has this period.

## 4. Undamped clearing fails; damping is a real repair

For odd `N`, the inherited uniform normalized error is

```text
epsilon_N=product_(d=1)^N x_d/(1-x_d).
```

Its denominator-cleared value is exactly

```text
D_N epsilon_N=2^(N(N+5)/2).                                (8)
```

Thus the undamped small-error estimate cannot give a small integer after
coefficient clearing. The clearing cost grows faster than that error shrinks.

This does **not** obstruct every modified selector. Define

```text
I_(N,k)(h)=h^k A_N(h),       k>=0.                         (9)
```

It has integer coefficients and target value `D_N`. On the first `N`
off-target shells it is zero, and at every remaining shell it is negative.
Since those shells are at most `x_(N+1)`,

```text
0 < -I_(N,k)(x_j)
  <= B_(N,k):=2^(N(N+5)/2) x_(N+1)^k,     j>N.             (10)
```

For odd `N>=3`, choose `k=(N+9)/2`. Because
`x_(N+1)<2^(1-N)`,

```text
B_(N,k) < 2^((9-3N)/2) <= 1,                             (11)
```

which tends to zero with `N`. The total degree is `(3N+9)/2`.
The normalized minorant `I_(N,k)/D_N=h^k Q_N` still has coefficient norm
below 512: multiplication by `h^k` only shifts its coefficient list.
Its least common coefficient denominator remains exactly `D_N`.

This successfully repairs the **coefficient-cleared smallness** step.
It does not repair the separate readout-integrality step.

## 5. The exact missing Apéry gate

Even an integer polynomial can have a nonzero expectation of magnitude less
than one. For a point mass at the valid off-target shell `x_2=16/25`,

```text
E[h^2(9h-8)] = -14336/15625,
0 < |E[h^2(9h-8)]| < 1.                                  (12)
```

The polynomial is integral coefficientwise; the readout is not an integer.
Clearing its actual denominator gives `-14336`, and smallness disappears.
Conversely a law supported on an annihilated shell gives an integral
readout equal to zero, so nonvanishing cannot be omitted either.

There is an all-length support-only obstruction. A rational function regular
at zero that takes values in a fixed lattice `(1/B)Z` at every `x_d` must be
constant. Its values converge to its value at zero. A convergent sequence
in this discrete lattice is eventually constant, and a rational function
equal to that constant at infinitely many distinct points is identically
constant. In particular, no nonconstant polynomial becomes uniformly
integer-valued on all shells after one fixed clearing factor. This also
rules out such an integrality claim derived solely from support for **all**
probability laws, since point masses test it. It does not rule out extra
identities for one specific actual law.

For the actual Collatz leaf law, define the moment readout

```text
J_N=E[I_(N,(N+9)/2)(h_m(j))].                             (13)
```

The known rooted leaves in section 6 supply positive mass at infinitely
many off-target shells. Therefore, **under the assumption p_m=0**, (10)
gives

```text
0 < -J_N <= B_(N,(N+9)/2) < 2^((9-3N)/2).                (14)
```

Here nonvanishing is paid independently of the selected source. A precise
sufficient arithmetic input would be positive integers `E_N` such that,
under that same assumption,

```text
E_N J_N is an integer,
E_N 2^((9-3N)/2) -> 0.                                   (15)
```

Then (14)--(15) would force a nonzero integer to have magnitude below one,
a contradiction, proving `p_m>0`. For example, the stronger hypothetical
bound `E_N<=2^N` would suffice. **No such readout lattice has been proved.**
Integer coefficients in (13), rationality of each separately known atom,
or the assumption `p_m=0` do not imply (15) for an infinite moment sum.

This is the lawful Apéry-style connection: it identifies the exact new
arithmetic lemma that would be needed, rather than treating a successful
approximation as the missing lemma. The forms use an increasing collection
of localized moments; they are not already known to reduce to a fixed list
of rational constants and the selected atom.

## 6. The raw actual moment sequence is not P-recursive

There is a second carrier boundary, independent of any conjecture. For a
fixed target `m`, define near zero

```text
F_m(z)=sum_(k>=0) H_k z^k
      =sum_(j>=0) p_j/(1-h_m(j)z).                        (16)
```

The interchange is absolutely justified for `|z|<1`. Group equal shell
values, writing their nonnegative total masses as `w_d`. Since `x_d->0`
and `sum w_d=1`, the series `sum_d w_d/(1-x_d z)` converges locally
uniformly on every compact subset of the plane avoiding its possible poles:
for all sufficiently large `d`, each denominator has magnitude at least
`1/2`, uniformly on that compact set. Thus it gives a meromorphic extension,
with a simple pole at `z=1/x_d` whenever `w_d>0`. Its residue is
`-w_d/x_d`, which cannot be cancelled by other shells.

Infinitely many of those masses are unconditionally positive. Namely,

```text
n_t=(4^(3t)-1)/3,       t>=1,
3n_t+1=2^(6t),          n_t=6j_t+3,
p_(j_t)=lambda(n_t)=W(n_t)=2/[(3t)(3t+1)]>0.              (17)
```

These are the inherited one-edge ROOT leaves. Their indices tend to infinity,
so for any fixed `m` they give infinitely many distinct distances and poles.
They do not presuppose that the selected source `6m+3` is rooted.

Suppose `F_m` satisfied a nonzero homogeneous linear differential equation
with polynomial coefficients, with highest derivative order `r`. At any
simple pole where its leading polynomial is nonzero, the term containing
`F_m^(r)` has pole order `r+1`; all lower-derivative terms have smaller pole
order and cannot cancel it. Hence the leading polynomial must vanish at
every one of the infinitely many poles, a contradiction. This includes
order zero. Therefore `F_m` is **not D-finite**.

Consequently `(H_k)` satisfies no nonzero finite-order recurrence with
polynomial coefficients in `k`. Indeed multiplying such a recurrence by
powers of `z` and summing replaces polynomial factors in `k` by powers of
`z d/dz`; finitely many initial terms contribute only a polynomial, which
can be annihilated by further differentiation. This would give the forbidden
homogeneous differential equation.

This excludes an Apéry-type P-recursive search for these **raw localized
moments**. It does not exclude a different carrier, a nonlinear recurrence,
an auxiliary-parameter construction, or special relations among selected
linear combinations. The finite generating polynomials `A_N` themselves
remain explicitly computable; they are not the function in (16).

### A different, entire carrier is available

The independently proved reciprocal product in
`collatz_reciprocal_selector_20261005.md` is

```text
R(z)=product_(d>=1) (1-x_d z)/(1-x_d),
R(1)=1,
B(z)=(1-z)R(z).                                          (18)
```

It has exactly the simple zeros needed to cancel the possible poles of
`F_m`. Therefore

```text
E_m(z)=B(z)F_m(z)                                         (19)
```

extends to an entire function, with

```text
E_m(1)=p_m,
E_m(1/x_d)=-w_d B'(1/x_d)/x_d.                            (20)
```

These formulas are limits of the local simple-pole expansions. They remain
valid when the corresponding mass is zero, in which case `F_m` is regular
and the product vanishes there. Since `B'(1/x_d)` is nonzero, the second
formula also recovers every aggregated off-target shell mass from this
entire carrier. It loses only the original left/right distinction within
a shell, which `H_k` had already lost.

Thus the raw-moment non-D-finite result is a carrier boundary, not a ban on
changing carriers. The entire function (19) is not thereby proved D-finite,
and its value at 1 is not thereby proved positive. Proving an independent
lower bound or a useful arithmetic recurrence for this new carrier remains
a separate obligation.

## 7. Validation and transfer boundaries

| Source | Map or operation | Retained fact | Missing coordinate |
|---|---|---|---|
| Normalized shell polynomial | `Q_N -> A_N=D_N Q_N` | Integer coefficients and all zeros | Target scale changes from 1 to D_N; readout need not be integral |
| Product denominator | Cyclotomic and prime-valuation grouping | Exact total clearing cost | Numerical cyclotomic factors overlap in primes |
| Integral polynomial | Multiply by `h^k` | Integer coefficients, zeros, strict tail sign | Higher moment orders; actual readout lattice still unproved |
| Known rooted leaf family | Positive shell masses in (16) | Infinite noncancelling poles and nonvanishing in (14) | Does not establish positivity of an arbitrary selected atom |

The standard-library script independently expands the integer numerator
and compares it with the prior rational-factor implementation. Its finite
universe includes:

* `N=1..25` for least denominators, contents, roots, bit-height bounds, and
  cyclotomic regrouping;
* `N=1..50` and every odd prime below 100 for (7);
* 25 normalized rational quadratic extensions at each `N=1,3,5,7,9`, plus
  integer unimodular shears;
* the damped tail at the first 20 remaining shells for odd `N=3..25`;
* 48 compressed explicit rooted-leaf/pole-address controls, the small
  noninteger hostile, and nine malformed-input controls.

No moment oracle is read. The infinite statements follow from the proofs,
not a census of the finite controls. Reproduce from the repository root:

```text
python -B 04-computation/experiments/collatz_shell_denominator_budget_20261005.py
python -B -O 04-computation/experiments/collatz_shell_denominator_budget_20261005.py
```

Both modes perform **2,465 exact checks** and agree with the saved output.
