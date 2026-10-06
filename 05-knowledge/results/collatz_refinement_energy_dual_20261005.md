# Arithmetic-shell refinement of source-atom moment certificates

**PROVED:** the exact integer-distance support supplies an increasing family of
signed atom minorants. Their power-basis coefficient norms are uniformly below
512, and their uniform tail errors decay quadratically in the exponent of the
degree. **FINITE-EXACT:** the same four moments give a positive floor under this
support constraint while admitting a zero-target completion on the relaxed
interval support. **CONDITIONAL:** a floor for an actual Collatz source requires
independently justified simultaneous moment intervals. No such new actual-source
measurement or universal positivity theorem is claimed here.

Script: `04-computation/experiments/collatz_refinement_energy_dual_20261005.py`.
Output: `05-knowledge/results/collatz_refinement_energy_dual_20261005.out`.

## 1. Inheritance and the missing coordinate

The closest mechanisms are:

* `collatz_localized_resolvent_floor_20261005.md`: for
  `h_m(j)=4t/(1+t)^2`, `t=2^(m-j)`, the signed selector
  `(9h-8)h^d` has coefficient norm 17 in the direct localized moments and
  error at most `8(16/25)^d`.
* `collatz_moment_localizer_feasibility_20261005.md`: optimization of
  `(h-8/9)P(h)^2/(1-8/9)`, including the unique actual-law optimizer, keeps
  the off-target interval `[0,8/9]`. Its finite PSD tests do not exclude
  continuum laws at points absent from the actual arithmetic support.
* SF4 in `collatz_source_refinement_floor_20261005.md`: bounded total energy
  of **exact evaluation representers** yields an atom floor. Its corrected
  hostile shows that bounded energy of normalized source-centered densities
  does not. Evaluation, approximation, and normalization must remain typed.
* THM-3674,
  `01-canon/theorems/THM-3674-sharp-successor-variance-drift-and-target-energy-tariff.md`:
  a lawful support zero strengthens a generic energy inequality. The transfer
  here is this method of retaining the lawful support, not its Fourier tensor
  or LRC conclusion.
* THM-4263,
  `01-canon/theorems/THM-4263-moving-multigraph-filtered-jet-and-finite-factor-density-transport.md`:
  fibre weights and conditional hazards must accompany refinement. No
  independence hypothesis is silently added here. THM-2370,
  `01-canon/theorems/THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile.md`,
  supplies the related warning that squared energy can lose signed target
  information.

The live comparison board is: **actual discrete support / interval relaxation /
signed moment floor / coefficient conditioning / exact evaluation energy**.
The new information is the first item. The experiment adds no moment, ROOT
certificate, or orbit observation to manufacture its gain.

Let `(p_j)_(j>=0)` be a probability law and fix a target index `m>=0`. Define

```text
x_d = 4*2^d/(2^d+1)^2,       d=0,1,2,...,
h_m(j)=x_|j-m|,
H_k=sum_j p_j h_m(j)^k.                                    (1)
```

Thus `x_0=1`, and the off-target support is the strictly decreasing sequence

```text
x_1=8/9, x_2=16/25, x_3=32/81, ... -> 0.                  (2)
```

There is no off-target shell in any open gap `(x_(d+1),x_d)`. Both sides of
the source may occupy the same shell; their masses are aggregated by `H_k`.
Only `j=m` gives value 1, so this aggregation preserves the selected atom.
For the Collatz application, `p_j=lambda(6j+3)` is the inherited normalized
leaf law. None of the proofs below assumes this selected atom is positive.

## 2. A support-aware refinement theorem

For any positive odd integer `N`, set

```text
Q_N(h) = product_(d=1)^N (h-x_d)/(1-x_d),
epsilon_N = product_(d=1)^N x_d/(1-x_d).                    (3)
```

On the support in (1),

```text
Q_N(h_m(j)) <= 1_(j=m) <= Q_N(h_m(j))+epsilon_N.             (4)
```

The target value is 1. The first `N` off-target shells are zeros. At every
remaining shell all `N` numerator factors are negative; `N` is odd, so
`Q_N<0`, with magnitude at most `epsilon_N`. This proves (4). The bound is
sharp as a supremum on the support closure, since `Q_N(0)=-epsilon_N`.
Zero is a limit shell, not an extra finite source atom.

Write `Q_N(h)=sum_(k=0)^N q_(N,k)h^k`. Then

```text
A_N = sum_(k=0)^N q_(N,k) H_k,
A_N <= p_m <= A_N+epsilon_N.                               (5)
```

If an independent arithmetic or measure estimate gives
`sum_(|j-m|>N) p_j <= t_N`, replace the error in (5) by
`epsilon_N*t_N`. The program accepts this as a separate premise; it does
not infer the tail cap from the source index or a finite certified bank.

For odd `N>=3`, there is genuine pointwise monotonicity:

```text
Q_(N+2)(h_m(j)) >= Q_N(h_m(j)),
A_(N+2) >= A_N,             lim_(N odd -> infinity) A_N=p_m. (6)
```

The two newly removed shells become zero. Below them,
`Q_(N+2)=Q_N R`, where `R` is the product of the two new normalized factors.
It lies between zero and one: both new `x_d` are below `1/2`, and each
factor's absolute value is at most `x_d/(1-x_d)<1`. Multiplying a negative
number by `R` increases it. The old zeros and target are unchanged.
Uniform convergence on the actual support follows from the next bound.
The qualification `N>=3` matters: `Q_3(x_10)<Q_1(x_10)`.

This is a discrete positive-cone extension of the old localizer. In fact,

```text
Q_(2s+1)(h) = (h-x_1)/(1-x_1)
  * product_(i=1)^s [(h-x_(2i))(h-x_(2i+1))]
                    /[(1-x_(2i))(1-x_(2i+1))].             (7)
```

Each paired gap factor is nonnegative at every possible shell, although it
is negative on the excluded interval between its roots. Replacing it by a
square would discard this lawful arithmetic distinction.

## 3. Refinement without direct-moment coefficient blowup

The exact error is

```text
epsilon_N = 2^(-N(N-3)/2)
            * product_(d=1)^N (1-2^-d)^(-2)
          <= (1024/81)*2^(-N(N-3)/2).                       (8)
```

For `N>=2`, the first two factors in the reciprocal product have product
`(1/2)(3/4)=3/8`. The remaining product is at least
`1-sum_(d>=3)2^-d=3/4`, using the finite product inequality and then its
limit. Thus the full product is at least `9/32`. The case `N=1` also meets
the bound. In particular the error decreases much faster than a fixed
geometric sequence in `N`.

All roots of `Q_N` are positive, so its coefficients alternate in sign.
Its exact coefficient norm is

```text
C_N = sum_k |q_(N,k)|
    = product_(d=1)^N (1+x_d)/(1-x_d) < 512.                (9)
```

Here is an explicit bound, not a numerical product estimate. The first
five factors give

```text
C_5 = 33835804361/95355225.
```

For `d>=6`, `x_d<=4*2^-d<=1/16`. Write the next norm factor as
`1+a_d`, where

```text
a_d=2x_d/(1-x_d) <= (128/15)*2^-d,
sum_(d>=6) a_d <= 4/15.
```

For nonnegative `a_d` with sum at most `s<1`, expanding products and bounding
each elementary symmetric sum by `s^k` gives `product(1+a_d)<=1/(1-s)`.
Therefore

```text
C_N <= (15/11)C_5 = 33835804361/69927165 < 512.              (10)
```

For `N<5`, simply use `C_N<=C_5`.

Suppose **simultaneously valid** moment intervals `l_k<=H_k<=u_k` have been
established independently. The exact rational floor is

```text
L_N = sum_(q_k>=0) q_k l_k + sum_(q_k<0) q_k u_k <= p_m.    (11)
```

Its upper counterpart plus `epsilon_N*t_N` is an atom upper bound. If all
interval widths are at most `omega`, then `A_N-L_N<=C_N*omega<512omega`.
Repeated packets can be intersected on shared moments; the running maximum
of all lower bounds and minimum of all upper bounds are nested. This
retains a valid older witness when a new noisy readout is weaker.

As a conditional accuracy statement, if `p_m>=eta>0`, choose odd `N>=3` with
`epsilon_N<=eta/4` and truthful widths at most `eta/2048`. Then (11) is at
least `eta/2`. The sufficient degree is `O(sqrt(log(1/eta)))`, with the
explicit integer inequality (8) in place of any asymptotic assumption.
This is a cost comparison, not an independently established `eta`.

**Oracle boundary.** The bound 512 applies to the direct localized moments
`H_0,...,H_N`. If these are reconstructed from the primitive resolvent
moments in `collatz_localized_resolvent_floor_20261005.md`, the inherited
expansion of `H_k` has coefficient norm `8^k`; a conservative combined
bound is then `512*8^N`. High-precision or independent direct moment access
has not been supplied for free. Passing the finite checker does not prove
the packet arises from the actual law.

## 4. A strict same-measurement gain and a decisive ghost

Consider the rational synthetic law

```text
mass 1/1000 at 1;
mass 333/1000 at each of x_1,x_2,x_3.                      (12)
```

It is a lawful integer-distance distribution: at `m=0`, use atom indices
`0,1,2,3`. Its moments through order three are

```text
H_0=1,
H_1=144377/225000,
H_2=206161417/455625000,
H_3=316191587057/922640625000.                              (13)
```

The new cubic gives `A_3=1/1000`. But the old degree-one interval localizer
`C[i,j]=(8/9)H_(i+j)-H_(i+j+1)`, `i,j=0,1`, is positive definite, and its
normalized optimum has floor

```text
-811338319/7894961000 < 0.                                 (14)
```

The limitation is exact, not a poor optimizer: the same four moments also
come from this positive law with **zero target mass**:

| Node | Mass |
|---|---:|
| `32/81` | `111782799/337280000` |
| `1/2` | `9/2125` |
| `16/25` | `1305719/3968000` |
| `8/9` | `107289/320000` |

The node `1/2` lies strictly between `x_3` and `x_2`, and therefore cannot
occur in (1). Indeed `Q_3(1/2)>0`; applying the arithmetic minorant to this
continuum ghost would be invalid. Thus an exact integer-support constraint
strictly improves the floor with the **same** `H_0,...,H_3`, independently
of any ROOT certificate for the target. This example is synthetic, not an
asserted measurement of an actual Collatz source.

## 5. What the energy does and does not say

Put `D_N=Q_(N+2)-Q_N` for odd `N>=3`. On the actual support, `D_N>=0`, and

```text
sum_(N odd>=3) D_N = 1_(j=m)-Q_3(h_m(j)).                  (15)
```

Consequently, for every law, including those with `p_m=0`,

```text
sum_N E[D_N^2]
 <= E[(sum_N D_N)^2] <= epsilon_3^2.                       (16)
```

There is no energy-to-positivity inference here. A point mass at shell four
has target mass zero, `Q_3(x_4)<0`, `Q_5(x_4)=0`, and all later values zero;
its total increment energy is the finite positive number `Q_3(x_4)^2`.
These are bounded approximants to the atom indicator, not SF4's normalized
evaluation representers of norm squared `1/c_N`. The signed expectation in
(11), rather than finiteness of (16), is what can pay a floor.

For the actual Collatz law, a separately justified positive floor for
`p_m=lambda(6m+3)=W(6m+3)` feeds the inherited B-step deadline or finite
superlevel receipt compiler. The actual/source-specific positivity problem
is not resolved by the existence or uniqueness of the polynomials `Q_N`.
A canonical first successful odd `N` can select one certificate after
truthful readouts make it positive; this rule supplies no existence proof.

## 6. Typed connection and finite validation

| Source | Target and map | Preserved predicate | Loss and required sidecar |
|---|---|---|---|
| Integer source indices `j` with law `p` | `h_m(j)=x_|j-m|` | Selected atom is exactly the value 1 | Left/right shell masses merge; retain fixed target `m` and integer-distance support |
| Direct moment packet through `N` | `sum q_(N,k)H_k` | Pointwise minorant gives an atom lower bound | Requires simultaneous truthful intervals; finite algebra does not authenticate them |
| Odd-degree arithmetic refinement | Multiply by two adjacent-shell factors | Minorants increase for `N>=3`; direct coefficient norm stays bounded | Continuum gap points would invalidate the sign; degree-one start is not monotone |
| Squared increment energies | Sum in (16) | Finite nonnegative refinement cost | Loses the signed readout; zero target remains possible |

The standard-library script has no orbit, ROOT receipt, weight-bank, or
external-oracle inputs. Its explicit finite universe is:

* source indices `0..6`, atom indices `0..40`, odd degrees `1..17`, comparing
  coefficient evaluation with an independently formed product;
* norm/error controls through odd degree 41;
* all 70 denominator-four laws on the target and first four shells;
* the exact same-moment gain and zero-target continuum ghost;
* six nested noisy packets, 14 malformed-input hostiles, and the missing-atom
  energy and degree-one monotonicity failures.

Reproduce from the repository root:

```text
python 04-computation/experiments/collatz_refinement_energy_dual_20261005.py
python -O 04-computation/experiments/collatz_refinement_energy_dual_20261005.py
```

Both modes perform **8,144 exact checks**. The computation checks finite
instances; the product/sign/tail proofs establish the all-length claims.
