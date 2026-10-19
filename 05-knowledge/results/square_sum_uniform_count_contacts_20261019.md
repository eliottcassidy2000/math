# Polynomial-degree contact heights uniform in the requested neighbor count

**Status: PROVED; independent proof audit passed.** This sharpens
only the dependence on the requested number of neighbors in the fixed-count
atlas. It does not claim polynomial complexity in the binary encoding length
of that number. No long-chain gluing, narrow-window, or other target is added.

## 1. Statement and inherited facts

Let `T=(a,b,C)` be an admissible primitive positive Pythagorean triple, meaning
`C²>4ab`, and let `epsilon in {-1,+1}`. Write

\[
 L=1+\lceil\log_2 C\rceil,\qquad \ell_r=\log_2(r+1).
\]

**Theorem.** A deterministic exact algorithm, on input `T,epsilon,r>=1`,
returns `r` pairwise distinct admissible primitive petals `U_i`, all different
from `T`, with positive rational certificates `rho_i` satisfying

\[
 \rho_i^2=1+\varepsilon f_T-\varepsilon f_{U_i},
 \qquad f_{(a,b,C)}=2ab/C^2.
\]

All coordinates and certificates have at most

\[
 O\!\left(r^8\ell_r^2(L+\ell_r)\right)                 \tag{1}
\]

bits, with an absolute effective implied constant. The total running time
is polynomial in `r` and `L`. In particular, an intentionally conservative
schoolbook bound is

\[
 O\!\left(r^{29}\ell_r^7(L+\ell_r)^3\right).           \tag{2}
\]

The required output count `r`, rather than only `log r`, occurs in this
bound. The total encoded output has at most `r` times the bound in (1).
For numerical heights one may also write

\[
 C_i\le C^{,O(r^8\ell_r^3)},                          \tag{3}
\]

because `C>=5`, `L=O(log C)`, and `ell_r>=1`. No small or optimal exponent
is asserted.

We inherit the following proved ingredients, without changing their scope:

1. [Fixed-count contact atlases](square_sum_fixed_count_contacts_20261019.md),
   Sections 2–3: low-cone normalization; the odd cusp sections and their
   first derivatives; the r-point signed recurrence lemma; exact acceptance,
   primitive normalization, and uniqueness of the petal represented by a
   spread.
2. [Polynomial-bit signed contact construction](square_sum_contact_output_height_20261019.md),
   Sections 4–6: the normalized division-polynomial recurrences, weighted
   degrees, the coefficient-norm bound, and the uniform orbit-height bound
   `h(Q+qP)=O(q²[L+log(|q|+1)])`, including ordinary exceptional branches.
3. The plus and minus bounded-section notes supply the explicit section
   coordinates used below. They are also displayed here.

The simpler global r-point recurrence already gives `q=O(r²L)` and hence
a polynomial-time r-neighbor algorithm with per-output bits
`O(r⁴L²[L+log(r+1)])`. The new contribution here is the polynomial-degree
numerical-height bound (3), retaining linear dependence on the input-bit
length in (1), with polynomial rather than unspecified dependence on r.
It controls the cusp atlas uniformly in `r`, including the cost of
constructing its polynomial coefficients. A bounded number of field
operations for fixed r alone would not establish this stronger result.

## 2. Polynomial input sections on a short model

Use the inherited low representative `t` with

\[
 1<t<\cot(5\pi/24)<4/3,\qquad
 f(t)=\frac{4t(t^2-1)}{(t^2+1)^2},\qquad H(t)\le\sqrt{2C}.
\]

Put

\[
 D=t^2+1,\quad A_+=t^2+2t-1,\quad A_-=-t^2+2t+1,
 \quad k=A_\varepsilon/D,\quad K=k^2.
\]

On `E_t:Y²=(X-K)(X²-1)`, write `P=(0,k)` and set

\[
 B_+=Q+P=\left(\frac{t^2-4t+1}{D},
                   -\frac{8t^2(t-1)}{D^2}\right),
\]
\[
 B_-=Q_-=\left(\frac{t^2-4t+1}{D},
                   -\frac{8t(t-1)}{D^2}\right).
                                                               \tag{4}
\]

Here `Q` is the positive-sheet plus anchor, while `Q_-` is the
negative-sheet minus anchor. For odd `d>=3`, the inherited cusp sections are

\[
 R_d^+=dB_+-P=dQ+(d-1)P,\qquad R_d^-=dB_- .            \tag{5}
\]

First set `x=D²X, y=D³Y`, then

\[
 \xi=9x-3A_\varepsilon^2,\qquad \eta=27y.
\]

The resulting short equation is

\[
 \eta^2=\xi^3+\alpha\xi+\beta,
\quad \alpha=-81D^4-27A_\varepsilon^4,
\quad \beta=-54A_\varepsilon^6+486A_\varepsilon^2D^4.   \tag{6}
\]

The images of `B_epsilon` have integer polynomial coordinates

\[
 \xi_B=9D(t^2-4t+1)-3A_\varepsilon^2,
\]
\[
 \eta_B=\begin{cases}
 -216t^2(t-1)D,&\varepsilon=+1,\\
 -216t(t-1)D,&\varepsilon=-1.
 \end{cases}                                         \tag{7}
\]

Also `P_s=(-3A_epsilon²,27A_epsilon D²)`. All these are fixed polynomials
with absolute degree and coefficient-norm bounds. With weights
`wt(xi)=2, wt(alpha)=4, wt(beta)=6, wt(eta)=3`, their ordinary degrees are
at most twice their weights. Work from now on in `h=t-1`, substituting
`t=1+h` in these fixed input polynomials. Their degrees and coefficient
norms still have absolute bounds.

The generic curve over `Q(h)` is nonsingular: its original roots `K,1,-1`
are distinct rational functions. Its singular special fiber at `h=0`
will only be used to evaluate explicit polynomials, never as a nonsingular
elliptic-curve input to a multiplication theorem.

## 3. Uniform univariate coefficient bounds with no polynomial reduction

Let `F_n` be the predecessor's normalized short-Weierstrass division
polynomial: `psi_n=F_n` for odd n and `psi_n=2 eta F_n` for even n. It proves

\[
 \log\max(1,\|F_n\|_1)\le4n^2\log(n+1),\qquad
 \operatorname{wtdeg}F_n=\begin{cases}n^2-1,&n\text{ odd},\\
 n^2-4,&n\text{ even}.\end{cases}                       \tag{8}
\]

Define the *evaluated polynomial*, not a reduced rational function,

\[
 G_n(h)=F_n(\xi_B(h),\alpha(h),\beta(h))\in\mathbb Z[h].
\]

Let `H>=1` bound the 1-norms of the three fixed input polynomials for both
ports. Every monomial of `F_n` has total degree at most `n²`, so coefficient
norm submultiplicativity gives directly

\[
 \|G_n\|_1\le\|F_n\|_1 H^{n^2},\qquad
 \deg G_n\le2n^2,
\]

and hence

\[
 \deg G_n=O(n^2),\qquad
 \log\max(1,\|G_n\|_1)=O(n^2\log(n+1)).               \tag{9}
\]

These are bounds on actual integer polynomial arrays before any rational
cancellation. A zero polynomial is harmless in the norm assertion.

For odd d, standard multiplication formulas become

\[
 \xi(dB)=\frac{\xi_B G_d^2-4\eta_B^2G_{d-1}G_{d+1}}{G_d^2},
 \qquad
 \eta(dB)=\frac{\eta_B G_{2d}}{G_d^4}.                 \tag{10}
\]

The second identity follows from `psi_{2d}=2 psi_d omega_d` and
`eta(dB)=omega_d/psi_d³`; equivalently it is
`eta(dB)=psi_{2d}/(2psi_d^4)`. These formulas, together with the definitions
of `omega_d` and `psi_d`, are in Sutherland's
[MIT 18.783, Lecture 5, Sections 5.5–5.6, Lemma 5.20 and Theorem 5.21](https://ocw.mit.edu/courses/18-783-elliptic-curves-fall-2025/mit18_783_f25_lec05.pdf).
They are applied on the nonsingular generic curve over `Q(h)`.

Each numerator and denominator in (10) has degree `O(d²)` and log 1-norm
`O(d² log(d+1))`. Subtract `P_s` once on the plus side by the ordinary
short-model addition formula, undo the model change, and apply

\[
 u=\frac{kY+X-K}{X(X-K)},\qquad
 z=\frac{2u^2X-K-2u-Ku^2}{k},\qquad v=\varepsilon u.   \tag{11}
\]

Use polynomial pairs with cross-multiplication for each rational operation;
do not take a polynomial gcd. There are an absolute number of operations
per d after (10). Thus the resulting unreduced integral pairs

\[
 v_d=N_d/E_d,\qquad z_d=Z_d/T_d                        \tag{12}
\]

still satisfy, for every constituent polynomial p,

\[
 \deg p=O(d^2),\qquad
 \log\max(1,\|p\|_1)=O(d^2\log(d+1)).                 \tag{13}
\]

Indeed multiplication adds degrees and log norms; addition costs only a
factor two after cross-multiplication. A fixed number of these operations
changes only absolute constants. In particular, this is not an appeal to
an unproved norm bound for arbitrary polynomial factors.

### All these denominator constants are nonzero

At h=0, for both ports, the short cubic and B are

\[
 \eta^2=(\xi-24)^2(\xi+48),\qquad B_0=(-48,0),
 \qquad P_{s,0}=(-12,216).
\]

The cubic derivative at `xi=-48` is `5184`, so
`G_1(0)=1` and `G_3(0)=-5184²` are nonzero. In the odd normalized
recurrence, the factor `16(xi³+alpha xi+beta)²` vanishes at B_0. The
surviving term is a signed product of smaller odd-indexed G's. Induction
therefore gives `G_d(0)!=0` for every odd d. Formula (10) consequently
specializes to `dB_0=(-48,0)`: its x correction has the factor `eta_B²`,
and its y numerator has the factor `eta_B`.

For plus subtraction, the two x-coordinates at zero are -48 and -12,
so its slope denominator is nonzero. Back on the original model, the
resulting R_d specializes to `(3,4)` for plus and `(-1,0)` for minus.
Both have `X(X-K)!=0`; also `D(0)=A_epsilon(0)=2`. Thus all raw denominator
polynomials in (12) have nonzero constant term. No cancellation at h=0 is
even needed in (12).

## 4. A dyadic cusp cutoff with polynomial bit length

The inherited formal-coordinate lemma says

\[
 v_d(1+h)=1+dh+O(h^2),\qquad z_d(1+h)=2\varepsilon+O(h).
                                                               \tag{14}
\]

For `d=3,5,...,2r+1`, construct polynomial pairs for

\[
 \frac{v_3-t}{h},\qquad
 \frac{v_{d+2}-v_d}{h}\ (3\le d\le2r-1),\qquad
 \frac54-v_{2r+1},                                    \tag{15}
\]

and include the denominators of every v_d and z_d. Form each difference
by cross-multiplication in the unreduced pairs (12). In the first two
kinds, the numerator has exactly one factor h, by (14) and the nonzero
denominator constants. Remove this single factor by shifting its integer
coefficient array. Nothing else is factored or canceled.

All numerators and denominators from (15) have nonzero integer constant
term. Their rational functions at zero have values `2,2,...,2,1/4`,
respectively. Choose signs separately so every polynomial under consideration
has positive constant term. Because each rational function has positive
value at zero, its numerator and denominator were changed by the same sign.
The denominator-only polynomials can be signed independently.

Let M be the largest 1-norm among these O(r) polynomials and define

\[
 b=\max\{4,\ 1+\operatorname{bitlength}(M)\},\qquad
 \eta_r=2^{-b}.                                       \tag{16}
\]

Their degrees are O(r²) and log norms are O(r² log(r+1)), by (13) and a
constant number of additional operations. Hence

\[
 b=O(r^2\ell_r),\qquad
 \eta_r\ge2^{-O(r^2\ell_r)}.                           \tag{17}
\]

For `0<=h<=eta_r<=1`, every signed polynomial satisfies

\[
 p(h)\ge p(0)-h\sum_{i\ge1}|p_i|
       \ge p(0)-\eta_r M\ge p(0)-\tfrac12
       \ge\tfrac12p(0)>0.
\]

It follows that, on `0<h<=eta_r`, all inverse denominators are nonzero and

\[
 1<t<v_3<v_5<\cdots<v_{2r+1}<5/4.                     \tag{18}
\]

The inherited monotonicity of f and `f(5/4)<1/2` produce r different
admissible neighbors, all different from the input. The quartic identity
gives the positive certificates `abs(z_d)/(v_d²+1)`.

For `eta_r<=1/16`, put `delta_r=f(1+eta_r)`. Directly,

\[
 \frac{f(1+h)}h=\frac{4(1+h)(2+h)}{((1+h)^2+1)^2}>1
 \quad(0<h\le1/16),
\]

since numerator >=8 and denominator <=`(41/16)²<8` even when
`1+h<=5/4`. Therefore `delta_r>=eta_r`, and the explicitly computable choice

\[
 J_r=b+2                                                \tag{19}
\]

satisfies `2^{J_r}delta_r>=4>5/2`. No root isolation or real integration is
required to compute this cutoff or J_r.

## 5. Uniform recurrence bound and output height

On `t>1+eta_r` the low-cone monotonicity gives `f(t)>delta_r`. The inherited
r-point recurrence argument applies with the positive-z anchor for both
ports; for minus this switches from the negative cusp sheet. Its guaranteed
spread-monotone admissible arc has invariant length greater than `1/20`.
Set

\[
 W_r=\frac57(4+J_r),\quad N_r=\lfloor20rW_r\rfloor,
 \quad M_r=2rN_r.                                     \tag{20}
\]

The candidate set `Q+qP`, with q nonzero even and `|q|<=M_r`, contains at
least r different successful spreads. Equations (17), (19), and (20) give

\[
 M_r=O(r^2J_r)=O(r^4\ell_r).                           \tag{21}
\]

Scan those candidates exactly using the original curve, ordinary group-law
branches, and the predecessor's acceptance test. Store the at most r
accepted spreads and compare each new one with them. The successful open
arc avoids the excluded inverse charts, so skipping those charts does not
invalidate termination. Positive-sheet certificates are `z/(u²+1)`.

The output-height predecessor, equation (14), gives every finite candidate
and every running multiple of P logarithmic height

\[
 O\bigl(M_r^2[L+\log(M_r+1)]\bigr)
   =O\!\left(r^8\ell_r^2(L+\ell_r)\right).             \tag{22}
\]

The fixed-degree inverse map, rational certificate, and primitive
normalization preserve this order. Repeatedly reduce rational numbers
after arithmetic operations; thus intermediate integers within a fixed
addition or inverse-map formula have the same asymptotic bit bound.

For the cusp branch, evaluate the raw polynomial pairs (12) at the rational
h=t-1. Their degrees and coefficient lengths give the smaller bound

\[
 O\bigl(r^2[L+\ell_r]\bigr)                            \tag{23}
\]

on output and certificate bits: clear the input denominator by degree,
use `h(h)=O(L)`, and bound by the coefficient 1-norm. Evaluating unreduced
pairs cannot increase this bound; reducing afterward can only decrease it.
Together (22) and (23) prove (1) and (3).

## 6. Explicit polynomial-time construction of all coefficient data

There is no symbolic-expression-efficiency assumption hidden above. Compute
all G_n up to `4r+2` for each port directly in `Z[h]`, substituting the
fixed polynomial inputs (6)–(7) into the predecessor's normalized polynomial
recurrences. Each recurrence uses a constant number of polynomial additions
and multiplications and smaller indices, so an increasing-index table works.
There are O(r) entries. Equation (9) bounds their degrees by O(r²) and
coefficient bits by O(r² ell_r). Each intermediate product has the same
orders, because a recurrence has only a fixed number of factors with indices
O(r). Addition and cancellation in coefficient arrays cannot increase the
1-norm beyond the triangle-inequality bound.

Dense arrays therefore suffice. Schoolbook polynomial convolution takes
O(r⁴) integer products per multiplication, on O(r² ell_r)-bit integers;
schoolbook integer multiplication gives a conservative O(r⁸ ell_r²) bit
cost per polynomial multiplication. O(r) such fixed-size recurrences,
section constructions, and cutoff polynomials cost O(r⁹ ell_r²) bits.
Their tables require O(r⁵ ell_r) bits. Only multiplication, addition,
coefficient inspection, and the exact removal of the known factor h are
used; there are no polynomial gcds or polynomial factorization calls.

Horner evaluation of O(r) polynomials of degree O(r²), followed by rational
reduction and normalization, also takes polynomial time in r and L.
A crude bound O(r⁹(L+ell_r)³) suffices, using the unamortized integer-gcd
bound O(B³) and operands B=O(r²(L+ell_r)).

The recurrence scan uses O(M_r) group additions and O(r M_r) rational
operations/comparisons including duplicate checks. With the operand bound
`B=O(r⁸ ell_r²(L+ell_r))`, the same O(B³) reduced-rational-operation
bound gives

\[
 O(rM_rB^3)=O(r^{29}\ell_r^7(L+\ell_r)^3),
\]

which dominates the stated preprocessing and cusp costs. This proves (2).
The memory needed by the coefficient tables, arithmetic, and retained
spreads is also polynomial; O(r B+r⁵ ell_r) bits suffice besides the
read-only input and returned outputs.

## 7. Scope and verification contract

This is a one-anchor, prescribed-port, r-neighbor theorem with exact
certificates. Its effective constants are absolute. In particular it
replaces the fixed-count note's arbitrary r-dependence for this task.
It makes no assertion that repeatedly applying this theorem along a growing
chain retains polynomial complexity in the chain length. It does not add
simultaneous constraints, prescribed narrow spread windows, or a graph
labeling conclusion beyond the inherited single-anchor contacts.

The proof relies on the proved fixed-count formal lemma and the proved
uniform division-polynomial bound. The companion checker supplies finite
implementation controls and does not replace these uniform proofs:

- Direct univariate integer-polynomial recurrences for n=1,...,14 on both
  ports verify the polynomial model, degree bound, substituted coefficient
  norm bound, and every odd denominator's nonzero cusp value.
- Exact first jets at d=3,5,7 on both ports verify all divisions are by units,
  the odd-multiple specializations, `v_d(0)=1`, `v_d'(0)=d`, and
  `z_d(0)=2 epsilon`. No polynomial reduction is used for these jets.
- Odd division-polynomial multiples are checked against independent ordinary
  group additions at h=1/64,1/1024,1/8, with exact contact certificates.
  The smallest h gives six primitive admissible outputs. The larger probes
  include failures of admissibility, controlling against treating the cusp
  construction as globally admissible.
- Four explicit signed-coefficient examples test the exact dyadic positivity
  cutoff and `J=b+2`; these are cutoff-lemma controls rather than a computed
  arbitrary-r atlas.

Normal and optimized Python runs are byte-identical, with **213 checks**:
94 exact-integer, 2 polynomial-symbolic, 53 exact first-jet, and 64
exact-rational checks. The largest evaluated division polynomial (n=14)
has degree 336; the coefficient 1-norms have 863 bits on plus and 785
bits on minus. These observations are finite controls, not asymptotic
estimates.

Run, with SymPy installed:

```bash
python3 04-computation/experiments/square_sum_uniform_count_contacts_20261019.py
python3 -O 04-computation/experiments/square_sum_uniform_count_contacts_20261019.py
```
