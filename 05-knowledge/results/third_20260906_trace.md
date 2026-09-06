# Sharp carried trace and norm jets, and a trace-convolution obstruction

**Status: PROVED in the stated all-height boundary scope;
[independent analytic and resultant audit accepted](third_20260906_trace_root_audit.md).** The all-height statements below are proved algebraically. The
root locations through height 7 are **FINITE-EXACT**, and universal trace
positivity and trace real-rootedness remain **OPEN**. No theorem ID or
external-priority claim is made.

The positive result is a sharp all-height refinement of the
[forced-factor theorem](third_20260906_laurent.md): its trace factor and
its final norm factor have exactly the predicted multiplicities. The trace
has a strictly positive first derivative, and deleting the lower carry
changes this derivative by less than one twelfth. Nevertheless, the same
deletion changes the norm's vanishing order and sign. The fully deflated
norm consequently has a four-periodic boundary sign. This rules out a
specific positive convolution of the trace polynomials at every height
congruent to 2 or 3 modulo 4.

These statements concern the polynomial continuation at a negative
parameter. They do not prove or refute noncancellation for an actual
support. The actual family has `x>=1`, whereas the boundary studied here
is `x=-h`.

## 1. Inheritance, exact search, and the retained board

The closest mechanism is the all-height carry degeneration in
[the first-round bundle](third_20260906_laurent.md), Sections 2–4.
The endpoint-33 and endpoint-39 positive characteristic certificates are
consumers of full actual normalization, rather than general trace
theorems. The older
[carry-corrected Hermite trace formulation](trinomial_trace_sign_empty_core_certificate_sep06.md)
already explains why a negative scalar trace is insufficient to make the
whole response negative. That distinction is retained here.

The canonical hostile is the
[actual full integer-power observer-cone separator](continuing1_20260906_laurent_cone_separator.md).
The corrected near misses are omission of the inverse carry and promotion
of individual real-rooted rows to a cross-row sign. The least-used sidecar
is the distinction between the first and last characteristic coefficient
at the same ramified parameter fiber. Targeted searches for trace
positivity, factorial/Newton traces, characteristic convolution, the
literal boundary, and the formulas below found no earlier exact-jet or
four-periodic trace-convolution statement. General trace and norm
operations are inherited.

The six-concept board is: the actual complete rows; their reversed
logarithm; the trace; the norm; the negative-parameter root cluster; and
positive polynomial convolution. The map takes the carried quotient
operator to its first and last characteristic coefficients. It preserves
the original coefficient normalization and the carry. It loses the
individual phase responses; the ramified root cluster is the needed
sidecar. The cheap tests were height 7, noninteger parameters, deletion of
the full carry, and a putative quadratic three-term recurrence.

The trace root pattern survived those finite height tests. The low-order
recurrence failed, and the norm jet supplied a structural obstruction to
the proposed positive convolution. The latter is the stopping reason for
this route; a larger root bank is not offered as a proof.

## 2. Rows, quotient, and an all-height logarithmic trace formula

Fix an integer `h>=1`, put `y=x+h`, and let `(a)_m` denote a falling
factorial. Write

```text
p(t)=sum_(j=0)^h
  (2h+1)! (y)_(h-j) / ((3h-3j)!(1+2j)!) * t^j,

q(t)=sum_(e=-1)^(2h)
  (2y)_(2h-e) / ((6h-3e)!(2+2e)!) * t^e.                 (1)
```

These are exactly the normalized full rows of the first-round bundle.
Initially work in `Q(y)[t]/(p)`. The same carry cancellation as in that
bundle makes the remainder `R=q mod p` a polynomial in `Q[y,t]`. Let
`T` be multiplication by `R` in the basis `1,t,...,t^(h-1)` and set

```text
det(zI-T)=z^h+sum_(k=1)^h c_(h,k)(y-h) z^(h-k),
T_h(y)=c_(h,1)(y-h),        N_h(y)=c_(h,h)(y-h).          (2)
```

The symbol `T_h` denotes a scalar trace polynomial; the matrix is `T`.
Use `q_e` and `p_j` for the coefficients in (1). Reverse `p`:

```text
P(u)=u^h p(1/u)=1+sum_(m=1)^h a_m(y)u^m,
a_m(y)=(2h+1)! (y)_m / ((3m)!(2h+1-2m)!),
ell_e(y)=[u^e] log P(u).                                (3)
```

Then, identically in `Q[y]`,

```text
T_h(y)=q_(-1)*p_1/p_0 - h*q_0
       +sum_(e=1)^(2h) e*q_e*ell_e.                     (4)
```

The first product in (4) is polynomial; its displayed ratio is only a
convenient expression in the fraction field. For an explicit finite sum,

```text
ell_e=sum_(l=1)^e (-1)^(l-1)/l
       *sum_(m_1+...+m_l=e, 1<=m_i<=h) product_i a_(m_i).
```

Equivalently, set `a_e=0` for `e>h` and compute
`ell_e=a_e-(1/e)sum_(j=1)^(e-1) j*ell_j*a_(e-j)`.
This is a lawful exact recurrence, but it is a signed recurrence. It is
not a positive convolution.

To prove (4), let `rho_i` be the generic roots of the monic `p`.
Newton's logarithmic identity gives
`sum_i rho_i^e=-e*ell_e` for `e>=1`, while
`sum_i rho_i^(-1)=-p_1/p_0`. Taking the negative trace of `q` gives (4).
The equality holds over the fraction field and then extends as a
polynomial identity. No simplicity assumption at a degenerate parameter
is required.

## 3. Positive trace slope and exact multiplicity for every height

Define positive rational constants

```text
d_h=(2h-1)!/(6h)!,
b_h=2(2h)!/(6h+3)!,
v_h=(3h)(3h-1)(3h-2)/6.
```

For every integer `h>=2`,

```text
T_h(0)=0,
T_h'(0)=h*d_h*(1-L_h),
L_h=4*v_h / ((h-1)(6h+1)(6h+2)(6h+3)),
0<L_h<1/(12(h-1))<=1/12.                               (5)
```

Thus the factor `x+h` in the trace has **exact multiplicity one for all
heights h>=2**, rather than only the heights checked in the first-round
certificate. In particular `T_h(y)>0` on some interval `0<y<epsilon_h`.
The deflated trace `U_h(y)=T_h(y)/y` is strictly positive at `y=0`.

Here is the proof, including the terms that disappear. Each `a_m`,
`1<=m<=h`, is divisible by `y`; hence `ell_e=O(y)` for every positive e.
For `e>h`, every composition in (3) has at least two parts, so
`ell_e=O(y^2)`. The coefficient `q_e` is divisible by y when
`e<2h`. Consequently every term in the sum in (4) is `O(y^2)`, including
the final term with e=2h. The derivative is entirely the competition
between the constant channel and the lower carry.

Directly from (1),

```text
q_0=-d_h*y+O(y^2),
q_(-1)=b_h*y+O(y^2),
p_1/p_0=v_h/(y-h+1).                                   (6)
```

For h>=2 the last denominator is nonzero at zero. Inserting (6) into
(4) gives `T_h'(0)=h*d_h-b_h*v_h/(h-1)`, which is (5).
Finally, `4v_h<18h^3` and
`(6h+1)(6h+2)(6h+3)>216h^3`, proving the strict bound.

The height-one exception is explicit:

```text
T_1(y)=(884y^2+123y+1)/90720.                            (7)
```

It has no zero at y=0. This exception is needed in the semiring statement
below: `T_1`, as well as every `U_h`, is positive near zero.

## 4. Sharp norm jet and the effect of deleting the carry

Set

```text
a_h=(-1)^(h-1)(h-1)!(2h+1)!/(3h)!.
```

For every integer `h>=1`, the full carried norm satisfies

```text
N_h(y)=(b_h^h/a_h)*y^(h-1)+O(y^h).                      (8)
```

Thus the factor `(x+h)^(h-1)` in the determinant has **exact
multiplicity h-1 for every height**. Let `q^+=q-q_(-1)t^(-1)` and define
`N_h^+` from its quotient multiplication operator in the same way. Then

```text
N_h^+(y)=d_h^h*y^h+O(y^(h+1)).                          (9)
```

For h>=2 the no-carry trace has derivative `h*d_h`, by (4) and (6).
The carry decreases this derivative by less than one twelfth, while
deleting it changes the norm's exact vanishing
order and changes its leading sign at every even height. Trace stability
does not control the norm.

For (8), use the convergent local root parameter, or equivalently formal
Puiseux series, `y=epsilon^h`. At y=0, `p(t)=t^h`; moreover

```text
p(epsilon*v)/epsilon^h = v^h+a_h+O(epsilon).
```

The h roots of `v^h+a_h` are nonzero and distinct. The implicit-function
theorem therefore labels the roots as
`t_i=epsilon*v_i(epsilon)`, where
`v_i(0)^h=-a_h`. Along each root, the full Laurent expression gives

```text
q(t_i)=(b_h/v_i(0))*epsilon^(h-1)+O(epsilon^h).
```

The inverse carry is the leading term. Multiplying, using
`product_i v_i(0)=(-1)^h*a_h`, and including the characteristic constant
term sign `(-1)^h` proves the leading coefficient in (8). The full
expression is already a polynomial in y, so its remaining exponents are
integers; the next order is at least h.

After removing the carry, (6) instead gives
`q^+(t_i)=-d_h*y+O(epsilon^(h+1))`: every positive-index term has higher
epsilon order, including the top channel of order `epsilon^(2h)`.
Their product proves (9), again using polynomiality to recover integer
y orders. This argument specializes the polynomially extended operator,
not an undefined inverse of t at y=t=0.

## 5. A four-periodic obstruction after full deflation

The inherited all-height factor theorem gives

```text
D_h(y)=product_(r=2)^h (y+r-h)^(r-1),
B_h(y)=N_h(y)/D_h(y) in Q[y].
```

Let `E_h=product_(r=2)^(h-1)(r-h)^(r-1)`, with empty products equal to
one. Formula (8) gives the exact nonzero endpoint

```text
B_h(0)=b_h^h/(a_h*E_h),
sign B_h(0)=(-1)^(h(h-1)/2).                            (10)
```

Consider the semiring of finite polynomials with nonnegative real
coefficients in the generators

```text
y, T_1(y), and U_j(y)=T_j(y)/y for integers j>=2.         (11)
```

For every `h congruent to 2 or 3 modulo 4`, the fully deflated determinant
`B_h` does **not** belong to this semiring. By (5) and (7), each generator
is nonnegative at zero. Every finite polynomial in those generators
with nonnegative coefficients is therefore nonnegative at zero,
contradicting (10). Allowing also the undeflated traces adds nothing,
because `T_j=y*U_j` for j>=2.

This is an all-height obstruction to a specified trace-to-exterior
positive convolution, including arbitrary finite products and sums. It
does not exclude a representation in the original positive coordinate
x, parameter shifts, signed terms, other exterior information, or a
representation valid only on the actual domain `y>=h+1`. It does not
claim an actual carried noncancellation counterexample.

## 6. Finite controls, a recurrence failure, and the remaining question

The [script](../../04-computation/third_20260906_trace.py) computes the
trace in two ways: quotient multiplication and the reversed logarithm.
For exactly `h=1,...,7`, all trace coefficients in y are nonnegative;
after removing the zero factor at h>=2, exact Sturm counts place all
remaining roots strictly in `(-1/4,0)`, with no repeated root. This is
**FINITE-EXACT evidence**, not an all-height real-rootedness theorem.
At the noninteger parameter `x=1/2` it also checks positive trace, all
characteristic coefficients strictly positive, and h simple negative
first-row roots, without claiming an actual integer support there.

The norm jets, full deflation, and no-carry comparison are independently
computed from exact multiplication determinants for `h=1,...,4`.
At the exact hostile `h=2, y=1/100000`, the trace is positive and the
norm is negative. The output records both rational values. This hostile
is outside the physical parameter domain and is used only against
claims on all `y>0`.

There is also no quadratic-coefficient three-term recurrence of the
following natural form, even with arbitrary signed scalar coefficients.
Let `V_h` be the monic normalization of `T_h`. Already at h=3 there are
no real constants a,b,c,d,e such that

```text
V_3=(y^2+a*y+b)*V_2+(c*y^2+d*y+e)*V_1.                 (12)
```

The six coefficient columns of
`y*V_2, V_2, y^2*V_1, y*V_1, V_1, V_3-y^2*V_2`, in degrees 0 through
5, have nonzero determinant. Reducing their rational entries modulo
101 (all denominators are invertible) gives

```text
 0  0  0  0  4  0
 0 11  0  4 88 33
11 93  4 88  1 46
93  5 88  1  0 98
 5  1  1  0  0 58
 1  0  0  0  0 81
```

Its determinant is 29 modulo 101. Thus the target is outside the span
of the five proposed recurrence columns over Q and over R. This
rejects (12), not every higher-order or differential recurrence.

The useful surviving question is whether `U_h` is positive, or even has
all its zeros in `(-1/4,0)`, for every height. The explicit signed
logarithmic formula makes that question precise. Neither that question
nor positivity of all characteristic coefficients on `x>=1` is proved
here. The new all-height facts are the sharp trace and norm jets, their
carry comparison, and the semiring obstruction (10)–(11).

Reproduce with:

```bash
python3 04-computation/third_20260906_trace.py
python3 -O 04-computation/third_20260906_trace.py
```

The stored output is
[third_20260906_trace.out](third_20260906_trace.out). Arithmetic is
rational throughout, checks are always active, and no repository
producer is imported. The first-round Laurent computational certificate
is unchanged by this bundle.

Normal and optimized runs pass 107 exact gates and byte-match. Raw LF
SHA256 values are:

```text
source 57da4695237e57ecbcc460d5173f1a47f59801ef4a87e08ddeb753ecaaa889e0
output 26995e55f5a2cbe53ca290d2d8705af4389b8274b8779638ec8ecfe326ad8660
semantic ab0029e4de0bc4abad2b43eab34b66b280800b7a34a063965135b1639892e5e9
```
