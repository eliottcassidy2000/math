---
id: THM-3576
title: "Higher-exponent Belyi Keller collision tower"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For every integer n at least two, an explicit rational Keller pair has a
  maximal polynomial-observable completion on the smooth surface
  c^n e=b(b-beta_n).  The resulting polynomial map from A2 is etale,
  generically of degree n(n-1)+1, noninjective, and misses exactly one target
  point.  Its physical symplectic form is exact and the constant function is
  a sum of two polynomial Poisson brackets, while every homogeneous and
  two-by-two or two-by-three weight cell is empty; every Darboux pair needs
  at least six homogeneous pieces.  The degree-d_n inverse carrier has
  cited monodromy A_(d_n) for odd n and S_(d_n) for even n; the n=3 member
  is the first nonsolvable case, of degree seven.  A single polynomial Darboux pair on any
  target would give a planar Jacobian counterexample; no such pair and no
  counterexample to JC(2) is claimed.
source: root, 2026-08-20
audit: >
  Independent hostile audit rederived the critical-value completion,
  maximal observable ring, exact image and fibre invoices, all homogeneous,
  two-by-two, and two-by-three weight obstructions including every small
  boundary and homogeneous side branch, and the cited Shabat/monodromy
  specialization.  Ordinary and optimized exact replays are byte-identical.
depends_on: []
related:
  - THM-3566-chebyshev-pell-odd-keller-collision-tower
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
external:
  - "Cameron, Kemp, Maslak, Melamed, Moy, Pham, and Wei, Shabat Polynomials and Monodromy Groups of Trees Uniquely Determined by their Passport, Involve 12 (2019), 791-812, Propositions 2 and 7, DOI 10.2140/involve.2019.12.791."
script: 04-computation/jc2_higher_exponent_belyi_keller_thm3576.py
output: 05-knowledge/results/jc2_higher_exponent_belyi_keller_thm3576.out
script_sha256: 3413997e87a787f4e147f0f8afc56192b7a6ba80229f70c280f3027882d61a45
output_sha256: d9b0ff5b7d10213cf403edeee6375c6069cd3c9855a864c481152f905688cab9
hash_basis: LF-normalized bytes
---

# THM-3576 -- higher-exponent Belyi Keller collision tower

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This theorem does **not** refute `JC(2)`.  It builds a new tower of exact
near-counterexamples whose boundary exponent, collision passport, and
monodromy differ from the exponent-two Chebyshev family.

All varieties are over `C`.  Fix an integer `n>=2` and put

```text
R(t)=1+t,
P_n(t)=sum_(j=0)^(n-1) binom(n-1,j)t^j/(nj+1),
B_n(t)=t P_n(t)^n.                                    (1)
```

## 1. The hypergeometric/Belyi identity

Coefficient comparison gives

```text
P_n+n t P_n'=R^(n-1),
B_n'=P_n^(n-1) R^(n-1).                               (2)
```

Equivalently,

```text
P_n(t)=integral_(s=0)^1 (1+t s^n)^(n-1) ds.           (3)
```

In particular `P_n(0)=1` and `P_n(-1)!=0`.  Define

```text
beta_n=B_n(-1)=-P_n(-1)^n,
W_n(t)=(B_n(t)-beta_n)/R(t)^n.                         (4)
```

Equation `(2)` shows that `W_n` is a polynomial and

```text
W_n(-1)=P_n(-1)^(n-1)/n !=0.                          (5)
```

Every root of `P_n` is simple: if a root had multiplicity `m`, the two
sides of `(2)` would give orders `nm-1` and `(n-1)m`, hence `m=1`.
The roots of `W_n` are also simple and disjoint from `P_nR`.  Indeed a
multiple root of `W_n` would be a critical point of `B_n`; `(2)` would put it
on `P_nR`, where the critical values are respectively `0,beta_n`, not a new
root of `W_n`.

Thus `B_n` is a two-finite-critical-value polynomial of degree

```text
d_n=n(n-1)+1.                                         (6)
```

Its passports over `0,beta_n,infinity` are

```text
(1,n^(n-1)),          (n,1^((n-1)^2)),          (d_n). (7)
```

This is exactly Adrianov's size-one Shabat family, rather than merely a
passport match.  In Cameron--Kemp--Maslak--Melamed--Moy--Pham--Wei,
Proposition 2, set their parameters `r=s=n,t=1` and then put `z=1+t`.
If `S_n(z)=sum_(k=0)^(n-1)(1/n)_k z^k/k!`, direct coefficient comparison
gives

```text
P_n(t)=S_n(1+t)/S_n(1),
B_n(t)=-F_(n,n,1)(1+t)/S_n(1)^n.                      (7a)
```

Their Proposition 7 therefore supplies the exact cited monodromy group

```text
Mon(B_n)=A_(d_n) for n odd,       S_(d_n) for n even. (7b)
```

The completion below is new data layered over that one-variable inverse
carrier; the theorem makes no novelty claim for `(1)--(7b)`.

## 2. A constant-Jacobian rational pair

On `A2_(x,q)` put

```text
t=x^n q,
a=q/R(t)^n,
c=x P_n(t)R(t),
b=B_n(t),
e=q W_n(t).                                           (8)
```

For arbitrary rational functions `H,G` of `t=x^nq`, direct differentiation
gives the master identity

```text
Jac_(x,q)(qH,xG)=-(t H G^n)'/G^(n-1).                 (9)
```

Taking `H=R^(-n)` and `G=P_nR`, equations `(2)` and `(9)` give

```text
Jac_(x,q)(a,c)=-1,
b=a c^n,
e=a(b-beta_n).                                        (10)
```

Although `a` is rational, `(b,c,e)` is polynomial and satisfies

```text
c^n e=b(b-beta_n).                                    (11)
```

## 3. The maximal polynomial-observable ring

The completion in `(11)` is not an arbitrary small choice:

```text
C[a,c] intersect C[x,q]=C[b,c,e]                      (12)
```

inside `C(x,q)`.  To prove this, use the torus weights

```text
wt(x)=wt(c)=1,             wt(q)=wt(a)=-n,             (13)
```

with `t,b` invariant.  A weight `-s<0` element of `C[a,c]` has the unique
form `c^(-s)f(b)`, where `b^m|f` and `m=ceil(s/n)`.  Pullback to `C[x,q]`
has a possible pole of order `s` along `R(t)=0`.  Since
`ord_R(B_n-beta_n)=n`, regularity is equivalent to
`(b-beta_n)^m|f`.  The two factors are coprime, so

```text
[b(b-beta_n)]^m divides f.                            (14)
```

Then `c^(-s)f(b)=c^(nm-s)e^m h(b)` belongs to the right side of `(12)`.
Conversely every such expression is polynomial.  Nonnegative weights are
already `c^rC[b]`.  Weight decomposition is exact in both torus-stable
rings, proving `(12)` without a cancellation gap.

## 4. Geometry, Poisson structure, and two-bracket collapse

Let

```text
Y_n=Spec C[b,c,e]/(c^n e-Sigma_n(b)),
Sigma_n(b)=b(b-beta_n).                                (15)
```

Since `beta_n!=0`, `Sigma_n` is squarefree and `Y_n` is smooth.  Localizing
at `c` gives `A1_b x Gm_c`; the two boundary lines yield

```text
Pic(Y_n)=Z,              Y_n^*=C^*,              chi_c(Y_n)=2. (16)
```

In particular `Y_n` is not `A2`.  Give it the Jacobian Poisson bracket

```text
{b,c}=c^n,
{c,e}=-Sigma_n'(b),
{b,e}=-n c^(n-1)e.                                    (17)
```

This bracket is everywhere symplectic.  On `c!=0` its inverse form is

```text
omega_n=db wedge dc/c^n.                              (18)
```

More generally, for any squarefree `Sigma`, choose `U,T in C[b]` with
`U Sigma-T Sigma'=1`.  On `c^n e=Sigma(b)` set

```text
alpha_0=Uce db-nTe dc-Tc de.                          (19)
```

Then `c^(n-1)alpha_0=db`, so

```text
d(alpha_0/(n-1))=omega_n,                             (20)
```

and rewriting `alpha_0` gives the exact polynomial identity

```text
1={((U+T')/(n-1))ce,b}+{-Te,c}.                       (21)
```

Thus raising the boundary exponent does not restore the standard
Danielewski cohomological obstruction.  The residual debt is again the
compression of two brackets to one.

### 4.1 The homogeneous and two-by-two cells are universally empty

For the grading `wt(b,c,e)=(0,1,-n)`, every homogeneous piece has the form

```text
c^r f(b),       Sigma^ceil(-r/n) divides f for r<0,   (21a)
```

and

```text
{c^r f,c^s g}=c^(r+s+n-1)(s f'g-rfg').               (21b)
```

A homogeneous scalar bracket needs `r+s=1-n`.  If both weights are negative,
both coefficients vanish on every squarefree arm, so their bracket also
vanishes there and cannot be one.  Thus one weight is nonnegative and, after
swapping, take `r>=0`.  Root divisibility eliminates every case except
`r=1,s=-n`; there the
coefficient is

```text
-n f' Sigma g-f(Sigma g)',                            (21c)
```

whose degree is `deg f+deg g+deg Sigma-1>=1` with nonzero leading
multiplier.  Thus no homogeneous scalar bracket exists.

For the complete two-by-two cell, first subtract scalar weight-zero pieces.
Each extreme output weight has a unique bracket owner and must vanish.  At a
simple arm a vanishing homogeneous bracket cannot mix strict opposite signs;
a zero-weight boundary forces the retained zero-weight coefficient to be
scalar and is therefore removable.  If the upper extreme were still
negative, both scalar-row summands would vanish at every arm.  Consequently
the only surviving sign routing has the cross-matched normal form

```text
P=c^(-R)f+c^(T-n+1)F,
Q=c^(-T)g+c^(R-n+1)G.                                 (21d)
```

Arm survival forces `R=n` or `T=n`.  If `R=n`, the low and high extreme
Wronskians give, for `T=nk` and `p=n(k-1)+1`,

```text
f=A h,       g=B h^k,       Sigma|h,
F=L K^p,     G=M K.
```

The scalar row factors exactly as

```text
(K h'+n hK')
(AM-kp LB K^(p-1)h^(k-1)).                            (21e)
```

The first factor has positive degree and leading coefficient
`(deg h+n deg K)lc(h)lc(K)`.  The case `T=n` is symmetric, so the complete
two-by-two cell is empty.  More generally, if one output is homogeneous,
its unique complementary component in the other output would itself give a
forbidden homogeneous scalar bracket.  Hence every Darboux pair on `Y_n`
needs at least two pieces in each output and at least five nonconstant
homogeneous pieces in total.  The stronger exponent-two two-by-three
obstruction in fact has the following all-exponent form.

### 4.2 The complete two-by-three cell is universally empty

Suppose one output has two weights `r_0<r_1` and the other has three.
The scalar row needs at least two complementary pairs, because one pair
alone would be a forbidden homogeneous scalar bracket.  Put
`delta=r_1-r_0`.  The two complementary weights differ by `delta`.
Translation of supports in `(21b)` shows that the third weight must extend
this length-two arithmetic progression at its lower or upper end.  At any
other position its two bracket rows form a disconnected zero component and
can be deleted, reducing to the forbidden two-by-two cell.

The same simple-arm sign lemma used above, including its inert zero-weight
boundary, puts the two-weight output in the form

```text
P=c^(-R)f+c^(T-n+1)F,   R,T>=1,   R+T-n+1>0.           (21f)
```

The two possible three-weight supports are therefore

```text
LOWER: {-R-2T+n-1, -T, R-n+1},
UPPER: {-T, R-n+1, 2R+T-2n+2}.                         (21g)
```

Here the small-weight boundaries are load-bearing.  At a simple root of
`Sigma`, if `R<n-1` then both coefficients in the first scalar pair have
negative weight and vanish; if `R=n-1` its derivative multiplier is zero;
if `R>n` then `(21a)` makes even the derivative of `f` vanish.  The only
survivor is `R=n`, where `f` may have a simple zero and the opposite weight
is one.  The identical argument for the second scalar pair gives `T=n`.
Thus all cases `1<=R,T<n`, including zero-weight boundaries, are discharged
and arm survival is equivalent to

```text
R=n                         or T=n.                    (21h)
```

It remains to close the exceptional ladders forced by `(21h)`.
The unique high extreme forces, modulo the removable zero boundary,
`T>=n` when `R=n` and `R>=n` when `T=n`: below those thresholds its two
weights have strict opposite signs, while equality is the already separated
`R=T=n` case.  Thus the integer parameters used below omit no small ladder.

For the lower support write

```text
Q=c^(-R-2T+n-1)g+c^(-T)g_1+c^(R-n+1)H.                (21i)
```

If `R=n` and `T=n`, the two extreme Wronskians and the shared row give

```text
f=A h^n,                 g=B h^(2n+1),
F=L K,                   H=M K,
g_1=h^n(D+(2n+1)LB hK/(nA)),             Sigma|h.      (21j)
```

The scalar row is divisible by `h^(n-1)(hK)'`, so it cannot be one.  If
`R=n` but `T!=n`, survival of the first scalar summand forces
`gcd(n,2T+1)=n`.  Hence `n` is odd and, for some `k>=1`,

```text
2T+1=n(2k+1),             p=T-n+1,
f=A h,                    g=B h^(2k+1),
F=L K^p,                  H=M K,
g_1=(2k+1)LB h^(2k)K^p/A,                 Sigma|h.     (21k)
```

The shared first-order equation also permits a formal homogeneous summand
`v_0` with `n h v_0'-T h'v_0=0`.  It would satisfy
`h^T=lambda v_0^n`.  But `(2T+1)=n(2k+1)` gives `gcd(n,T)=1`, so UFD
comparison would make `h` an `n`th power, contradicting its simple zeros on
the squarefree arms.  Hence `v_0=0`, and `(21k)` is the complete polynomial
solution.

Its scalar row has the nonunit factor

```text
K h'+n hK'.                                             (21l)
```

Finally take `T=n,R>n` and put `d=gcd(R,n+1)`.  The lowest and highest
Wronskians give

```text
f=A h^(R/d),       g=B h^((R+n+1)/d),
F=L K,             H=M K^(R-n+1).
```

The complete shared-row solution is

```text
g_1=C h^((n+1)/d)K+v_0,
C=(R+n+1)LB/(RA),
v_0^d=lambda h^n.                                      (21m)
```

Here `gcd(d,n)=1`.  A nonzero polynomial `v_0` therefore forces
`h=mu J^d` and `v_0=nu J^n`; at every arm its order is at least `n`, so it
cannot supply the simple zero needed by the scalar row.  The particular
summand could be simple only if `d=n+1` and the arm order of `h` were one.
But then `ord(f)=R/(n+1)<ceil(R/n)`, contradicting `(21a)`.  Thus `g_1`
is never simple on an arm; the other scalar summand contains `f`, whose
order is at least two because `R>n`.  The scalar row vanishes at every arm,
so the lower support is empty.

For the upper support write

```text
Q=c^(-T)g+c^(R-n+1)G+c^(2R+T-2n+2)H.                  (21n)
```

If `R=n`, arm survival in the lowest Wronskian forces `T=nk`.  Put

```text
p=n(k-1)+1, q=nk+2, d=gcd(p,n+1),
m=p/d,             ell=(n+1)/d,             nu=q/d.
```

Then all extreme and shared rows have the exact solution

```text
f=A h,       g=B h^k,       F=L K^m,       H=M K^nu,
G=nu AM hK^ell/(mL)+G_0,                    Sigma|h,
K'G_0-dKG_0'=0.                                           (21o)
```

If `G_0=0`, the scalar row has factor `dKh'+nhK'`.  Otherwise UFD
comparison gives `K=lambda J^d,G_0=mu J`, and the factor is
`Jh'+nhJ'`.  Each has degree at least one, with nonzero characteristic-zero
leading coefficient.  The dual case `T=n` forces `R=nk` and has

```text
f=A h^k,       g=B h,       F=L K,       H=M K^q,
q=n(2k-1)+2,
G=D K^(n(k-1)+1)+(qAM/L)h^kK^(q-1),                   (21p)
```

whose scalar row again has factor `Kh'+nhK'`.

These cases exhaust `(21g)--(21h)`.  Therefore no constant bracket has
support sizes at most `(2,3)` or `(3,2)` on any smooth squarefree
`c^n e=Sigma(b)` with `n>=2` and `deg Sigma>=2`.  Every Darboux pair lies
in one of

```text
(2,>=4),                   (>=4,2),                   (>=3,>=3),
```

and uses at least six nonconstant homogeneous pieces.  Degree-one `Sigma` is
the sharp failure boundary: eliminating `b` makes the surface a polynomial
plane, and `{c,e}=-Sigma'` is already a homogeneous nonzero constant.

## 5. The polynomial etale collision map

The polynomial observables in `(8)` define

```text
Phi_n:A2_(x,q) -> Y_n.                                (22)
```

On the dense chart, `(10)` and `(17)` show that `Phi_n` is anti-Poisson:

```text
Jac_(x,q)(F(Phi_n),G(Phi_n))=-{F,G}(Phi_n).           (23)
```

Both sides are polynomial, so `(23)` holds everywhere.  Nondegeneracy of
the source and target brackets makes `Phi_n` etale.

The function-field degree is `deg B_n=d_n`; after a root `t` of
`B_n(t)=b` is chosen, `x=c/(P_nR)` and `q=t/x^n` are forced.  Therefore
`Phi_n` has generic degree `d_n`.

Its complex-point image is exact:

```text
Phi_n(A2(C))=Y_n(C) minus {(beta_n,0,0)}.              (24)
```

For `c!=0`, a good root of `B_n(t)=b` avoids `P_nR`: at `b=0` use `t=0`,
at `b=beta_n` use any root of `W_n`, and elsewhere every root is good.
The central arm `b=c=0` is filled using `t=0` and `W_n(0)=-beta_n`.
On the side arm `b=beta_n,c=0`, the root `t=-1` and `(5)` fill every
`e!=0`; `e=0` there would force `q=0` and contradict `t=-1`.

Every central-arm point with `e!=0` has exactly

```text
1+n(n-1)=d_n                                           (25)
```

distinct preimages: one at `x=0`, and `n` choices of `x` over each of the
`n-1` roots of `P_n`.  Hence `Phi_n` is noninjective.

## 6. Exact nonproperness and fibre invoice

In the coordinates `(t,x)`,

```text
b=B_n(t),          c=xP_n(t)R(t),          e=tW_n(t)/x^n. (26)
```

Because `B_n` is proper in `t`, an escaping source sequence with bounded
image has either `x->infinity` and `P_nR->0`, or `x->0` and `W_n->0`.
The roots are simple and disjoint, so every limiting point is realized.  The
Jelonek set is exactly

```text
N(Phi_n)={b=0,e=0} union {b=beta_n,e=0}
         union {b=beta_n,c=0}.                         (27)
```

The fibre sizes expose the lost sheets:

```text
c!=0:             d_n / 1 / (n-1)^2
                   at generic / b=0 / b=beta_n;
central arm:      d_n for e!=0, 1 at its origin;
side arm:         n for e!=0, 0 at the omitted origin. (28)
```

Etaleness removes ramification from the surface map, but not the nonproper
boundary modes.

## 7. Alternating/symmetric monodromy and the degree-seven carrier

The exact Shabat identification `(7a)` and the cited monodromy theorem `(7b)`
show that the inverse carrier is already the full alternating group for every
odd `n` and the full symmetric group for every even `n`.  This sharply
separates the tower from a radical or hidden Kummer construction.

For `n=3`,

```text
P_3=(2t^2+7t+14)/14,
beta_3=-729/2744,
W_3=(8t^4+60t^3+258t^2+557t+729)/2744.               (29)
```

The degree-seven cover has passports

```text
over 0:             (3,3,1),
over beta_3:        (3,1,1,1,1),
over infinity:      (7).                              (30)
```

As an independent check of the cited all-scale result, its monodromy group is
exactly `A_7`: transitivity in prime degree makes the group primitive; side
inertia is a literal 3-cycle, so Jordan's theorem puts `A_7` inside the
group.  Every branch permutation in `(30)` is even, giving the reverse
inclusion.  This is the smallest nonsolvable carrier in the tower, unlike
the dihedral/Kummer inverse of the exponent-two Chebyshev family.

## 8. Counterexample reduction, ledger, and boundary

If `P,Q in C[Y_n]` satisfied

```text
{P,Q}=1,                                               (31)
```

then `(P(Phi_n),Q(Phi_n))` would be a polynomial etale noninjective
endomorphism of `A2`, hence a counterexample to `JC(2)`.  Equations `(12)`
and `(21)` sharpen the search: the natural polynomial-observable ring is
already maximal, and the constant is already a sum of two brackets.  What is
missing is one polynomial Darboux pair on the non-`A2` target.

```text
source:             the rational pair (a,c) on A2,
target:             Y_n and its maximal observable ring C[b,c,e],
map:                critical-value completion (8),(22),
preserved:          constant Jacobian, symplectic bracket, collision degree,
destroyed:          proper sheet data at P_nR and W_n,
needed sidecar:     a one-bracket Darboux compression on Y_n,
cheapest hostile:   n=3, where even A7 monodromy still leaves a non-A2 target.
```

The theorem makes no literature-novelty claim, no assertion that `(31)` is
solvable, and no conclusion about `JC(2)`.  Repeated-root `Sigma`, exponent
one, positive characteristic, and non-weight-equivariant completions are
outside its scope.

## 9. Exact verification contract

The companion uses exact rational symbolic arithmetic.  For `2<=n<=8` it
checks the polynomial/ODE and critical-value identities, the exact Shabat
change of variables, structural maximal-intersection exponent invoices,
Jacobian and surface identities, the two-bracket formula, and complete
two-by-two and two-by-three factorization controls.  It independently checks every displayed
`n=3` coefficient, all Poisson-generator pullbacks, and the
discriminant/passport.  These finite rows control the formulas; the image,
fibre, nonproperness, and universal quantifiers come from the explicit
valuation and torus-grading proofs above rather than finite extrapolation.
