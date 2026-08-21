---
id: THM-3576
title: "Higher-exponent Belyi Keller collision tower"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; PENDING INDEPENDENT AUDIT.
  For every integer n at least two, an explicit rational Keller pair has a
  maximal polynomial-observable completion on the smooth surface
  c^n e=b(b-beta_n).  The resulting polynomial map from A2 is etale,
  generically of degree n(n-1)+1, noninjective, and misses exactly one target
  point.  Its physical symplectic form is exact and the constant function is
  a sum of two polynomial Poisson brackets, while every homogeneous and
  two-by-two weight cell is empty.  The n=3 member has degree seven and
  alternating monodromy A7.  A single polynomial Darboux pair on any
  target would give a planar Jacobian counterexample; no such pair and no
  counterexample to JC(2) is claimed.
source: root, 2026-08-20
audit: PENDING INDEPENDENT HOSTILE AUDIT
depends_on: []
related:
  - THM-3566-chebyshev-pell-odd-keller-collision-tower
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
script: 04-computation/jc2_higher_exponent_belyi_keller_thm3576.py
output: 05-knowledge/results/jc2_higher_exponent_belyi_keller_thm3576.out
script_sha256: PENDING
output_sha256: PENDING
hash_basis: LF-normalized bytes
---

# THM-3576 -- higher-exponent Belyi Keller collision tower

**PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; PENDING INDEPENDENT AUDIT.**
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

A homogeneous scalar bracket needs `r+s=1-n`.  After swapping, take `r>=0`.
Root divisibility eliminates every case except `r=1,s=-n`; there the
coefficient is

```text
-n f' Sigma g-f(Sigma g)',                            (21c)
```

whose degree is `deg f+deg g+deg Sigma-1>=1` with nonzero leading
multiplier.  Thus no homogeneous scalar bracket exists.

The complete two-by-two cell has cross-matched normal form

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
`(deg h+n deg K)lc(h)lc(K)`.  The case `T=n` is symmetric.  Hence every
putative Darboux pair on `Y_n` needs at least five nonconstant homogeneous
pieces.  The stronger exponent-two two-by-three obstruction does not
automatically transfer to `n>2`.

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

## 7. The degree-seven alternating carrier

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

Its monodromy group is exactly `A_7`.  Transitivity in prime degree makes
the group primitive; side inertia is a literal 3-cycle, so Jordan's theorem
puts `A_7` inside the group.  Every branch permutation in `(30)` is even,
giving the reverse inclusion.  This is a nonsolvable carrier, unlike the
dihedral/Kummer inverse of the exponent-two Chebyshev tower.

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

The companion uses exact rational symbolic arithmetic.  It checks the
universal polynomial/ODE, divisibility, Jacobian, surface, critical-value,
maximal-intersection, and fibre invoices for `2<=n<=8`.  It independently
checks every displayed `n=3` coefficient, all Poisson-generator pullbacks,
the discriminant/passport, seven collision branches, and the three
nonproperness modes.  These finite rows control the formulas; the universal
quantifiers come from the coefficient, valuation, and torus-grading proofs
above.
