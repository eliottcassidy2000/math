---
id: THM-2796
title: "Balanced response Stieltjes--Pade normal form and complete one-double-zero chamber"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, with a
  FINITE-EXACT dessin census through degree eight.  Every balanced
  square-potential response has an exact constant-numerator logarithmic
  derivative, a Stieltjes first integral, moment annihilation through the
  penultimate index, and a forced first square defect.  The full e<=1
  passport chamber is classified in all degrees.  This is a response-layer
  theorem, not chart entry, degree-at-least-26 closure, JC(2), or DC(2).
source: root/balanced-response-stieltjes-pade-2026-07-28
depends_on:
  - THM-2784-nonsplit-response-square-potential-divisor-and-infinity-classification
related:
  - THM-2245-degree-fourteen-spectral-quartic-discriminant-reduction
  - THM-2247-nonsplit-terminal-quartic-degree-fourteen-closure
  - THM-2760-exact-prefix-even-faber-flux-gcd-and-smooth-boundary-exclusion
  - THM-2778-all-degree-complete-chosen-sheet-split-exact-prefix-closure
script: 04-computation/jc2_balanced_response_stieltjes_pade_classification_thm2796.py
output: 05-knowledge/results/jc2_balanced_response_stieltjes_pade_classification_thm2796.out
script_sha256: 30b3a490d4529a1394ff770c3d9e49ffca4de2340f2b7533e5d3c1d121406a48
output_sha256: 0c7023ae6f89fd494b3514402cc553d330671c94fafc2947218db8eb50507f85
independent_script: 04-computation/jc2_balanced_response_stieltjes_pade_hostile_audit_thm2796.py
independent_output: 05-knowledge/results/jc2_balanced_response_stieltjes_pade_hostile_audit_thm2796.out
independent_script_sha256: bc1dd327b56aef6c3e5cc2bf421a3e4193dc8242888d7b5fec75200445375ede
independent_output_sha256: 3cccaeb56d3458d1182b802fcf0c1f0b79f76202398a5e3ea53e0c41368fc6f0
hash_basis: LF-normalized bytes
---

# THM-2796 -- balanced response Stieltjes--Padé normal form and the complete one-double-zero chamber

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The dessin
census is **FINITE-EXACT through `N=8`**; the factor, moment, square-contact,
and `e<=1` classification theorems are all-degree proofs.

This theorem sharpens only the balanced response branch of THM-2784.  It
does not close the nonsplit quartic chart, prove chart entry, or identify the
source polynomial of an arbitrary Keller map.

## 1. Strongest result

Work over an algebraically closed field of characteristic zero.  Fix
`kappa!=0` and a balanced nonconstant solution

```text
F=V G^2,                       2VG'+V'G=2kappa,
```

from THM-2784.  Thus `F` has equally many finite zeros and finite poles,
say degree `N`, and tends to a nonzero value `mu` at infinity.  Use the
balanced passport notation

```text
zero partition       2^e 1^s,
pole partition       p=(p_1,...,p_h),
third partition      (r,1^(N-r)),

N=s+2e=sum_j p_j,
r=s+e+h-1=N-e+h-1,
N-r=e-h+1.
```

Let the distinct geometric points be

```text
alpha_i : simple zeros of F and simple roots of V,
gamma_l : double zeros of F away from V=0,
beta_j  : poles of F and high roots of V.
```

Put

```text
S=product_i(x-alpha_i),        deg(S)=s,
E=product_l(x-gamma_l),        deg(E)=e,
T=product_j(x-beta_j),         deg(T)=h,
D=product_j(x-beta_j)^(p_j),
T_j=T/(x-beta_j),
M=S E T,                       deg(M)=r+1.
```

All four displayed polynomials are monic.  Then there are nonzero constants
`v,g`, with

```text
lambda=vg,                     mu=vg^2,
```

such that the entire solution is

```text
V=v S D T^2,
G=g E/(D T),
F=mu S E^2/D.                                      (1)
```

The response ODE is **equivalent** to the single polynomial identity

```text
2 S T E' + E T S' - E S sum_j p_j T_j = C,          (2)
C=2kappa/lambda in C*.                              (3)
```

Conversely, **within the displayed balanced degree constraints**
`N=sum_j p_j=s+2e` (so in this nonconstant chamber `h>=1` and `r>=1`),
pairwise-disjoint squarefree `S,E,T`, positive `p_j`, and a nonzero constant
identity `(2)` construct a balanced square-potential solution by choosing
`v,g` with `vg=2kappa/C` and using `(1)`.  The potential is genuinely
nonsplit precisely when

```text
s>0, or at least one p_j is odd.                     (4)
```

Thus `(2)` is a complete factor-data classification, not merely a necessary
condition.

### Proof of the factor normal form

The valuation ledger in THM-2784 gives the factors of `V` and `F` in `(1)`.
At an `alpha_i`, `F'` is a unit; at a `gamma_l`, `F'` has a simple zero; and
at a pole `beta_j` of order `p_j`, `F'` has pole order `p_j+1`.  At every
other finite point, the square-potential equation says `F'` is a unit.
Therefore `G=F'/(2kappa)` has exactly the divisor asserted in `(1)`.
The identity `F=VG^2` gives `mu=vg^2`.  Substitution into
`2VG'+V'G=2kappa` gives `(2)`--`(3)`.  The same calculation backwards proves
the converse.  Finally, the squareclass of `V` is represented by
`S product_j(x-beta_j)^(p_j)`, proving `(4)`.

## 2. Exact Padé identity and moment annihilation

Define the signed integral residue measure

```text
nu = sum_i delta_(alpha_i)
     +2 sum_l delta_(gamma_l)
     -sum_j p_j delta_(beta_j).
```

Its Cauchy transform is the logarithmic derivative of `F`.  Identity `(2)`
is exactly

```text
F'/F
 = S'/S + 2E'/E - sum_j p_j/(x-beta_j)
 = C/M.                                               (5)
```

This is an **exact rational identity**, not a Padé analogy.  The initially
degree-`r` numerator of the left side cancels all the way to a nonzero
constant.  Equivalently, if

```text
m_k = sum_i alpha_i^k
      +2 sum_l gamma_l^k
      -sum_j p_j beta_j^k,
```

then expansion at infinity gives

```text
m_0=m_1=...=m_(r-1)=0,            m_r=C!=0.          (6)
```

Indeed `m_0=s+2e-N=0` is the balanced degree identity.  If `q` is the first
index with `m_q!=0`, then the monic degree-`r+1` polynomial `M` makes
`M F'/F` have leading term `m_q x^(r-q)`.  It is a nonzero constant exactly
when `q=r`, proving `(6)`.

In Padé language, `(5)` is the exact `[0/(r+1)]` constant-numerator
representation of the logarithmic derivative at infinity, and `(6)` is its
order-`r` cancellation condition.  No convergence or asymptotic theorem is
being used.

The cheapest coefficient test for any proposed balanced response `F_*` is
therefore:

```text
find monic M of degree r+1 and C in C*
such that                    M F_*' - C F_* = 0.      (7)
```

After denominators are cleared, `(7)` is linear in the coefficients of `M`
and in `C`.  It avoids root extraction, dessin enumeration, and construction
of `U=sqrt(V)`.

## 3. Exact Heine--Stieltjes-type operator and Bethe equations

For fixed `S,T,p`, set

```text
B=T S' - S sum_j p_j T_j,
L_1(E)=2ST E'+B E.
```

Then `(2)` is the inhomogeneous first-order equation

```text
L_1(E)=C.                                             (8)
```

Differentiating gives the exact homogeneous second-order operator

```text
L_2(E)
 =2ST E''
  +[2(ST)'+B]E'
  +B'E
 =0.                                                  (9)
```

The coefficient degrees in `(9)` have the Heine--Stieltjes pattern, with
locked Van Vleck term `-B'`.  This is an exact operator identification:
`L_2=d/dx composed with L_1`.  It is more constrained than the generic
Heine--Stieltjes problem because the accessory polynomial is not free, and
the required solution must also have the nonzero first integral `(8)`.
No external Heine--Stieltjes classification is being invoked here.

At a root `gamma_l` of `E`, equation `(9)` yields the exact Bethe equations

```text
2 sum_(m!=l) 1/(gamma_l-gamma_m)
+(3/2) sum_i 1/(gamma_l-alpha_i)
+sum_j (1-p_j/2)/(gamma_l-beta_j)
=0.                                                   (10)
```

They are a differentiated shadow of the stronger common-value
interpolation laws obtained directly from `(2)`:

```text
lambda E(alpha_i)T(alpha_i)S'(alpha_i) = 2kappa,
lambda S(gamma_l)T(gamma_l)E'(gamma_l) = kappa,
-lambda p_j E(beta_j)S(beta_j)T_j(beta_j) = 2kappa.  (11)
```

Thus every marked root is tied to the same nonzero scale.  Equations
`(10)`--`(11)` are the most economical root/residue constraints to add after
the two Faber flux equations have specified candidate coefficient
functions.

## 4. Exact square contact at infinity

The factorization has a second consequence which is invisible in the bare
passport.  From `(1)`,

```text
V = v M^2 (mu/F).                                    (12)
```

Since `(5)` gives

```text
F/mu = 1-(C/r)x^(-r)+O(x^(-r-1)),
```

one obtains

```text
V/(vM^2)=1+(C/r)x^(-r)+O(x^(-r-1)),                 (13)

deg(V-vM^2)=r+2,
LC(V-vM^2)=vC/r=2kappa v/(r lambda) != 0.            (14)
```

Both `V` and `vM^2` have degree `2r+2`.  Hence their coefficients from
degrees `2r+2` down through `r+3` agree: exactly the top `r` coefficients.
The first defect, in degree `r+2`, is prescribed and nonzero.

This interfaces exactly, but on a different axis, with the repository's
polynomial exact-square-prefix chart.  The inherited chart has
`P=H^2+L` as a polynomial identity in the fibre variable `z`; `(14)` gives
an order-`r` square contact for the base coefficient `V` in the source
variable `x`.  It does **not** turn `V` into a square and does not move the
branch to the chosen-sheet split chart.  In particular:

```text
original exact prefix:       exact in z;
new response square contact: truncated in x, with a forced nonzero defect.
```

Any independent flux specialization forcing
`deg(V-vM^2)<=r+1` kills the balanced response immediately.  More generally,
`V=vM^2+W`, `deg(W)=r+2`, `LC(W)=2kappa v/(r lambda)` is a coefficient
parameterization that removes the top `r` variables before a flux
elimination.

## 5. Complete all-degree classification for `e<=1`

Since `h<=e+1`, the cases below exhaust every balanced passport having at
most one extra double zero.

### 5.1 `e=0`: the cyclic hostile is the whole chamber

Here necessarily

```text
h=1, s=N, p=(N), r=N.
```

Translate the high root to zero and write `y=x-beta`.  Equation `(2)` is

```text
yS'-NS=C.
```

For monic `S`, every intermediate coefficient must vanish, so

```text
S=y^N+a,                   a!=0, C=-Na.              (15)
```

Consequently every map in this chamber is affine/target equivalent to

```text
F_N=4(1-y^(-N)),
V_N=y^(N+2)(y^N-1)/N^2,         kappa=1.             (16)
```

The map is a function of `y^N`; its monodromy is cyclic `C_N`.  Thus the
hostile family already installed in THM-2784 exhausts the natural `e=0`
subchamber, rather than representing one isolated construction.

### 5.2 `e=h=1`: one map, and it is symmetric-monodromy

Now `N>=2`,

```text
s=N-2, p=(N), r=N-1.
```

Normalize the pole to `0`, the double zero to `1`, and `F(infinity)=1`.
Since

```text
G=g(x-1)/x^(N+1),
```

the primitive and the double-zero value uniquely force

```text
B_N(x)=x^N-Nx+N-1,
Q_N(x)=B_N(x)/(x-1)^2,

F_N(x)=B_N(x)/x^N,
V_N(x)=4x^(N+2)Q_N(x)/(N^2(N-1)^2),      kappa=1.   (17)
```

The only common root of `B_N` and `B_N'=N(x^(N-1)-1)` is `x=1`, where
`B_N''(1)!=0`.  Hence `Q_N` is squarefree and disjoint from `0,1`.
For `N=2`, `Q_N=1` and `V_N=x^4` is split.  For every `N>=3`, `V_N` is
genuinely nonsplit.

The branch cycles are an `N`-cycle and an adjacent transposition, so the
monodromy is `S_N`.  The smallest genuinely noncyclic member is `N=3`:

```text
F=(x-1)^2(x+2)/x^3,
V=x^5(x+2)/9,

passport:      (2,1), (3), (2,1),
monodromy:     S_3.                                  (18)
```

No transitive group in degree one or two is noncyclic, so `(18)` is
globally the smallest genuinely noncyclic balanced response.

### 5.3 `e=1,h=2`: explicit reciprocal-monomial chord family

Write the unordered pole partition as `(d,b)`, with

```text
d+b=N,                       1<=d<=floor(N/2).
```

Normalize the two poles to `0,1` and put

```text
D_(N,d)(x)=x^d(x-1)^b,
gamma=d/N,
c=D_(N,d)(gamma)=(-1)^b d^d b^b/N^N.                (19)
```

The residues of

```text
(x-gamma) dx/[x^(d+1)(x-1)^(b+1)]
```

at both poles vanish if and only if `gamma=d/N`.  Equivalently,

```text
d[-1/(N D_(N,d))]/dx
 =(x-d/N)/[x^(d+1)(x-1)^(b+1)].                     (20)
```

The unique normalized map and potential are therefore

```text
F_(N,d)=1-c/D_(N,d),
Q_(N,d)=(D_(N,d)-c)/(x-d/N)^2,

V_(N,d)
 =4x^(d+2)(x-1)^(b+2)Q_(N,d)/(c^2N^2),   kappa=1.  (21)
```

The only off-pole critical point of `D_(N,d)` is `d/N`, and it is simple.
Thus `D_(N,d)-c` has an exact double root there and every other root is
simple.  This proves directly that `Q_(N,d)` is squarefree.  Every member of
`(21)` is genuinely nonsplit.

Combinatorially, the zero transposition is a chord of cyclic distance `d`
in the `N`-cycle supplied by the third branch point.  Distances `d` and
`N-d` are conjugate, giving exactly `floor(N/2)` maps.  If

```text
g=gcd(N,d),                      m=N/g,
```

the conjugates of the chord generate `S_m^g` on the `g` residue classes,
while the `N`-cycle rotates those classes.  Hence the monodromy order is

```text
g (m!)^g.                                            (22)
```

It is `S_N` exactly when `g=1`; otherwise it is imprimitive.  This gives a
complete all-degree classification, explicit maps, and exact monodromy for
the second `e=1` chamber.

## 6. Exact bounded census beyond the classified chamber

An independent permutation-triple enumeration fixes the third inertia as
the standard `r`-cycle, enumerates involutions with `e` transpositions,
tests transitivity, and quotients by the full centralizer of the third
inertia.  Through degree eight it gives:

| `N` | all dessins | genuine nonsplit | distinct genuine passports | genuinely noncyclic |
|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 0 |
| 2 | 3 | 2 | 2 | 0 |
| 3 | 3 | 3 | 3 | 2 |
| 4 | 7 | 6 | 6 | 5 |
| 5 | 10 | 10 | 9 | 9 |
| 6 | 26 | 25 | 18 | 24 |
| 7 | 49 | 49 | 24 | 48 |
| 8 | 127 | 126 | 44 | 125 |

The first failure of “passport determines map” occurs at `N=5`: the
passport

```text
e=2, h=2, p=(4,1)
```

has two nonisomorphic dessins (both monodromy order `20`).  This is the
sharp stopping boundary for a passport-only classification.  From `e=2`
onward one must retain a Nielsen/dessin coordinate or solve the Stieltjes
system `(2)`; the integer partition alone discards real moduli-free but
combinatorially distinct maps.

The census is finite exact evidence only through `N=8`.  The factor theorem,
moment theorem, square-contact theorem, and complete `e<=1` classification
are all-degree proofs.

## 7. Relation to THM-2760's Jacobi/simple-root carrier

There is a useful common meta-pattern but no proved transport map.

THM-2760 uses the Euler operator in the Faber residue parameter to show
that a common zero of two residual flux polynomials is exactly a repeated
root of one carrier polynomial.  Its Schur--Szego/Jacobi factorization then
proves that carrier squarefree.  Here, `L_2=(d/dx)L_1` turns the balanced
response into a Stieltjes/Bethe root system, while the nonzero first integral
`L_1(E)=C` supplies the common-value and simple-root carrier.  In both
arguments, a simultaneous condition is converted into a discriminant or
root-collision condition.

The correspondence stops there:

```text
THM-2760 variable/operator:    Faber residue r and theta=r d/dr;
present variable/operator:     source x and d/dx;

THM-2760 squareclass:          chosen-sheet split;
present squareclass:           genuinely nonsplit;

THM-2760 root theorem:         Jacobi + Schur--Szego;
present root theorem:          nonzero first integral, plus explicit
                               reciprocal-monomial geometry for e<=1.
```

No substitution from `E,S,T` to THM-2760's Jacobi polynomial `H_k` has been
proved, and THM-2760's simple-negative-root theorem cannot be imported into
this nonsplit chamber.  The exact new maps are `(5)`, `(8)`--`(9)`, and
`(12)`--`(14)`; the link to THM-2760 is presently methodological analogy.
A plausible next niche is to ask whether fixed low-`h`, higher-`e`
instances of `(9)` admit hypergeometric/Jacobi realizations, but that is
open.

## 8. The Faber quotient already is the Padé denominator

In the nonsplit quartic notation, write

```text
q^2=T_F=A_src^2/V,                R_Q=q K.
```

The response square function is then **exactly**

```text
F_resp=R_Q^2=T_F K^2.                              (23)
```

This is the concrete map from the Faber carrier to the present theorem.  It
is not an analogy.

There is a sharper identity than an auxiliary linear search for `M`.  Since
`q=A_src/U`, `R_Q=UG`, and `U^2=V`,

```text
A_src K=U R_Q=V G=lambda M.                          (24)
```

The retained Keller one-form then gives

```text
(T_F K^2)'/(T_F K^2)=2kappa/(A_src K)=C/M,
V T_F K^2=(A_src K)^2=lambda^2 M^2.                 (25)
```

Thus the direct response carrier is `P_resp=A_src K`.  In every balanced
survivor it is a nonzero scalar multiple of the squarefree polynomial
`M=SET`, has degree `r+1`, and its marked roots carry logarithmic residues
`(+1,+2,-p_j)`.  The cheapest test is therefore:

```text
derive K=R_Q/q;
form P_resp=A_src K;
test polynomiality, degree r+1, squarefreeness and support factorization;
test the residue alphabet and balance;
test the forced degree-(r+2) square defect (14);
only then normalize the remaining spectral curve.
```

THM-2245's degree-fourteen formula is a retrospective exact control for
`(24)--(25)`, not a live target: THM-2247 already closes degree fourteen,
and later canon closes the inherited nonsplit degree-twenty-two branch.  The
live program begins at degree twenty six: derive its `K`, then apply the
`P_resp` tests above before another spectral normalization.  That program is
**OPEN**.  It does not by itself prove that an abstract spectral solution
comes from polynomial `A_src,V`, nor does it supply the original fibre exact
prefix.

## 9. Exact controls

Primary artifacts:

```text
04-computation/jc2_balanced_response_stieltjes_pade_classification_thm2796.py
05-knowledge/results/jc2_balanced_response_stieltjes_pade_classification_thm2796.out

04-computation/jc2_balanced_response_stieltjes_pade_hostile_audit_thm2796.py
05-knowledge/results/jc2_balanced_response_stieltjes_pade_hostile_audit_thm2796.out
```

The script has no Python `assert` node.  It performs:

- the full permutation-triple census through `N=8`;
- exact group closure and genuine cyclicity by element order;
- the all-degree chamber claims in bounded combinatorial controls;
- symbolic cyclic and one-pole families through `N=10`;
- all `25` two-pole chord maps through `N=10`;
- the factor ODE, square-potential equation, exact Padé identity,
  differentiated Heine--Stieltjes operator, Newton-sum moments,
  interpolation remainders, nonsquare boundary, and exact first square
  defect; and
- normal/optimized/stored transcript identity.

The symbolic moment checks use Newton identities, not radical root
extraction, so irreducible roots are not silently omitted.  Group cyclicity
is tested by an element of full group order, not by the weaker condition of
abelianness.

Reproduce with

```bash
python3 04-computation/jc2_balanced_response_stieltjes_pade_classification_thm2796.py
python3 -O 04-computation/jc2_balanced_response_stieltjes_pade_classification_thm2796.py
python3 04-computation/jc2_balanced_response_stieltjes_pade_hostile_audit_thm2796.py
python3 -O 04-computation/jc2_balanced_response_stieltjes_pade_hostile_audit_thm2796.py
```

The primary transcript reports `1485` exact gates and `FAILED CHECKS: NONE`.
The independent engine rederives the factor identity, hostile missing-balance
witness, moment/defect indexing, all three `e<=1` families, chord monodromy,
the first duplicate passport, and `(24)--(25)`.  Both companions have zero
Python `assert` nodes; ordinary, optimized, and stored transcripts agree.

QED.
