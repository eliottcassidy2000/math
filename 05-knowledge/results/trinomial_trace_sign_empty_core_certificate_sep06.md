# Carry-corrected trinomial signs as a Hermite trace-form test

**Status: PROVED reformulation using a classical mechanism; REFUTED
uncorrected canonical negativity; FINITE-EXACT named controls.** The new
application is to the complete first and doubled trinomial factorial rows,
with the carry factor retained. General negative definiteness, general
higher-channel two-rung coprimality, and the sharp arbitrary-support return
bound remain **OPEN**.

## 1. Inheritance and provenance

The source is
[THM-4436 — complete factorial-row simple negative roots](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md).
The inherited carry hostile is `(-13,1,8)`; the row `(-15,1,9)` below refines
it to three first channels. The corrected near miss is inferring a cross-mass
sign from individual real-rootedness;
the [companion signed-remainder note](trinomial_root_sign_empty_core_returns_sep06.md)
gives both a logical counterexample and the actual carry correction. Its
old adjacent-pair controls, determinant zero at `(-4,3,10)` and positive at
`(-4,5,14)`, concern failed pair compression; the complete rows still pass.

Targeted repository searches for `Hermite trace`, `Tr_A`, `trace form` and
`trace-form` recovered
[THM-2842 — positive-cone Vandermonde observability, Section 3](../../01-canon/theorems/THM-2842-ordered-positive-cone-vandermonde-multiplier-observability.md)
(weighted quotient traces at simple real nodes), and
[THM-3519 — three-coordinate discriminant class, Section 7](../../01-canon/theorems/THM-3519-level-three-sporadic-keller-three-coordinate-primitivity-and-common-discriminant-class.md)
(trace Gram matrices and basis congruence). These are prior operations;
the abstract trace-form identity is not new.

For explicit external provenance, Le–Safey El Din,
[*Solving parametric systems of polynomial equations over the reals through
Hermite matrices*, arXiv:2011.14136v2, Section 4.1, Definition 5](https://arxiv.org/pdf/2011.14136v2#page=9),
defines the quotient multiplication-trace construction. Nathanson,
[*The Hermite–Sylvester criterion for real-rooted polynomials*,
arXiv:1911.01745v2, Theorem 1 and its proof](https://arxiv.org/pdf/1911.01745v2#page=2),
gives the root-power matrix and sum-of-squares argument. Both primary PDFs
were read. The weighted factorization needed here is derived directly below;
no general Hodge–Riemann, Turán or total-positivity implication is imported.

Concurrent incoming work at `origin/main` commit `e5a45df05f`,
[root geometry, Section 6](overnight3_20260906_moments_root_geometry.md),
already retains the same raw carry factor and records a **FINITE-EXACT**
bank of 221 supports and 1015 first-root signs. This companion credits that
independent recovery; its application is the rational trace certificate,
not a newly claimed general raw-sign theorem or a larger sign census.

## 2. The exact carry correction

Let the primitive support be `(-a,b,c)`, `a>=1`, `1<=b<c`, and suppose its
first support return has at least three channels. In the canonical notation
of THM-4436 and the
[proved semigroup classification](synthesis_20260905_moments_trinomial.md), put

```text
g=gcd(a+b,a+c), A=(a+b)/g, B=(a+c)/g, delta=B-A,
v=(x,Bh+r,z), 0<=r<B, 0<=z<A, h>=2,
d=(delta,-B,A), tau=alpha^delta*gamma^A/beta^B,
M=alpha^x*beta^(Bh+r)*gamma^z.
```

The complete first polynomial `P=P_g` has degree `h`, positive integer
coefficients and `CT(f^g)=M P(tau)`. For the complete doubled polynomial
`Q=P_(2g)` with its own canonical left anchor, define

```text
epsilon_r=floor(2r/B), epsilon_z=floor(2z/A),
h2=2h+epsilon_r+epsilon_z,
r2=2r-B*epsilon_r, z2=2z-A*epsilon_z,
x2=2x-delta*epsilon_z.
```

Thus `deg Q=h2` and the second anchor is `v2=2v-epsilon_z*d`. In particular

```text
CT(f^(2g))/M² = tau^(-epsilon_z) Q(tau).              (1)
```

The actual positive Laurent row supplies the necessary admissibility:
`b=A*x-delta*z>0`, hence `x>delta*z/A`. Therefore
`x2=2x-delta*floor(2z/A)>0`. The lower carry cannot produce a negative
count here. An arbitrary abstract factorial row from THM-4436 need not
satisfy this inequality; its doubled tuple requires a separate check.

The lower carry changes the sign on the negative ray; the upper carry also
changes the degree and coefficients and must remain in the full row.
Changing the first anchor to `v+k*d` multiplies the quotient in `(1)` by
`tau^(-2k)`, which is positive at every negative root. Its sign is therefore
independent of the first-channel anchor. Canonical `Q` without the factor
in `(1)` does not have that property across supports.

For the explicit three-channel support `(-15,1,9)`, `g=8`, the first tuple
`(A,B,h,r,z,x)` is `(2,3,2,0,1,1)`, so `epsilon_z=1`. Exactly

```text
P=56+560tau+56tau²,
Q=16+10920tau+400400tau²+1681680tau³+720720tau⁴+8008tau⁵,
Q mod P=-47087024-466126752tau.
```

The remainder zero `rho=-2942939/29132922` is right of the vertex `-5` of
`P/56`, and `(P/56)(rho)=23910838225/848727144258084>0`. It is right of
both roots. Its negative slope therefore makes `Q` positive at both roots,
refuting uncorrected canonical negativity. In contrast,

```text
tau^(-1)Q mod P=47087024tau+4743488
```

is negative at both. This control is not claimed minimal, and it does not
refute the corrected candidate.

## 3. Exact rational certificate for the corrected predicate

Divide `P` by its positive leading coefficient and form the rational algebra

```text
K=Q[tau]/(P),               dim_Q K=h.
```

Since `P(0)>0`, the class of `tau` is invertible. Thus the unique
degree-below-`h` residue

```text
R=tau^(-epsilon_z)*Q mod P
```

exists over `Q`, computed by polynomial inversion and remainder. Let `Tr_K`
mean the trace of multiplication on this finite algebra. Define the
symmetric rational matrix, indexed by `0<=i,j<h`,

```text
H_R[i,j]=Tr_K(R*tau^(i+j)).                            (2)
```

THM-4436 gives distinct negative real roots `lambda_1,...,lambda_h` of `P`.
After extending scalars to `R`, evaluation identifies `K` with `R^h`.
For `V[l,j]=lambda_l^j`, the literal matrix identity is

```text
H_R=V^T diag(R(lambda_1),...,R(lambda_h)) V.            (3)
```

The Vandermonde matrix is invertible. Consequently

```text
tau^(-epsilon_z)Q(tau)<0 at every P-root
  iff H_R is negative definite
  iff (-1)^k det(H_R[0:k,0:k])>0 for k=1,...,h.       (4)
```

The last equivalence is the real symmetric leading-principal-minor
criterion, obtained by successive completion of squares with pivot ratios.
All tests in `(4)` are rational; no floating-point root approximation is
needed. A companion multiplication matrix or Newton sums computes `(2)`.

The weaker noncancellation predicate is merely

```text
gcd(P,Q)=1 iff det(H_R)!=0.                            (5)
```

Indeed `(3)` makes its determinant the positive square `det(V)^2` times
the product of the nonzero root values. Determinant nonvanishing forgets
the individual signs. An anchor change gives
`H_(tau^(-2k)R)=D^T H_R D`, where `D` is multiplication by `tau^(-k)`;
it therefore preserves inertia directly over the quotient algebra.

This is the connection contract: complete factorial rows and their two
carries map to a weighted trace Gram matrix. The map preserves each
rootwise corrected sign and exact common-root failure. Passing to just
the determinant loses sign information; passing to individual polynomial
real-rootedness loses the entire matrix `H_R`. The required sidecar is the
full weighted form or its complete signed leading minors.

## 4. Named exact controls and the uniform quadratic model

The [standalone producer](../../04-computation/trinomial_trace_sign_empty_core_certificate_sep06.py)
and [output](trinomial_trace_sign_empty_core_certificate_sep06.out) retain
three named supports, not a census. They enumerate both full channel rows,
check literal Laurent multiplication and every earlier empty moment, and
compare multiplication-matrix traces with an independent Newton recursion.
The signed leading minors in `(4)` are:

```text
(-15,1,9), g=8, (epsilon_r,epsilon_z)=(0,1):
  461383264,
  587632760217600.

(-21,1,12), g=11, (epsilon_r,epsilon_z)=(0,1):
  752350432547692,
  21236578848830718128306804/3,
  942083106929885721679784497035400/27.

(-25,1,14), g=13, (epsilon_r,epsilon_z)=(1,1):
  248262381336934549/309375,
  103041355107681072187206907429/2994750000,
  464540893300260027432902715300343770182359/1452901465125000000.
```

All are positive: the corrected second value is strictly negative at every
first root for these three controls. This is **FINITE-EXACT**, including
the two-carry cubic case, not a general negative-definiteness theorem.

The [companion width-four proof](trinomial_root_sign_empty_core_returns_sep06.md)
proves the same predicate uniformly for `(-4,g-4,2g-4)`, odd `g>=5`, using
its remainder `C_g+D_g*tau` and shifted-positive polynomials `H(g),J(g)`.
Our independent recursive reduction of `tau^j` recovers its remainder and
the equivalent quotient trace/norm certificate:

```text
Tr_K(R)=-g(g-2)(g-1)(2g-3)(2g-1)J(g)/15120<0,
Norm_K(R)=g²(g-3)(g-2)²(g-1)²(2g-3)²(2g-1)²
             *(5g-7)H(g)/914457600>0.
```

The two root values are real, have positive product and negative sum, so
both are negative. The shifted coefficients in that proof establish these
signs for every real `g>=5`; the first-return interpretation requires odd
integral `g`. The independent analytic and final-text audit of that
companion proof passed. This strengthens a sign predicate on an already
closed smaller-endpoint-four family, not the existing width-eight theorem.

## 5. Reproduction and stopping object

```sh
python3 -B 04-computation/trinomial_trace_sign_empty_core_certificate_sep06.py
python3 -B -O 04-computation/trinomial_trace_sign_empty_core_certificate_sep06.py
```

All 37 explicit gates pass; normal and optimized outputs are byte-identical.
The semantic digest is
`2efc40c00ef00ae73583ebd3dbe813e5fbfc4536a1a5e937809c32d4f312e169`.
Raw SHA-256:

```text
source f9d6978cd378a7f9f2023baa07c6109a5fa36aff9ac69f52fe43cac5fe07efc4
output 484617dcc91fb4dc177db6f3a61952c6063d181dfe0f6a9389738d5462ffed09
```

Root independently audited the quotient/Newton construction, carry
admissibility and Vandermonde congruence: **PASS**. The returns referee
independently reviewed the criterion and replayed all three listed minor
arrays and all 37 gates: **PASS**.

The first symbolic test generated by this criterion is now **PROVED**:
[the carry family `(-15,2g-15,3g-15)`](trinomial_width15_empty_core_returns_sep06.md),
`g>=8`, `gcd(g,15)=1`, has negative quotient trace and positive norm for
every parameter. Its quadratic first row and degree-five second row retain
one lower carry; both independent proof audits passed. The next precise
test is the cubic family `(-21,2g-21,3g-21)`, integral `g>=11`,
`gcd(g,21)=1`, whose degree-seven doubled row has the same lower carry.
All three signed leading minors are required there. That cubic family and
general negative definiteness remain **OPEN** at this checkpoint.
