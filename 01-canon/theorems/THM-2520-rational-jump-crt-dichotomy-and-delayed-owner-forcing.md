---
id: THM-2520
title: "Rational-jump CRT dichotomy and delayed-owner forcing"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Let a periodic rational
  step response have endpoint conductor D=13^nu d with gcd(d,13)=1, and
  aggregate its true jumps modulo d.  Once the Perron depth clears nu, this
  finite prime-to-13 current is exactly the jump current of the Perron
  response, up to scaling and a CRT permutation.  It vanishes exactly when
  the Perron response is constant and exactly when all twelve last-digit
  Fourier ladders vanish; otherwise every ladder is nonzero.  An exact
  cosecant/Parseval identity gives a common energy floor.  Moving any
  positive BV future owner R further steps beyond the collision then makes
  all twelve owner-weighted collision colours strictly positive at an
  explicit O(13^(-R)) threshold.  The zero branch can contain a nonconstant
  Boolean balanced inverse-fibre tile; for an integer-valued response it is
  a constant-multiplicity inverse-fibre multisection, so a prime-to-thirteen
  denominator in the mean excludes it immediately.  An undelayed owner can
  still miss a nonzero ladder.  This proves a finite jump-current reduction,
  not that every live LRC response has nonzero current, an oriented loop, a
  C_7 owner-clock section, a row exclusion, or LRC(14).
source: codex-2026-07-27-rational-jump-crt-dichotomy
depends_on:
  - THM-2518-perron-inverse-branch-owner-word-cospan-recovery
  - THM-2519-last-digit-collision-drift-and-k13-dirichlet-boundary
related:
  - THM-2441-septimal-ancestry-event-period-collapse
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2517-cubic-anchored-spectrum-flt3-and-three-time-boolean-lift
script: 04-computation/lrc14_rational_jump_crt_thm2520.py
output: 05-knowledge/results/lrc14_rational_jump_crt_thm2520.out
script_sha256: c6a6bbd10d4e0e9003eeb424101894b633655a3926773bd7a266c1b5dc12f69c
output_sha256: 6071f203e593cefde442bdd7722c9cb8d24cc7e1ba0572868f9474b82430c4aa
independent_referee: 04-computation/lrc14_rational_jump_crt_owner_forcing_thm2520.py
independent_output: 05-knowledge/results/lrc14_rational_jump_crt_owner_forcing_thm2520.out
independent_referee_sha256: 3f1deeb09402a559f7ea29de59d024a134691f27e0ee38b08ec1be6890f83ab9
independent_output_sha256: 347722718137ec5e54c065807d59a13d65b3ce9fe177e5c2bfc47645fe20db06
typed_probe: 04-computation/lrc14_typed_owner_endpoint_current_thm2520_probe.py
typed_probe_output: 05-knowledge/results/lrc14_typed_owner_endpoint_current_thm2520_probe.out
typed_probe_sha256: 02d247bbe42f26b753916afeda3d83dd4e48b32f74d17b41cb2214abc54c3ec1
typed_probe_output_sha256: e14f767c1e89b6b91c5491e367a0444781e01bc5e91ae252e8df39bea0df05a4
hash_basis: working-tree bytes (LF)
---

# THM-2520 -- endpoint jumps decide delayed collision colour

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2518 replaces free moment arms by lawful inverse-branch cospans.
THM-2519 shows that their last collision digit is a `K_13` conditional
variance, but a bulk mixed moment can lie in its zero branch.  For rational
step responses, that branch is not diffuse or asymptotic.  It is a finite
signed current on the part of the endpoint grid coprime to thirteen.

The exact reduction is

```text
rational endpoint jumps
  -> aggregate modulo the 13-free conductor d
  -> zero: a balanced inverse-fibre tile and constant Perron response;
  -> nonzero: all twelve last-digit ladders survive
              and every sufficiently delayed owner sees them.             (1)
```

The theorem does not assert which branch contains every live LRC packet.
It replaces that analytic question by one finite CRT current.

## 1. The prime-to-thirteen jump current

Let `F:T->Q_(>=0)` be a periodic step function.  Use its true net jumps

```text
Delta_j=F(x_j+)-F(x_j-),                 x_j=b_j/D,             (2)
```

after combining repeated endpoints.  Choose a common endpoint conductor

```text
D=13^nu d,                         gcd(d,13)=1.                 (3)
```

Minimality of `D` is not required.  What is required is that (2) use the
true net jump at each grid point and that (3) expose the full `13`-primary
part of the chosen grid.

The **coprime jump current** is the rational `C_d`-vector

```text
C_t=sum_(j: b_j=t mod d) Delta_j.                              (4)
```

Periodicity gives `sum_t C_t=0`.  Let

```text
L-1>=nu,                 M=13^(L-1),
U=13^(L-1-nu) mod d,     h=P_M F,                              (5)
```

where

```text
(P_MF)(y)=1/M sum_(r=0)^(M-1)F((y+r)/M).                      (6)
```

As an equality of periodic jump measures,

```text
dh=1/M sum_(v in C_d) C_(U^(-1)v) delta_(v/d).                (7)
```

Indeed, one branch of (6) crosses `x_j` precisely when
`y=Mx_j mod 1`.  Equation (5) clears the complete `13`-primary endpoint
denominator, so only `b_j mod d` remains; multiplication by `U` merely
permutes that grid.

Use the Fourier convention

```text
Fhat(n)=integral_T F(x)exp(-2 pi i n x)dx,
Chat(b)=sum_t C_t exp(-2 pi i bt/d).                           (8)
```

Distributional integration by parts and `hhat(k)=Fhat(kM)` give, for
`k!=0`,

```text
hhat(k)=Fhat(kM)=Chat(Uk)/(2 pi i kM).                         (9)
```

Thus the current is not merely a necessary endpoint statistic: it is the
complete high-frequency Perron numerator.

## 2. CRT all-or-all ladder dichotomy

Let `zeta=exp(2 pi i/13)` and, for `a in F_13`, put

```text
H_a(y)=1/13 sum_(u in F_13)zeta^(-au)h(y+u/13).              (10)
```

Then

```text
H_a(y)=sum_(k=a mod 13)hhat(k)exp(2 pi i ky).                 (11)
```

For any fixed `a!=0`, the congruence class `k=a mod 13` runs through every
class modulo `d`, because `(13,d)=1`.  Equation (9) and finite Fourier
inversion on `C_d` therefore prove the equivalent conditions

```text
C=0;

h=integral_T F;

H_a=0 for one prescribed a!=0;

H_a=0 for every a!=0.                                        (12)
```

Consequently,

```text
C!=0  iff  H_a!=0 for all twelve a in F_13^*.                (13)
```

This is a literal CRT mechanism.  It does not use rational Galois
propagation: a single nonzero `d`-frequency of `C` can be paired with each
prescribed nonzero residue modulo thirteen.

There is an equally useful position-space description.  If `h` jumps at
`t/d`, then the thirteen translates in (10) jump at

```text
t/d-u/13,                         u in F_13.                   (14)
```

The map `(t,u)->13t-du mod 13d` is bijective.  Hence no two translated
jumps can cancel at one point.  This explains geometrically why every
character sees a nonzero jump whenever `C!=0`.

The branch type is stable for every `L-1>=nu`; increasing `L` only scales
the current by `1/M` and permutes it by a power of thirteen modulo `d`.
Its labelled form is periodic in `L` with period dividing `ord_d(13)`.

## 3. Exact ladder energy and a quantitative floor

Put

```text
Q=13d.
```

For `a!=0` and `b in C_d`, let `r_(a,b)` be the unique representative in
`{1,...,Q-1}` satisfying

```text
r_(a,b)=a mod 13,                    U r_(a,b)=b mod d.        (15)
```

Grouping (11) by residue classes modulo `Q` and using

```text
sum_(n in Z) 1/(r+nQ)^2
 =pi^2/Q^2 csc^2(pi r/Q)                                      (16)
```

gives the exact identity

```text
||H_a||_2^2
 =1/(4M^2Q^2) sum_(b in C_d)
    |Chat(b)|^2 csc^2(pi r_(a,b)/Q).                          (17)
```

The `b=0` term is harmless: `Chat(0)=sum_t C_t=0`, while `a!=0` keeps
`r_(a,0)` away from zero.  Unnormalized finite Parseval and
`csc^2>=1` yield the common floor

```text
||H_a||_2^2
 >= [sum_t C_t^2]/(4*13^2*M^2*d)                for all a!=0. (18)
```

This exact norm is stronger than a nonvanishing assertion.  It invoices the
`M^(-2)` cost of asking a late inverse branch to remember its last digit.

Let

```text
B=||F||_infinity,                       V=Var(F).              (19)
```

The jump formula also gives

```text
||H_a||_infinity<=B,
Var(H_a)<=Var(h)<=V/M,
Var(|H_a|^2)<=2BV/M.                                         (20)
```

## 4. A sufficiently delayed common owner sees every colour

Let `G:T->[0,1]` be a positive BV future owner--word factor, and write

```text
rho=integral_T G>0,                       W=Var(G).            (21)
```

For `R>=0`, set

```text
K=13^(R+1).                                                    (22)
```

The fixed-last-digit collision table with the owner placed `R` steps after
the collision is

```text
B_u^[L,R]
 =1/M sum_(e=0)^(M-1) integral_T
    G(13^(L+R)x)F(x)F(x+(u+13e)/13^L)dx.                     (23)
```

The exact THM-2519 reduction, with the later common factor retained, is

```text
B_u^[L,R]
 =integral_T G(Ky)h(y)h(y+u/13)dy.                           (24)
```

Because `K/13` is an integer, `G(Ky)` is invariant under every
`y->y+u/13`.  Therefore its normalized `u`-Fourier coefficients are exact
weighted norms:

```text
Bhat^[L,R](a)
 :=1/13 sum_u B_u^[L,R] zeta^(-au)
 =integral_T G(Ky)|H_a(y)|^2dy
 >=0.                                                        (25)
```

There is also an exact fixed-delay all-or-all law.  Away from its finite
jump set, the thirteen numbers `h(y+u/13)` are rational.  If one nontrivial
`H_a(y)` vanishes, then the degree-twelve rational polynomial

```text
sum_(u=0)^12 h(y+u/13)X^u
```

vanishes at a primitive thirteenth root.  Irreducibility of `Phi_13` forces
all thirteen coefficients to be equal, so every nontrivial `H_a(y)`
vanishes.  Since the weight in (25) is nonnegative,

```text
at every fixed R, either all twelve Bhat^[L,R](a) vanish
or all twelve are strictly positive.                          (25a)
```

This pointwise form needs rational response values, but no rationality of
the owner weight.

The periodic BV covariance estimate and (20) give

```text
|Bhat^[L,R](a)-rho||H_a||_2^2|
 <=BVW/(6MK).                                                 (26)
```

Put

```text
S_C=sum_t C_t^2.
```

If `C!=0`, then all twelve coefficients in (25) are strictly positive as
soon as

```text
K > 2*13^2*d*M*B*V*W/(3*rho*S_C).                            (27)
```

When `W=0`, the error is zero and no delay threshold is needed.  Since
arbitrarily large `R` are available in either parity class, (27) is
compatible with either septimal orientation convention.

If `F` and `G` are rational step functions with rational endpoints, every
entry of (23) is rational.  In the lawful response setting, expand `F` into
its packet factors before the old-root sums.  Then (25) is simultaneously a
Boolean-before-Fourier construction, an analytic positive norm, and an exact
rational collision table.

In the zero branch, (12) makes (24) completely rigid:

```text
C=0  ->  B_u^[L,R]=rho (integral_T F)^2
          for every u and every R.                            (28)
```

Equations (27)--(28) are the promised delayed-owner dichotomy.

Equivalently, with `F` fixed,

```text
C=0
 iff every positive BV owner at every later delay has zero collision drift;

C!=0
 implies every positive BV owner has all twelve colours
         at every sufficiently large delay.                   (28a)
```

For the reverse direction in the first line, take the constant owner
`G=1`; a zero collision drift is the sum of the twelve unweighted ladder
norms, so (12) forces `C=0`.

### The constant branch is an inverse-fibre multisection

There is a cheaper necessary test when the response is lattice-valued.
Suppose that, for some `lambda>0`,

```text
F(x) in lambda Z                         almost everywhere.  (28b)
```

If `P_MF=A` is constant, then almost every inverse fibre has the same
lattice weight:

```text
M A/lambda
 =sum_(r=0)^(M-1) F((y+r)/M)/lambda
 in Z.                                                        (28c)
```

For a nonnegative integer-valued response this is literally a
constant-multiplicity inverse-fibre multisection.  In particular,

```text
A in (lambda/M)Z.                                             (28d)
```

If `F` is integer-valued and `A=p/q` is reduced, the zero-current branch
forces `q|M`.  Since `M` is a power of thirteen, a prime-to-thirteen factor
of `q`, or a thirteen-primary denominator deeper than `M`, proves `C!=0`
without enumerating a single endpoint.  This condition is only necessary:
a mean in `M^(-1)Z` does not force the fibre multiplicities to be constant.

## 5. Sharp boundaries and hostile controls

Each hypothesis has a concrete job.

### The depth must clear the `13`-primary conductor

Let

```text
F(x)=1_[0,1/2)({13x}).
```

Its endpoint conductor is `26=13*2`, and its aggregate coprime current is
nonzero.  At `L=1`, however, `h=F` is `1/13`-periodic, so every nontrivial
`H_a` vanishes.  The theorem begins at `L=2`, exactly as (5) requires.

### The constant branch can be nonconstant and Boolean upstairs

For

```text
F=1_[0,1/13),                       L=2,
```

the response is nonconstant but

```text
P_13F=1/13,                         C=0.                      (29)
```

Hence every collision entry is `rho/169`.  More generally, a Boolean
response can choose the same number of active inverse branches over every
coprime base cell.  The zero branch is an exact balanced fibre tile, not
pointwise constancy of the original response.

### A nonzero current need not meet the undelayed owner

Take

```text
F=1_[0,1/14),                       L=1,
G=1_[13/14,1).
```

The coprime current is nonzero, but over the support of the immediate future
owner all thirteen predecessors miss `F`; every weighted colour in (25) is
zero at `R=0`.  The extra mixing delay in (27) is therefore load-bearing.

### Use true jumps, not a presentation ledger

Repeated endpoint descriptions must first be combined into the net jump in
(2).  Artificial subdivision creates zero jumps and cannot create a charged
current.  Changing a valid common conductor only refines or permutes the
same conclusion once its full `13`-part is cleared.

### Rational response values are sharp for fixed-delay all-or-all

Choose a small rational base interval `I`.  On its thirteen disjoint inverse
images, prescribe the real step values

```text
2+2 cos(2 pi a_0r/13),                       r in F_13,
```

and take the future weight supported on `I`.  This fibre has Fourier support
only at `+-a_0`, so the owner sees two positive colours and ten zeros.  Thus
the rational `Phi_13` argument behind (25a) cannot be extended to arbitrary
real step values.  The eventual delayed conclusion from the jump CRT and
mixing still holds without this pointwise shortcut.

## 6. Lawful cospan scope and the remaining LRC object

In (23), the two histories first coalesce at depth `L`; the factor
`G(13^(L+R)x)` is one genuine common owner `R` steps later.  If

```text
G=E_j (Q o T^K_0),
```

the entire block must be shifted together, and its terminal word occurs at
depth `L+R+K_0`.  Expanding each response density before its root sum makes
every summand a literal Boolean product on the inverse-branch fibre product.
The distinguished old deep leg and its diagonal zero survive exactly as in
THM-2518.

What THM-2520 closes is precise:

```text
nonzero rational-step last-digit ladder
  -> some positive common future owner sees all twelve colours.             (30)
```

What remains is equally precise.  For each live response packet, form (4).
Either its current is nonzero and (27) applies, or its late Perron response
is the balanced constant branch (28).  The latter must be excluded or used
combinatorially; bulk cubic or square charge alone does not exclude it.

The live test is sparse: choose one common endpoint denominator, list only
actual boundary jumps, accumulate them by residue modulo `d`, and test the
resulting rational vector for zero.  There is no need to expand the full
`13^nu` inverse tower.  Different response cells can still have different
currents, and signed mixed-table combinations can cancel; a marked-cell
current is not automatically intrinsic ANOVA drift.

For the lawful response densities used here,

```text
F_j=sum_r h_(j,r),
```

each summand is Boolean, so `F_j` is integer-valued.  Reduce its mean first:
the denominator test (28c)--(28d) excludes the balanced branch whenever the
reduced denominator has a prime-to-thirteen factor.  Only the remaining
pure-thirteen cases require the sparse current computation.

### Canonical typed-row positive control

The exact typed probe applies this two-stage test to the canonical
THM-2309/2334/2349 numerical row

```text
(H,q_1,...,q_5,c_1,c_2,c_3)
 =(1,14,27,40,53,66,13,2197,742586),                         (31)
```

with owner `c_1=13`, word clock `169`, and the three lawful word strata
`Q_a,Q_b,Q_ab`.  The common grid has conductor
`D=13^8*61,704,720`, and all `26` deep endpoints align with the `c_3` word
boundary.  Before and after summing the old deep factor, respectively, the
six current `(support,L1)` pairs are

```text
Q_a:  (14,376112),  (14,694328),
Q_b:  (22,3867136), (22,7412064),
Q_ab: (4,347100),   (4,665252).                              (32)
```

The mean-denominator screen already proves nonzero current in five cases.
The deep-summed `Q_ab` mean has pure denominator `13^7`, so the screen is
sharp there, but its explicit four-residue current is still nonzero.  Thus
all three typed word strata have nonzero current both before and after the
old deep sum, and THM-2520 supplies every collision colour after a common
sufficiently late owner delay in this positive-control instance.

The constructor explicitly does **not** identify (31) with one of the `165`
scalar-cover rows.  Equation (32) is therefore a canonically typed lawful
signal and a test of deep compatibility, not a row exclusion.  The scalar
cover identification, mixed-table cancellation, and Boolean emission remain
open.

This theorem does not orient `u` against `-u`, turn the positive norm into a
signed owner-loop current, or choose the free `C_7` owner clock of THM-2517.
The next joint carrier should retain the seven disjoint Latin cells `P_c`
before replacing them by the scalar union `sum_c P_c`.  Their clock character

```text
Z=sum_(c in C_7) zeta_7^c P_c                                (31)
```

has its own finite jump current.  A nonzero current would retain a
`C_7`-typed ladder before the positive norm is taken; a zero current would
mean exact clock balance on every late inverse fibre.  The norm in (25)
forgets that phase, so a lawful polarization/emission sidecar is still
required.  This is a sharper next object, not a claimed `C_91` closure.

## 7. Exact companion

Run

```bash
python3 04-computation/lrc14_rational_jump_crt_thm2520.py
python3 -O 04-computation/lrc14_rational_jump_crt_thm2520.py
python3 04-computation/lrc14_rational_jump_crt_owner_forcing_thm2520.py
python3 -O 04-computation/lrc14_rational_jump_crt_owner_forcing_thm2520.py
python3 04-computation/lrc14_typed_owner_endpoint_current_thm2520_probe.py
python3 -O 04-computation/lrc14_typed_owner_endpoint_current_thm2520_probe.py
```

Both executions match the stored transcript.  The exact `Fraction`
companion checks:

- `640` Perron jump-current cells and `128` mean identities at two legal
  depths on deterministic rational `13*5` profiles;
- a nonconstant Boolean balanced tile whose two late Perron profiles are
  exactly constant `1/13`;
- all `65` shifted-grid CRT points and all `120` source-by-coprime ladder
  addresses at two levels; and
- the strict BV/Parseval invoice at owner scale `13^3`, all twelve positive
  collision toothpick energies, and all twelve exact cyclotomic
  nonvanishing tests.

The independent referee uses a separate finite-grid implementation.  It
checks `144` conductor/depth cases, `540` direct Perron cells, `1,728` CRT
coverage rows, `144` variation contractions, `74/70` constant/charged
branches, the pure-`13` hostile, the prime-to-`13` control drift `42/169`,
three lattice-fibre invoices, and the independent safe delay constant `20/3`
obtained from `pi^2<10`.  The typed sparse-endpoint probe checks the six
currents in (32), all `26` deep endpoint alignments, and both denominator
controls.  All normal/optimized runs match their stored transcripts.

Independent audits rederived (7), (9), every CRT quantifier, the exact
cosecant identity (17), the Parseval normalization, the variation invoice,
the threshold (27), all three hostile boundaries, and the physical timing
`G(13^(L+R)x)`.  No live row is removed.  LRC(14) remains open. **QED.**
