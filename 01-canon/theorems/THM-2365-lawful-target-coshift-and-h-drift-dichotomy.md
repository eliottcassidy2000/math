---
id: THM-2365
title: "Lawful target co-shift and H-drift dichotomy"
status: >
  PROVED CANDIDATE UNDER INDEPENDENT AUDIT. The full lawful
  endpoint/deep thirteen-shift tensor H(r,s,t) is nonnegative and
  vanishes on the diagonal plane r=t. Its finite transform B(a,b,h)
  has target vector (b,a+h), every target-vector line sums to zero,
  and its nonzero-target drift energy is exactly the squared distance
  from the circulant subspace H(r,s,t)=G(r-t). If that drift is
  positive, at least 1/13 of it has a nonzero deep colour; absolute
  endpoint and deep-mode collapse then gives an exact 91-unit m and
  ordinary frequency X with a nonzero fixed-triangle target fibre.
  For the rational LRC interval sets, all twelve nonzero untwisted
  deep colours survive. Delayed-word drift converges quantitatively
  to positive word mass times bare-owner drift at rate O(1/R).
  Thus positive bare drift lands a target at every sufficiently large
  delayed clock; the exact residual law is H_E(r,s,t)=G_E(r-t).
  The extracted triangle need not be the prior marked nonzero-total
  triangle. No scalar row is excluded and LRC(14) remains open.
source: codex-2026-07-25-lawful-target-coshift
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2354-deep-shift-comb-cover-and-grouped-unit-current
related:
  - THM-2343-deep-comb-affine-target-catalyst
  - THM-2353-target-plaquette-holonomy-and-cross-axis-slice
  - THM-2364-anchored-corner-forces-mixed-deep-blocker-colour
script: 04-computation/lrc14_lawful_target_coshift_h_drift_thm2365.py
output: 05-knowledge/results/lrc14_lawful_target_coshift_h_drift_thm2365.out
script_sha256: b71508ac71cffc427c8f01722309bd49a3d9a3dca7aa66468fab31b6d83fb425
output_sha256: 9df085bbf9493278f4de45d4ee8c09c2d1732c9258b91646465f3c927624fda4
hash_basis: working-tree bytes (LF)
---

# THM-2365 -- H-drift is exactly departure from the inverse line

**PROVED CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2354 gives a positive factor-coloured deep current but forgets the
endpoint phases needed by the target quotient. Restoring those phases
does not merely introduce an uncontrolled cancellation: it produces an
exact finite tensor with one rigid null direction.

The resulting statement has two faces.

```text
non-circulant H-drift
  -> an exact 91-unit deep leg with a nonzero target fibre;

zero H-drift
  <=> H(r,s,t)=G(r-t), independent of the other target shift.
                                                               (1)
```

The second line is the precise functional form of the remaining
inverse-character hostile.

## 1. Lawful target coordinates

Put `p=13` and let `c=c_3` be the deepest blocker coordinate. Use the
THM-2309 owner packet modulo `p`. Its relation kernel splits as

```text
K_p=L_p direct_sum <e_a,e_c>,

K_p/L_p isomorphic to F_p^2,                         (2)
```

where `a` is the other target blocker. Choose quotient-dual lifts

```text
eta,ell in L_p^perp/<w>
```

normalized by

```text
eta_c=0,                    ell_c=1.                 (3)
```

The two target coordinates of a relation `r` are

```text
q(r)=(eta.r,ell.r).                                 (4)
```

This is the abstract version of the concrete lawful basis

```text
eta=e_(c_2)-e_(q_1),
ell=e_(c_3)-e_(q_2)
```

in THM-2353. Raw ambient blocker coordinates are not substituted for
these quotient-dual differences.

Use THM-2334's present owner factors `I_i`, transported word factors
`J_i`, and clock

```text
R=13^k.
```

On the circle, with variable `x`, put

```text
E_(s,t)(x)
 =product_i I_i(
    w_i x-(s eta_i+t ell_i)/p
  ),

W_(s,t)(x)
 =product_i J_i(
    R w_i x-R(s eta_i+t ell_i)/p
  )
 =W(x),                                             (5)

F_(s,t)=E_(s,t)W.
```

The equality in (5) is exact because `p|R` and every base factor is
one-periodic. Thus the target action moves the present and bare endpoint
packets lawfully while leaving the delayed word unchanged.

Let

```text
d(y)=1_(||y||<1/14),

Delta_r(x)=d(c x-r/p),               r in F_p,       (6)

H(r,s,t)
 =integral_T F_(s,t)(x)E_(s,t)(x)Delta_r(x) dx
 =integral_T F_(s,t)(x)Delta_r(x) dx.               (7)
```

The second equality uses `F_(s,t) subset E_(s,t)`. Therefore

```text
H(r,s,t)>=0.                                        (8)
```

The `c`-coordinate factor of `E_(s,t)` is the complement

```text
1-d(c x-t/p),
```

by (3). Hence, up to null endpoints,

```text
H(t,s,t)=0                 for every s,t.           (9)
```

This diagonal-plane zero is exact at the indicator boundary. Poisson
smoothings need not have zero product at finite smoothing parameter.

## 2. The three colours and the actual target

Let `zeta=exp(2*pi*i/p)` and define the normalized transform

```text
B(a,b,h)
 =1/p^3 sum_(r,s,t)
    H(r,s,t)zeta^(a r+b s+h t).                    (10)
```

Use atomic left-present, word, right-present, and deep indices

```text
u,beta,v,m.
```

The relation address is

```text
r_full=u+R beta+m e_c-v.                           (11)
```

The minus translations in (5)--(6) give an atomic phase

```text
zeta^(
 -m r-s eta.(u-v)-t ell.(u-v)
).
```

Thus character orthogonality in (10) selects

```text
m=a,
eta.(u-v)=b,
ell.(u-v)=h                     mod p.              (12)
```

Because `p|R`, the word digit is target-neutral:

```text
eta.(R beta)=ell.(R beta)=0             mod p.
```

Equations (3), (4), and (11)--(12) show that the actual target vector of
`B(a,b,h)` is

```text
q=(b,a+h).                                          (13)
```

The factor colour `a` and endpoint colour `h` are separately
gauge-dependent; their sum in (13) is the preserved target coordinate.

## 3. Exact target-null circulation

Fix a target vector `(b,q_2)`. Summing (10) over all decompositions
`q_2=a+h` gives

```text
sum_a B(a,b,q_2-a)
 =1/p^2 sum_(s,t)
    H(t,s,t)zeta^(b s+q_2 t)
 =0.                                               (14)
```

The first equality is character orthogonality in `r-t`; the second is
(9). Therefore the full lawful co-shift does not transfer THM-2354's
factor-coloured mass after the deep colour is forgotten. It organizes
that mass into an exact **target-null circulation**.

At the untwisted endpoint,

```text
G(r)=H(r,0,0).
```

Consequently THM-2354's deep transform is the marginal

```text
C(a)
 =1/p sum_r G(r)zeta^(a r)
 =sum_(b,h)B(a,b,h).                               (15)
```

Thus a nonzero `C(a)` is distributed among target vectors, but its
target-vector sum cancels by (14).

## 4. H-drift is distance from one circulant subspace

Let the target-action average be

```text
(P H)(r,s,t)
 =1/p^2 sum_(u,v in F_p)
    H(r+v,s+u,t+v).                                (16)
```

It is the orthogonal projection, for normalized counting measure, onto
the functions invariant under

```text
(r,s,t)->(r+v,s+u,t+v).
```

Those functions depend only on `r-t`. In Fourier space, `P` retains
exactly

```text
b=0,                       a+h=0,
```

which by (13) is the zero-target line. Hence the exact drift energy

```text
D_H
 :=1/p^3 sum_(r,s,t)|H(r,s,t)-(P H)(r,s,t)|^2

 =sum_(
    a,b,h;
    (b,a+h)!=0
  ) |B(a,b,h)|^2.                                  (17)
```

In particular,

```text
D_H=0
 iff H(r,s,t)=h_0(r-t) for some h_0
 iff H(r,s,t)=G(r-t).                              (18)
```

The last equality evaluates at `(s,t)=(0,0)`.

On each nonzero target line, (14) and Cauchy--Schwarz give

```text
|B(0,b,q_2)|^2
 <=12 sum_(a!=0)|B(a,b,q_2-a)|^2.
```

After summing all nonzero targets,

```text
sum_(
  (b,q_2)!=0;
  a!=0
 ) |B(a,b,q_2-a)|^2
 >=D_H/13.                                         (19)
```

Therefore

```text
D_H>0
```

forces a coefficient which has both a nonzero target and a nonzero deep
colour.

The constant `1/13` is sharp for the line-sum relation alone: take the
twelve `a!=0` entries equal and the `a=0` entry equal to minus their sum.

The energy in (17) is also a positive physical finite-difference
observable. If

```text
(T_(u,v)H)(r,s,t)=H(r+v,s+u,t+v),
```

then the Hilbert-space group-average identity gives

```text
D_H
 =1/(2p^2) sum_(u,v)
    1/p^3 sum_(r,s,t)
      |H(r,s,t)-T_(u,v)H(r,s,t)|^2.                (19a)
```

No Fourier phase must be estimated to decide whether it is positive.

There is an even cheaper sufficient test. Define

```text
S(s,t)=sum_r H(r,s,t).
```

The thirteen-shift successor identity gives

```text
S(s,t)
 =integral_T F_(s,t)(x)(2-d(13c x)) dx.            (19b)
```

If (18) holds, `S` is constant. Therefore any variation among these
`169` nonnegative successor-weighted masses already forces `D_H>0` and
the target landing below. Constancy of `S` is only a first hostile test,
not a proof of zero drift.

## 5. Extraction of an exact fixed triangle

The conclusion after (19) is not only a formal finite colour. Two
absolute-convergence facts let us descend to exact integer frequencies.

For each fixed target twist and exact `m`, endpoint Parseval gives

```text
sum_X
 L^(s,t)(X)conjugate(E^(s,t)(X+m c))

 =F_(s,t)^hat(-m c),                               (20)
```

and the sum on the left is absolutely convergent by Cauchy--Schwarz.
Moreover

```text
sum_m
 |d_hat(m)F_(s,t)^hat(-m c)|

 <=||d||_2 ||F_(s,t)||_2<infinity.                 (21)
```

There are only `p^2` endpoint twists. We may therefore first separate
the absolutely convergent sum in (21) into residue classes `m=a mod p`,
then take the finite target transform, and finally use the absolutely
convergent `X`-sum in (20).

More explicitly, with `A_(X,m)(q)` denoting the fixed-triangle target
fibre of THM-2334, including the deep coefficient,

```text
B(a,b,h)
 =sum_(m=a mod p) sum_X
    A_(X,m)(b,a+h),                                 (21a)
```

with absolute `m`-then-`X` meaning supplied by (21) and (20).

If the coefficient forced by (19) is nonzero, some exact

```text
m=a mod p
```

has a nonzero target fibre, and for that `m` some ordinary frequency
`X` has

```text
A_(X,m)(q)!=0                    for q!=(0,0).       (22)
```

Because `a!=0`, this `m` is not divisible by thirteen. The centred
danger coefficient `d_hat(m)` vanishes at every nonzero multiple of
seven. Every live `m` in (22) therefore satisfies

```text
gcd(m,91)=1.                                        (23)
```

Equation (22) is a genuine fixed-`(X,m)` target-fibre landing in the
THM-2334 relation current. It need not be the preselected THM-2349
triangle whose untwisted total current is nonzero. That same-triangle
allocation remains a separate obligation.

## 6. Rationality upgrades one colour to all twelve

In the LRC application, `F` and every `D_(c,r)` are finite unions of
intervals with rational endpoints. Hence

```text
G(r) in Q.
```

We also have

```text
G(0)=0,                    sum_r G(r)>0.             (24)
```

For any `a!=0`, `zeta^a` is a primitive thirteenth root. If `C(a)=0`,
the rational polynomial

```text
sum_(r=0)^12 G(r)X^r
```

would be divisible by

```text
Phi_13(X)=1+X+...+X^12.
```

Both have degree at most twelve, so all thirteen values `G(r)` would be
equal, contradicting (24). Therefore

```text
C(a)!=0                   for every a!=0.           (25)
```

There is no row-uniform magnitude floor for each individual colour in
(25); the quantitative `1/156` maximum from THM-2354 remains the safe
uniform statement.

If `D_H=0`, (18), (15), and (25) give the complete residual spectrum

```text
B(a,b,h)=0 unless (b,h)=(0,-a),

B(a,0,-a)=C(a)!=0              for every a!=0.      (26)
```

Thus the zero-drift branch is exactly the full inverse-character line,
not an unspecified cancellation.

## 7. Delayed words reduce to a bare-owner drift test

The delayed clock supplies a further simplification. Let `Q` be the
fixed positive terminal word on the base circle, let

```text
F_R(x)=E(x)Q(Rx),
```

and write `H_R` for (7). Let `H_E` be the same tensor with `Q(Rx)`
deleted. For each shift triple put

```text
A_(r,s,t)(x)=E_(s,t)(x)Delta_r(x).
```

Every such function and `Q` has bounded variation. Set

```text
V_A=max_(r,s,t) Var(A_(r,s,t)),

V_Q=Var(Q).
```

Fourier covariance gives

```text
H_R(r,s,t)-mu(Q)H_E(r,s,t)
 =sum_(n!=0) Q_hat(n)A_(r,s,t)^hat(-nR).
```

Using

```text
|f_hat(n)|<=Var(f)/(2*pi*|n|)
```

twice and `sum_(n!=0)n^(-2)=pi^2/3` yields the uniform bound

```text
max_(r,s,t)
 |H_R(r,s,t)-mu(Q)H_E(r,s,t)|
 <=V_Q V_A/(12R).                                  (27)
```

Since `I-P` is an orthogonal contraction,

```text
sqrt(D_(H_R))
 >=mu(Q)sqrt(D_(H_E))-V_Q V_A/(12R).               (28)
```

Consequently,

```text
D_(H_E)>0
```

implies `D_(H_R)>0` whenever

```text
R>V_Q V_A/(12 mu(Q)sqrt(D_(H_E))).                 (29)
```

Conversely, if `D_(H_R)=0` along an unbounded sequence of clocks, then
(27)--(28) force

```text
D_(H_E)=0,

H_E(r,s,t)=G_E(r-t).                               (30)
```

The converse of the last implication is not claimed: even when the bare
owner is circulant, a delayed word correction may create positive drift.

## 8. Consequence on the 165-row frontier

THM-2349 supplies a positive shallow owner `E_j` on every one of the
`165` first-depth-one rows. THM-2305 partitions the fixed residual
`R_j` into three fixed terminal words. Choose one positive word
`Q_(j,sigma)` once and for all (the largest has mass at least
`eta/3`). BV mixing with this fixed word holds at every sufficiently
large clock, so the clock may be chosen arbitrarily large. Sections
1--7 therefore give the exact rowwise alternative:

```text
bare-owner H-drift is positive
  -> at every sufficiently large delayed clock,
     some exact 91-unit m and X have a nonzero target fibre;

bare-owner H-drift is zero
  -> H_E(r,s,t)=G_E(r-t), independent of s.         (31)
```

The second branch is a finite-torus X-ray or toothpick law: all endpoint
motion transverse to the deepest relative shift disappears, and only
the difference `r-t` remains. It is the all-frequency physical
nonnegative analogue of THM-2343/2353's fixed-triangle inverse-character
boundary; the two are not identified term-by-term.

No intrinsic binary relation is present in this tensor. A tournament
encoding of its shift cells would discard amplitudes, the target-action
gauge, and the circulant projection, so it is not a faithful operation
here.

The theorem does not yet prove that bare drift is positive on every
coefficient row, nor that a word correction breaks every circulant bare
row. It does not land the prior marked triangle, impose the all-coordinate
`91`-unit mask, exclude any scalar profile, or prove LRC(14). The scalar
ledger remains `165`.

The next cheapest decisive object is now exact:

```text
the 13^3 rational overlap masses H_E(r,s,t),
their orbit averages P H_E,
and the first non-circulant correction along R=13^k.              (32)
```

## 9. Exact companion

The dependency-free companion verifies, over exact rational/Gaussian
integer arrays:

- the `13^3` target-action projection and its orbit formula;
- diagonal-plane target-line cancellation;
- the sharp `1/13` unit-deep energy inequality;
- non-circulant positive and circulant hostile controls;
- the complete inverse-line spectrum of a circulant table;
- prime-cyclotomic rational rigidity for all twelve colours; and
- the constants in the BV covariance estimate.

Run

```bash
python3 04-computation/lrc14_lawful_target_coshift_h_drift_thm2365.py
python3 -O 04-computation/lrc14_lawful_target_coshift_h_drift_thm2365.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_lawful_target_coshift_h_drift_thm2365.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Independent audit is pending. QED.
