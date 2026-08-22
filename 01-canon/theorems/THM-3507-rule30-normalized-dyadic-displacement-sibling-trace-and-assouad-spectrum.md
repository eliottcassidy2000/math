---
id: THM-3507
title: "Rule 30 normalized dyadic displacement, sibling trace, and Assouad spectrum"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Normalizing each dyadic
  phase displacement at its first
  live state bit gives an odd 2-adic unit.  The ordered sum of its two sibling
  sheets is exactly the next normalized displacement, multiplied by the next
  innovation-gap power of two.  Exact local covering ratios identify the
  Assouad and lower-Assouad dimensions of the seed-orbit closure with the
  upper and lower Banach densities of the innovation word.  No Rule 30 prize
  consequence is claimed.
source: root/rule30-normalized-displacement/2026-08-16
audit: >
  PASS (2026-08-16), independent proof, scope, and adversarial replay audit.
  The auditor rederived the all-phase sibling trace, section-gap indexing,
  borrow transducer, hard seam scalar, exact local covers, Banach-density
  formulas, no-111 bound, and uniform-perfectness equivalence.  An independent
  width-25 coordinate path reproduced all finite counts.  Ordinary, optimized,
  and stored transcripts are byte-identical.
depends_on:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
  - THM-3500-rule30-dyadic-section-cut-defect-and-cross-depth-valuation-carrier
  - THM-3501-rule30-universal-cover-green-potential-and-slack-holonomy-seam
  - THM-3503-rule30-odometer-ultrametric-regrading-and-orbit-closure-dimensions
related:
  - THM-2538-anchored-transverse-gain-and-common-ancestry-arrival-boundary
  - THM-3488-rule30-inward-slack-monicity-and-parity-cartier-ramification
  - THM-3502-rule30-four-fifths-staircase-entropy-and-sixteen-fifths-compiler
script: 04-computation/rule30_normalized_displacement_assouad_thm3507.py
output: 05-knowledge/results/rule30_normalized_displacement_assouad_thm3507.out
script_sha256: 08823b3be627348812f51d55f7c0325393bc18ba0f6525d5b58415d65f948bda
output_sha256: 9cb0b757e9f4b252d93cc6cb141653e7990443244cc461bbda1f14f5687711ba
hash_basis: raw bytes
---

# THM-3507 -- Rule 30 normalized dyadic displacement, sibling trace, and Assouad spectrum

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-3500 identifies the dyadic Mealy cut defect with a physical phase-current
derivative.  THM-3501 retains an event-graded lift of its leading bit.
THM-3503 identifies innovation depths with the radial distortion of the
odometer embedding.  The missing cross-depth object is the complete
**oriented 2-adic displacement**, normalized at its first live coordinate.

## 1. Inheritance and conventions

Write

```text
Phi(x)=x xor ((2x) or (4x)),    R_t=Phi^t(1).          (1)
```

Let `P_w` be the seed period modulo `2^w`, and put

```text
epsilon_w=log_2(P_(w+1)/P_w),
E_w=log_2 P_w,
kappa_1<kappa_2<...={w:epsilon_w=1},
v_m=kappa_(m+1),
q_m=2^m,             d_m=v_(m+1)-v_m.                (2)
```

Thus `P_(v_m)=q_m` and `P_(v_m+1)=2q_m`.  THM-3458 and
THM-3503 give a homeomorphism

```text
iota: Z_2 -> X=closure{R_t:t>=0},
iota(t+1)=Phi(iota(t)),
nu_2(iota(t+q_m)-iota(t))=v_m                         (3)
```

for every `t in Z_2`.  Differences below are ordinary signed 2-adic
differences, not XOR.

The inheritance pass is: THM-3500's lossless dyadic cut defect is the closest
mechanism; equal XOR holonomy with opposite signed changes is the hostile;
the repair is to normalize the full ordered displacement; and the least-used
sidecar is the ordinary subtraction borrow state.

## 2. Exact normalized sibling trace

Define

```text
U_m(t)=(iota(t+q_m)-iota(t))/2^(v_m).                 (4)
```

Equation (3) makes `U_m(t)` an odd 2-adic unit.  Telescoping the two ordered
child intervals gives

```text
2^(v_m)(U_m(t)+U_m(t+q_m))
 =iota(t+2q_m)-iota(t)
 =2^(v_(m+1))U_(m+1)(t).
```

Therefore, for **every** phase `t`,

```text
boxed:
U_m(t)+U_m(t+q_m)=2^(d_m)U_(m+1)(t),                 (5)

d_m=nu_2(U_m(t)+U_m(t+q_m)).                         (6)
```

In particular `U_m(t+q_m)=-U_m(t) mod 2^(d_m)`, and the congruence fails
modulo `2^(d_m+1)`.  The innovation gap is a uniform sibling-cancellation
depth, not a statistic at the seed alone.  At the marked origin,

```text
U_m(0)=(R_(2^m)-1)/2^(v_m).                          (7)
```

The first delayed lift supplies a physical hostile:

```text
R_4=401,   R_8=102849,   v_2=4,   v_3=6,
U_2(0)=25, U_2(4)=6403,
25+6403=4*1607=2^(d_2)U_3(0).                        (8)
```

## 3. Boolean current, signed unit, and borrow

Put

```text
H_j^(q)(t)=bit_j(iota(t+q))+bit_j(iota(t)) in F_2.    (9)
```

For `q=q_m`, `j>0`, and `0<=t<q`, THM-3500 proves

```text
H_j^(q)(t+1)+H_j^(q)(t)=(rho_q Delta_(j-1))(t)       (10)
```

and `v_(m+1)=1+min{k>=0:M_0(Delta_k)=1}`.  Since `Delta_k=0`
through `k=v_m-1`,

```text
boxed:
d_m=min{r>=1:M_0(Delta_(v_m+r-1))=1}.                (11)
```

For `m>=1`, the first-carry criterion is

```text
d_m=1
 iff M_1(p_(v_m-2)^(q_m))=1
 iff U_m(t)+U_m(t+q_m) is nonzero mod 4.              (12)
```

The passage from Boolean holonomy to `U_m` needs more than `H`.  For a fixed
phase, from state coordinate `v_m` upward, put

```text
x_j=bit_(v_m+j)(iota(t)),
h_j=H_(v_m+j)^(q_m)(t),
y_j=x_j xor h_j.                                     (13)
```

The states agree below `v_m`, so ordinary subtraction starts with borrow
`c_0=0`.  Its exact two-state transducer is

```text
u_j=h_j xor c_j,
c_(j+1)=((1-y_j)x_j) or (c_j(1-h_j)).                (14)
```

The bits `u_j` are exactly the digits of `U_m(t)`.  Formula (14) follows by
testing whether `y_j-x_j-c_j` is negative.  Thus the Boolean current needs
the ordered base state and borrow to recover the signed unit.

The cheapest hostile uses zero-extended one-bit words: `(x,y)=(0,1)` gives
`+1`, whereas `(1,0)` gives `-1 mod 4`, although both have `h_0=1`.  Physically,
THM-3500's `Delta_4=0011` has even parity before `Delta_5=0111` determines
`d_2=2`.

At a hard innovation depth `k=v_m<q_m=P_k`, THM-3501's seam obeys

```text
Omega_a(1)=H_k^(q_m)(a)=U_m(a) mod 2=1,   0<=a<k.    (15)
```

It is an event-graded lift of the leading unit digit, not of the higher signed
digits.  This is THM-2538's common-ancestry lesson: form the ordered pair
before taking its XOR or marginal shadow.

## 4. Exact local coverings

Let `B_w(x)` be the ambient closed state ball of points agreeing with `x`
modulo `2^w`, and let `N_(w+L)(A)` be the least number of precision-`w+L`
state balls covering `A`.  Under `iota`, every nonempty `B_w(x) intersect X`
is one phase coset modulo `P_w`.  Refining it to precision `w+L` splits it
into the phase cosets modulo `P_(w+L)`.  Hence, uniformly in `x`,

```text
boxed:
N_(w+L)(B_w(x) intersect X)
 =P_(w+L)/P_w
 =2^(E_(w+L)-E_w).                                  (16)
```

## 5. Assouad spectra and uniform perfectness

Equation (16) gives

```text
dim_A X
 =lim_(L->infinity) (1/L) sup_(w>=1)(E_(w+L)-E_w),   (17)

dim_L X
 =lim_(L->infinity) (1/L) inf_(w>=1)(E_(w+L)-E_w).   (18)
```

Arbitrary radii lie between consecutive powers of two; rounding changes the
logarithmic scale and covering depth by at most one.  If `A_L` and `a_L` are
the displayed supremum and infimum numerators, then
`A_(L+M)<=A_L+A_M` and `a_(L+M)>=a_L+a_M`.  Fekete's lemma supplies the
limits.  They are the upper and lower Banach densities of `epsilon`.

THM-3463 proves that `epsilon` contains no `111`.  Every length-`L` window
contains at most `ceil(2L/3)` innovations, so

```text
boxed: dim_A X<=2/3.                                 (19)
```

The abstract word `(110)^infinity` makes this sharp from no-`111` alone.

The distance spectrum from every point is `{2^(-v_m):m>=0}`.  Consequently

```text
boxed:
sup_m d_m<infinity
 iff dim_L X>0
 iff X is uniformly perfect.                         (20)
```

If `d_m<=G`, every state-width interval of length `L` contains at least
`L/G-O(1)` innovations, so `dim_L X>=1/G`; the same fact fills every
sufficiently small annulus (the constant `2^(-(G+1))` is safe for closed
balls).  Conversely, an unbounded gap gives arbitrarily long
innovation-free state-width intervals and empty metric annuli, forcing
`dim_L X=0` and failure of uniform perfectness.

This is stronger than `dim_H X>0`: THM-3503 shows that positive Hausdorff
dimension controls average gaps `v_m=O(m)`, while rare individual gaps may be
unbounded.  The quantifier in (18) is load-bearing.  In innovation-index
language the **state-scale span** `v_n-v_m` must tend to infinity; merely
taking `n-m->infinity` misses one-step gaps of unbounded length.

## 6. Next certificate target

The sharp next object is a finite, physically closed, owner/borrow-aware
renormalization graph whose actual path carries edge cost `d_m`.  It must
retain ordered child orientation, borrow, calibrated phase origin, sufficient
slack/Cartier owner data, and an overflow state until boundedness is proved.

If a finite graph and potential `psi` satisfy

```text
d(e)+psi(head(e))-psi(tail(e))<=C,                    (21)
```

then telescoping gives `v_m<=Cm+O(1)` and `dim_H X>=1/C`.  A finite edge
maximum `G` additionally gives `dim_L X>=1/G` and uniform perfectness.  This
is a max-plus dual of THM-3502's Perron subeigenvector: Perron weights control
path counts; (21) controls accumulated lift cost.  THM-3502's spurious
overlap path is the mandatory hostile: a certificate on an overgraph is sound,
whereas a bad spurious cycle only refutes that quotient.

## 7. Preservation, loss, and scope

| map | preserved | destroyed | sidecar |
|---|---|---|---|
| displacement -> `U_m` | all signed digits after first live scale | common `2^(v_m)` | `v_m` |
| `U_m` -> leading XOR bit | first holonomy bit | orientation and higher digits | base state + borrow |
| phase cylinder -> state ball | exact local covers | marked phase | physical origin |
| innovation word -> dimensions | upper/lower Banach densities | marked center bits | owner/borrow dynamics |

No tournament is intrinsic: the siblings are ordered interval boundaries.

Run

```bash
python3 04-computation/rule30_normalized_displacement_assouad_thm3507.py
python3 -O 04-computation/rule30_normalized_displacement_assouad_thm3507.py
```

and compare both byte-for-byte with the stored output.  The companion has no
`assert` gates.  It checks the width-24 quotient, every available sibling
phase, independent borrow decoding, section-defect valuations, all finite
local covers, (8), and dense/sparse/burst-gap controls.

The universal statements are proved above; the finite universe is an
implementation control.  This theorem proves no bound on actual `d_m`, no
positive lower or Hausdorff dimension, no center nonperiodicity, balance, or
prediction lower bound, and no Rule 30 prize.  No literature novelty is
claimed.
