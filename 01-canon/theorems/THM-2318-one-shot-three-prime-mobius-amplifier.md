---
id: THM-2318
title: "One-shot three-prime Möbius amplifier"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Let two
  nonzero nonnegative real step functions have disjoint supports. For any
  prime ell distinct from 7 and 13, a fully minimal collision cell
  (s,1,1) in the commuting Perron directions (13,7,ell) isolates a common
  frequency 13^(s-1)n with gcd(n,91 ell)=1 and
  n<=(91 ell)J_AJ_B-1. Thus the two one-shot directions occur in the
  square-free Möbius selector but not in the landed frequency. For
  THM-2306's owner-normalized word current, the frequency is
  c 13^(s-1)n and n<=8(91 ell)S^2-1. The exact THM-2299
  multiplier-four carrier with ell=65537 has a fully minimal cell
  (3,1,1), so it simultaneously lands at thirteen-adic grade three with
  a whole outside multiplier coprime to 91. This removes THM-2313's
  grade-versus-seven-valuation tradeoff but does not select the root
  residue, target gain, or an incident pair edge. No scalar row is
  excluded and LRC(14) remains open.
source: codex-2026-07-25-one-shot-mobius-amplifier
depends_on:
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2306-owner-normalized-disjoint-support-and-first-collision-shell
  - THM-2313-biprime-pareto-collision-frontier-and-91-unit-current-shell
related:
  - THM-2293-quadratic-root-energy-raises-the-ancestry-shell
  - THM-2302-same-label-expiration-dichotomy-and-pure-terminal-shell-no-go
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
script: 04-computation/lrc14_one_shot_mobius_amplifier_thm2318.py
output: 05-knowledge/results/lrc14_one_shot_mobius_amplifier_thm2318.out
script_sha256: eb5523dbe8111bd779899b80e308bfe5fc36f7f3483ec2436094faef2447d155
output_sha256: de2578cdc9ffea2b1dc99010d14bef0534eff9d6101ec5c2ed23e792267c3cf9
hash_basis: working-tree bytes (LF)
---

# THM-2318 -- use primes as one-shot differences, not permanent clocks

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2313 introduced a second Perron direction in order to remove multiples
of seven from the residual Fourier core. Its exact hostile frontier exposed
a sharp tradeoff. The grade-three corner `(3,7)` carries an added factor
`7^6`, while the corner `(8,1)` carries no added factor of seven but has the
wrong thirteen-adic grade.

That tradeoff comes from asking one direction to do two jobs:

```text
seven must appear in the Möbius selector,
and repeated seven-pushes must supply all auxiliary spatial mixing.       (1)
```

These jobs can be separated. Keep seven for exactly one Boolean difference,
and use a second prime `ell` for exactly one additional difference. At a
fully minimal cell `(s,1,1)`, both one-shot primes occur in the selector,
but their exponents in the predecessor scale are zero:

```text
directions:       13, 7, ell;
minimal cell:     (s,1,1);
shell scale:      13^(s-1) 7^0 ell^0;
residual selector:
  gcd(n,13*7*ell)=1.                                 (2)
```

The top push still uses the large product `7 ell` to create collision.
The landed frequency does not. This is the **one-shot Möbius amplifier**.

## 1. Three commuting collision directions

On `T=R/Z`, write

```text
P_a h(y)=(1/a)sum_(r=0)^(a-1)h((y+r)/a),            (3)

(P_a h)_hat(n)=h_hat(an),
P_aP_b=P_(ab)=P_bP_a.                               (4)
```

Let `A,B` be nonzero nonnegative real step functions with

```text
AB=0 almost everywhere,                             (5)
```

and let `J_A,J_B` be their numbers of nonzero jumps. Fix a prime `ell`
different from `7,13`. For `a,b,d>=0`, put

```text
A_(a,b,d)=P_13^a P_7^b P_ell^d A,
B_(a,b,d)=P_13^a P_7^b P_ell^d B,

I_(a,b,d)=integral_T A_(a,b,d)B_(a,b,d),            (6)

Omega={x in N^3:I_x>0}.                             (7)
```

As in THM-2313, `Omega` is an upper set. If two nonnegative functions
`U,V` have a common-positive set of positive measure, the image of that set
under `x -> px mod 1` still has positive measure. Every image point has a
preimage at which both functions are positive, so `P_pU` and `P_pV`
overlap positively. This applies separately to `p=13,7,ell`.

Each coordinate axis eventually enters `Omega`, because

```text
||P_p^r h-integral h||_2 ->0                       (8)
```

for every `h in L^2(T)`. Thus the Pareto frontier is finite. The theorem
needs the following particularly useful frontier shape.

> **One-shot corner hypothesis.** For some `s>=1`, the cell `(s,1,1)` is
> Pareto-minimal in `Omega`.

Since `Omega` is upper, it is enough to check the three codimension-one
predecessors:

```text
I_(s,1,1)>0,

I_(s-1,1,1)=I_(s,0,1)=I_(s,1,0)=0.                 (9)
```

Every other vertex in the lower cube is dominated by one of the three
faces and therefore also has zero collision.

## 2. The full cube derivative

Parseval gives the absolutely convergent series

```text
I_(a,b,d)
 =sum_(n in Z)
   A_hat(13^a 7^b ell^d n)
   conjugate(B_hat(13^a 7^b ell^d n)).              (10)
```

At the predecessor scale `13^(s-1)`, expand the three-prime selector:

```text
[gcd(n,13*7*ell)=1]

 =(1-[13|n])(1-[7|n])(1-[ell|n])

 =1-[13|n]-[7|n]-[ell|n]
   +[91|n]+[13ell|n]+[7ell|n]
   -[91ell|n].                                      (11)
```

The corresponding eight collision terms are

```text
 I_(s-1,0,0)
-I_(s,  0,0)
-I_(s-1,1,0)
-I_(s-1,0,1)
+I_(s,  1,0)
+I_(s,  0,1)
+I_(s-1,1,1)
-I_(s,  1,1).                                      (12)
```

Under (9), the first seven terms vanish. Hence

```text
sum_(gcd(n,91ell)=1)
 A_hat(13^(s-1)n)
 conjugate(B_hat(13^(s-1)n))

 =-I_(s,1,1)<0.                                    (13)
```

The negative sign is the parity of a three-dimensional corner. More
important than the sign is the scale:

```text
seven and ell occur in the selector (11),
but neither occurs in the Fourier argument in (13).                  (14)
```

This is not obtainable from a two-dimensional `(s,t)` corner with `t>1`.

## 3. Endpoint-Prony landing

Put

```text
U=P_13^(s-1)A,
V=P_13^(s-1)B.                                      (15)
```

They have at most `J_A,J_B` jumps. For `q>0`,

```text
C_q=(2*pi*q)^2 U_hat(q)conjugate(V_hat(q))          (16)
```

is an exponential sum on at most

```text
L<=J_AJ_B                                           (17)
```

endpoint-difference nodes. Equation (13), together with conjugation for
negative frequencies, shows that at least one positive unit-residue
progression

```text
q=r+(91ell)h,             gcd(r,91ell)=1            (18)
```

is not the zero sequence. A nonzero `L`-node exponential sum cannot vanish
at `L` consecutive indices. Testing `h=0,...,L-1` supplies a positive
integer `n` such that

```text
gcd(n,91ell)=1,
1<=n<=(91ell)J_AJ_B-1,                              (19)

A_hat(13^(s-1)n)!=0,
B_hat(13^(s-1)n)!=0.                                (20)
```

The ceiling in (19) is exact bookkeeping:

```text
(91ell-1)+(91ell)(L-1)=91ell L-1.                  (21)
```

Now use THM-2306's owner-normalized source/current pair

```text
f=1_QP_13^k1_E,
A=P_c f,
B=P_c1_E,                                           (22)

J_A<=4S,
J_B<=2S.                                           (23)
```

Equations (19)--(20) become a common original source/current frequency

```text
N=c13^(s-1)n,                                      (24)

f_hat(N)!=0,
(1_E)_hat(N)!=0,                                   (25)

gcd(n,91ell)=1,
n<=8(91ell)S^2-1.                                  (26)
```

If

```text
c=13^lambda 7^gamma w,
gcd(w,91)=1,                                       (27)
```

then

```text
nu_13(N)=lambda+s-1,
nu_7(N)=gamma.                                     (28)
```

Unlike THM-2313's interior `(s,t)` frequency, (24) adds no seven-adic
valuation. The auxiliary primes are finite-difference coordinates, not
permanent clock factors.

## 4. Exact grade-three amplifier on the hostile carrier

Use THM-2299's exact multiplier-four local carrier and THM-2306's normalized
supports:

```text
B=(1/13)1_F,
F centers {1,15}/16, half-width epsilon;

A=(1/2197)1_C,
C centers {7,9}/16, half-width 169epsilon;

epsilon=10^(-12).                                  (29)
```

Take the Fermat prime

```text
ell=65537=2^16+1.                                  (30)
```

Direct trial division through `256` verifies primality. Moreover

```text
ell congruent 1 (mod 16),
gcd(ell,91)=1.                                     (31)
```

For a scale

```text
M=13^a 7^b ell^d,
```

with `b,d in {0,1}`, multiplication by `ell` fixes both center sets modulo
sixteen, while multiplication by seven swaps the source and current center
patterns up to signs. The cross-gap is therefore

```text
g_a=3/8 if a is even,
g_a=1/8 if a is odd.                               (32)
```

The sum of the pushed support radii is `170M epsilon`; at all four
load-bearing cells neither interval wraps. Positive overlap is exactly

```text
170M epsilon>g_a.                                  (33)
```

The top cell and its three faces have the following exact integer ledger.
Multiplying (33) by `8*10^12` gives the displayed comparison.

| cell | scale `M` | left side `1360M` | right side | collision |
|:---|---:|---:|---:|:---:|
| `(3,1,1)` | `1007893523` | `1370735191280` | `1000000000000` | yes |
| `(2,1,1)` | `77530271` | `105441168560` | `3000000000000` | no |
| `(3,0,1)` | `143984789` | `195819313040` | `1000000000000` | no |
| `(3,1,0)` | `15379` | `20915440` | `1000000000000` | no |

Thus `(3,1,1)` is fully Pareto-minimal. In fact it is the only positive
cell in the complete target box

```text
0<=a<=3,        0<=b,d<=1.                         (34)
```

The owner is `c=13`. Equations (24)--(26) therefore give

```text
N=13^3 n,
gcd(n,13*7*65537)=1.                               (35)
```

In particular,

```text
nu_13(N)=3,
nu_7(N)=0,
gcd(N/13^3,91)=1.                                  (36)
```

This is exactly the simultaneous arithmetic that the two-dimensional
hostile frontier did not contain:

```text
THM-2313 corner (3,7):
  correct grade, but added 7^6;

THM-2313 corner (8,1):
  no added seven, but wrong grade;

THM-2318 corner (3,1,1):
  correct grade and a whole 91-unit outside multiplier.              (37)
```

Both normalized hostile supports have four jumps. Since

```text
91*65537=5963867,
```

the local endpoint-Prony bound is

```text
5963867*4*4-1=95421871.                             (38)
```

The global THM-2306 form is

```text
n<=47710936 S^2-1.                                 (39)
```

The larger numerical bank is the cost of selecting a complete
three-prime unit class without truncating Fourier tails.

## 5. What the amplifier does and does not land

For the hostile profile `(1,3,5)`, every THM-2293 shell-graph vertex has
thirteen-adic grade three. Equation (36) now places the exact
word-current/source atom at that same grade and gives it a residual
multiplier coprime to `91`.

This closes a genuine representational loss:

```text
prescribed clock + exact word current
  -> correct shell grade
  -> no added factor of seven
  -> complete 91-unit residual.                    (40)
```

It does **not** identify the current atom with a shell-graph vertex or edge.
Three independent coordinates remain.

1. **Root residue.** Equation (35) only says

   ```text
   n mod 13 is nonzero.
   ```

   It does not force THM-2293's selected character `kappa`. The hostile
   pair character `4`, for example, would require `n congruent 4 mod 13`,
   which (13) does not select.

2. **Pair-edge incidence.** The current coefficient is

   ```text
   (1_QP^k1_E)_hat(N),
   ```

   while a THM-2293 edge consists of two bare source coefficients whose
   difference is `m c_3`. A unit outside multiplier of the correct
   arithmetic type is not yet the proved relation multiplier `m`.

3. **Target gain and phase.** The provisional marked-target analysis in
   THM-2315 separates a terminal word's support fibre from its exact mixed
   gain. Independently, audited THM-2303 and THM-2299 show that component
   base phase can cancel. The real scalar collision cube retains neither
   coordinate.

There is also no universal claim that every LRC word current has a corner
`(b-lambda+1,1,1)` for this fixed `ell`. The theorem proves the conditional
one-shot mechanism and an exact hostile positive control. A future global
use must either force such a corner, choose an amplifier prime from the
actual component geometry, or replace positivity by a lawful
residue-sensitive sidecar.

## 6. Connection and loss ledger

```text
source:
  THM-2313's collision antichain and three-prime divisor-lattice hint;

new operation:
  use seven and ell once each, so their product creates the top collision
  while both remain at exponent zero in the predecessor shell;

preserved:
  source owner, original source coefficient, exact word-current
  coefficient, exact thirteen grade, unchanged seven valuation, and a
  residual core coprime to 91ell;

new exact object:
  a fully minimal Boolean collision cube with odd Möbius sign;

destroyed or unselected:
  prescribed root residue, pair-edge incidence, exact target gain,
  component base phase, and a global proof that every live row has the
  required one-shot corner;

cheapest decisive probe:
  ell=65537 on THM-2299's exact two-interval carrier, where the top and all
  three faces are separated by large rational margins.               (41)
```

The structural lesson is broader than this carrier: an auxiliary prime can
be used as a **difference direction** without becoming a valuation of the
landed atom. This is the correct way to add mixing power when the target
requires the auxiliary prime to stay absent from the final multiplier.

No scalar profile is eliminated, the hostile carrier remains local rather
than a global LRC cover, and LRC(14) remains open.

## 7. Exact verification

The companion uses only integer and `Fraction` arithmetic. It exhausts all
twenty upper subsets of the Boolean three-cube, verifies the eight-term
Möbius selector and odd-corner sign, proves the exact endpoint-Prony
ceiling, trial-divides `65537`, checks every cell in (34), records the four
strict integer inequalities above, and audits `7912` sample landed
frequencies for exact thirteen/seven valuations and `91`-unit residuals.

Reproduce with

```bash
python3 04-computation/lrc14_one_shot_mobius_amplifier_thm2318.py
python3 -O 04-computation/lrc14_one_shot_mobius_amplifier_thm2318.py
```

Every load-bearing check raises explicitly, so optimized mode executes the
same audit. QED.
