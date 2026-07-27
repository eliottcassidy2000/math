---
id: THM-2581
title: "The b-word depth-five owner-clock host and reflection-breaking mixedness"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On the canonical typed row, wholly inside the lawful sigma={b}, K=2,
  r=5 first-collision fibre, the common-base collision-root/owner-clock
  array satisfies both zero-margin laws after centring, is nonzero in all
  91 cells, has maximal centred column rank 6, and has all 1,638 cyclotomic
  2x2 minors nonzero.  Its reflected central-cell difference is a nonzero
  odd drift with every nontrivial F_13 Fourier colour.  This puts the owner
  clock and the theta=t-2u sidecar on one word fibre, but does not perform
  that contraction, construct a temporal intertwiner or physical target
  current, exclude a row, or prove LRC(14).
source: b-r5-owner-clock-2026-07-28
depends_on:
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2575-word-depth-collision-law-and-owner-clock-host-array
  - THM-2577-all-clock-word-depth-collision-law-and-support-image-mechanism
related:
  - THM-2512-lawful-interaction-cut-bundle-transplant-and-replica-dichotomy
  - HYP-9032-the-transplant-trichotomy-rehoming-the-91-stalk-laws
script: 04-computation/lrc14_b_r5_owner_clock_host_thm2581.py
output: 05-knowledge/results/lrc14_b_r5_owner_clock_host_thm2581.out
script_sha256: 38c8eac0925d119331ef248fa7197f61bd985166c56704575cec18c361543647
output_sha256: 429dbc796f73b35a03d9223b09e089d83049d8fadf3231ba017df5490f41dcb9
engine_sha256: 99deb8891f8be1e218e7b4c5a7bad7a918ddc1da37b05b53ccbbc3ca7b8c7998
hash_basis: working-tree bytes (LF)
---

# THM-2581 -- the owner clock lives on the depth-five `b` fibre

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2575 produced two adjacent facts on different word fibres: a
`13 x 7` owner-clock host on the `sigma={a}`, `r=3` packet and an open
`theta=t-2u` sidecar on the `sigma={b}`, `r=5` packet.  THM-2577 proved that
the latter depth is an all-clock support law, but warned that multiplying the
two facts would repeat the common-base error of MISTAKE-281.

The host can instead be rebuilt directly on the `b` packet.  It is stronger
there: its former central reflection degeneracy becomes an oriented odd
drift, every `2 x 2` minor is nonzero, and the centred array has the largest
rank compatible with its owner-clock zero sum.

## 1. The same-fibre object

Retain the canonical typed row and licensed owner

```text
w=(1,14,27,40,53,66,13,13^3,2*13^5),      j=1.               (1)
```

Let `E=E_1`, let `Q=Q_{1,{b}}`, and take the canonical packet clock `K=2`:

```text
f=1_Q P_(13^2)1_E,

d=13^5,

U=P_d f,                         V=P_d1_E.                    (2)
```

THM-2577 gives the first-collision index `r=5`; THM-2471 (47) gives the
unique affine deep-root invariant

```text
theta=t-2u mod 13.                                             (3)
```

Here the invariant has a particularly concrete form.  In THM-2471's Boolean
stalk coordinate

```text
X_(u,a)=(y+u+13a)/(13d),
```

the deepest speed is `C=2d`, so

```text
C X_(u,a)=2(y+u)/13                         modulo 1.          (3a)
```

Write `z={2y}` and `epsilon=floor(2y) in {0,1}`.  In the oriented deep-root
chart `(z+t)/13`, equation (3a) says

```text
t=2u+epsilon mod 13,              theta=epsilon.              (3b)
```

Thus the affine sidecar is an endogenous half-circle label, not a declared
gauge.  Its transformation laws also expose its exact scope:

```text
(u,t)->(u+h,t+2h):          theta fixed,

t->t+q with u fixed:        theta->theta+q,

(u,t)->(-u-1,-t-1):         theta->1-theta.                   (3c)
```

The last line is oriented chart reflection away from the null chart
boundaries; equivalently `2theta-1` changes sign.  In particular, this gives
an absolute binary reference but not a target-translation-invariant full
`F_13` reference.  The theorem below does not yet insert the signed
half-circle weight into the host contraction.

Disintegrate the collision root exactly as in THM-2471:

```text
A(y,u)=U((y+u)/13),           F(y,u)=V((y+u)/13),

c_s(y)=sum_(u in F_13) A(y,u+s)F(y,u).                        (4)
```

The native owner-clock cells from THM-2449 are

```text
H_ell={y: frac(13y) in [(2ell-1)/14,(2ell+1)/14) mod 1},
                                                       ell in F_7. (5)
```

Define

```text
C_ell(s)=integral_(H_ell)c_s(y)dy,

W(k,ell)=1/169 sum_(s in F_13)C_ell(s)zeta_13^(-ks).           (6)
```

This order is load-bearing: the owner-cell indicator sits **inside the same
`y`-integral** as the collision-root density, before the root DFT and before
either marginal.  Thus (6) is a common-base object on the `b`, `r=5` fibre;
it is not a product of the earlier `a`-word host with the later depth law.

## 2. Exact reconstruction and first-collision gates

All sets and profiles are reconstructed on the half-open grid

```text
T_DEN=182 lcm(w)=297836897838480.                              (7)
```

The input anchors are

```text
#intervals(E)=57072,      measure(E)=1882176/28589561,

#intervals(Q)=131762,     measure(Q)=4839079319/190921088358,

integral f=35505957232/16132831966251.                         (8)
```

The exact engine computes (6) twice:

1. extract all thirteen root-chart profiles in `y`, multiply each ordered
   pair, and integrate it separately over the seven cells;
2. use the independent substitution `x=(y+u)/13` to compute

   ```text
   C_ell(s)=13 integral_T U(x)V(x-s/13)
       1_{frac(169x) in [(2ell-1)/14,(2ell+1)/14)} dx.         (9)
   ```

The two routes agree on all `91` cells and all `13` global row sums.  A
separate hostile toy with `p=3`, five clock cells, wraps, and unequal pieces
agrees with definitional midpoint integration on all `15` entries.

Every `C_ell(s)` is an integer over the common denominator

```text
D_C=169 d^2 T_DEN
   =6939029398456584394954868880
   =2^4 3^3 5 7^2 11 13^18 53.                               (10)
```

The exact first-collision laws are

```text
C_ell(0)=0                              for every ell,

sum_(s,ell) C_ell(s)/D_C
 =169 I_5,

I_5=48602521488933856/337437093630814766589.                  (11)
```

The first line is the cellwise form of the pre-collision disjointness; the
second reproduces the independently stored THM-2575 collision value.

## 3. Margins, saturation, and maximal rank

Put

```text
J(k)=sum_ell W(k,ell),

W^c(k,ell)=W(k,ell)-J(k)/7.                                  (12)
```

Then exactly

```text
sum_ell W(k,ell)=J(k),                 J(0)=I_5,

sum_k W(k,ell)=C_ell(0)/13=0,

W(-k,ell)=conj(W(k,ell)).                                    (13)
```

All twelve `J(k)`, `k!=0`, are nonzero.  In every one of the seven owner
cells all twelve nonzero root colours survive.  Moreover

```text
W^c(k,ell)!=0                         for all 91 pairs,

sum_ell W^c(k,ell)=0                  for every k,

sum_k W^c(k,ell)=0                    for every ell.           (14)
```

Thus `W^c` satisfies HYP-9032's (T1), and its `k`-independence branch is
nonvacuously excluded as in (T2).  Its exact coordinate floor is

```text
D=184240653122424862557594
 =2*3^2*7^3*11*13^15*53,

min nonzero absolute zeta-basis coordinate
 =76/11823941286254964867 >=1/D.                              (15)
```

This supplies (T3) in the stated cyclotomic basis.  Over `Q(zeta_13)`,

```text
rank W=7,

rank W^c=6.                                                   (16)
```

The second rank is maximal because the first identity in (14) gives one
column relation.  Reduction at the primitive thirteenth root `18 mod 79`
already gives ranks `7` and `6`, so (16) is also independently certified.

## 4. The reflection law and the oriented drift

The sets `E,Q` are invariant under `x -> -x`, hence `U,V` are even.  The
change `(y,u)->(-y,-u)` in (4), together with `H_ell -> H_{-ell}`, proves

```text
C_ell(s)=C_{-ell}(-s).                                       (17)
```

The two central reflected owner cells are `ell=3,4`.  Define their oriented
difference

```text
Delta(s)=C_4(s)-C_3(s).                                      (18)
```

Equation (17) makes `Delta` odd.  In integer numerator coordinates over
`D_C`, its exact functional form is

```text
(0,
 2371250083334441400,
-2372727465429080760,
 0,0,
 7537663748160,
 13943390382112070880,
-13943390382112070880,
-7537663748160,
 0,0,
 2372727465429080760,
-2371250083334441400).                                      (19)
```

In particular,

```text
Delta(1)=47817205/139928299245620886 !=0.                    (20)
```

Its mean vanishes and every nonzero Fourier colour is nonzero:

```text
sum_s Delta(s)=0,

sum_s Delta(s)zeta_13^(-ks)!=0       for every k!=0.          (21)
```

Since

```text
W(k,4)-W(k,3)=1/169 sum_s Delta(s)zeta_13^(-ks),              (22)
```

the orientation lost by summing the owner clock is retained in every
nontrivial collision-root colour.  This is the array-level odd sidecar
compatible with the `theta=t-2u` branch (3).  It does **not** by itself
evaluate the THM-2512 cut contraction or identify a target-active current.

The separate `sigma={a}`, `r=3` fibre is a sharp hostile control.  Replaying
the identical engine there gives

```text
Delta_a(s)=0 for every s,

C^a_3(s)=C^a_4(s),                                           (23)
```

which is exactly the earlier central-column degeneracy.  Thus the drift in
(19) is not a generic artifact of half-open cells, the root DFT, or the
integration engine; it appears after moving the whole host to the lawful
depth-five word fibre.

## 5. Complete minor certificate

For all choices

```text
0<=k_1<k_2<=12,       0<=ell_1<ell_2<=6,                     (24)
```

the determinant

```text
W(k_1,ell_1)W(k_2,ell_2)
-W(k_1,ell_2)W(k_2,ell_1)                                   (25)
```

is nonzero in `Q(zeta_13)`.  This is exactly

```text
C(13,2) C(7,2)=1638/1638                                    (26)
```

nonzero minors.  Two independent exact certificates are stored:

- canonical reduction of every integer numerator polynomial modulo
  `Phi_13` finds no zero;
- evaluation at primitive thirteenth roots in `F_79` and `F_131` proves
  nonvanishing without using cyclotomic canonicalization.  The first prime
  certifies `1625` minors; the second certifies the remaining `13`.

Neither prime divides the common denominator.  The hostile `a`, `r=3`
replay instead has exactly `90` zero minors, all accounted for by (23).

At the two other finitely censused lawful clocks `K=3,4`, the `b`, `r=5`
host again has `1638/1638` nonzero minors.  This is a robustness control, not
an all-clock mixedness theorem; the theorem quantifier remains `K=2`.

## 6. Consequence and sharp boundary

The former word-fibre mismatch is gone on this typed row:

```text
sigma={b}, K=2, r=5 packet
  -> common collision-root service
  -> native owner-clock host W
  -> maximal centred F_7 rank and every root-colour drift
  -> theta=t-2u available on the same packet.                 (27)
```

What remains is substantive.  The array `W^c` is cyclotomic-valued rather
than the rational THM-2512 defect `d_A`; its nonvanishing is not yet a
target-active Boolean event.  No theorem here performs the slaved
`theta=t-2u` contraction, constructs the off-diagonal temporal intertwiner,
identifies semantic endpoint arrival, transports the construction over the
other `164` rows, or excludes even this typed row.  The LRC(14) ledger remains
`165`.

## 7. Exact companion

Run

```bash
python3 04-computation/lrc14_b_r5_owner_clock_host_thm2581.py
python3 -O 04-computation/lrc14_b_r5_owner_clock_host_thm2581.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_b_r5_owner_clock_host_thm2581.out.
```

The companion imports only the previous exact route-doubled interval/profile engine
from `lrc14_base_only_bridge_opus_20260728.py`; its hash is recorded in the
frontmatter.  It rebuilds both word fibres from their Boolean combs and runs
all positive and hostile controls described above.

The independent audit rederived the common-base normalization and both
marginals, checked the `theta=floor(2y)` law and its gauge/reflection action,
and verified that the modular minor cover proves nonvanishing in
`Q(zeta_13)`.  It replayed the ordinary and optimized scripts, stored-output
comparison, hashes, and bytecode compilation.  The auditor confirmed that
the theorem stays on the `{b}`, `K=2`, `r=5` fibre and does not identify the
retained root with THM-2365's still-missing quotient-dual target endpoint.

**QED.**
