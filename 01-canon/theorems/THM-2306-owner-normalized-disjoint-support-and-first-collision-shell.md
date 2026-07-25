---
id: THM-2306
title: "Owner-normalized disjoint support and the first collision shell"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. On every
  positive canonical blocker-word handoff at prescribed expiration,
  multiply the exact arrival density by its current word and normalize both
  that density and the original exclusive source by the source-owner comb.
  The two normalized functions lie on opposite sides of D_1. Their exact
  negative Fourier covariance forces one frequency shared by the original
  source and the word-restricted current at an owner multiple
  c_j q, with 1<=q<=8S^2-1. More sharply, push both normalized densities
  until their supports first meet with positive mass. The preceding
  collision-free scale gives a negative covariance on one exact
  thirteen-adic shell, hence a shared atom
  c_j 13^(r-1)n with 13 not dividing n and
  1<=n<=104S^2-1. Its valuation is exactly lambda_j+r-1. The collision
  depth is finite and has an explicit BV/mass bound, but is not uniformly
  one. The current coefficient belongs to the product arrival density, not
  to the bare Perron arrival, and the owner-multiple bank is not the
  prescribed shallow-pair subgroup. No scalar row is excluded and LRC(14)
  remains open.
source: codex-2026-07-25-owner-normalized-first-collision
depends_on:
  - THM-2296-prescribed-expiration-return-or-bounded-ancestry-resonance
  - THM-2305-canonical-blocker-word-handoff-hypergraph
related:
  - THM-2267-static-owner-coverage-is-flag-and-transition-holonomy-is-a-cut-kernel
  - THM-2276-shallow-owner-residue-aligned-crossing
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2302-same-label-expiration-dichotomy-and-pure-terminal-shell-no-go
  - THM-2303-terminal-component-phase-current-and-defect-rank
script: 04-computation/lrc14_owner_normalized_first_collision_thm2306.py
output: 05-knowledge/results/lrc14_owner_normalized_first_collision_thm2306.out
script_sha256: 5cbaf6c1a4c7691eceacaac63b3ee26b63bc4776527f4a433c02a341048f4c80
output_sha256: 7fc14ce95f487d60fe1b84c5732605904fa0f5c55660ecffee0fb8220ca3b9a2
hash_basis: working-tree bytes (LF)
---

# THM-2306 -- normalize an edge by its owner, then wait for first collision

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The blocker-word graph in THM-2305 has an honest pairwise observable on a
pure edge and an honest directed hyperedge on a fork. It still cannot compose
two edges merely because their owner labels agree: the incoming and outgoing
subsets of the middle owner can be disjoint.

There is, however, a canonical operation on every single positive handoff.
Normalize both ends by the **source owner itself**. The original source then
lives in `D_1`, while the exact word-restricted arrival lives in its
complement. This gives a signed identity, not just energy:

```text
positive current word
  -> owner-normalized source and current have disjoint support
  -> their zero-frequency product must be cancelled by nonzero modes
  -> one exact owner-multiple source/current atom is shared.             (1)
```

Pushing the two normalized densities by thirteen eventually makes them
overlap. At the first scale where this happens, the negative covariance is
concentrated on one exact thirteen-adic valuation shell. This first
collision depth is the lawful gluing sidecar that the bare owner graph lacks.

## 1. The exact word-restricted arrival

Use THM-2296 and THM-2305. The selected exclusive owner is `j`, with

```text
E=E_j,
c=c_j=13^lambda u,              13 does not divide u,
k=lambda+1,
P=P_13.                                                     (2)
```

Let `Q=Q_(j,sigma)` be any one of the three canonical terminal words

```text
sigma in {{a},{b},{a,b}}
```

having positive prescribed return

```text
rho_Q
 :=integral_Q P^k 1_E
  =measure(E intersection T^(-k)Q)
  >0.                                                       (3)
```

On THM-2305's selected large word in the quantitative return branch,

```text
rho_Q>=rho_j/3,

strict:
 rho_Q>=39002430583/160481782761300;

repeated-first:
 rho_Q>=13560199813/160481782761300.               (4)
```

Define the exact **word-restricted arrival density**

```text
f=1_Q P^k1_E.                                       (5)
```

The product in (5) is load-bearing. It simultaneously retains:

- ancestry multiplicity at the prescribed clock;
- the current point, rather than a pushed target;
- the complete pure or double blocker word `sigma`.

It has mass

```text
integral f=rho_Q.                                   (6)
```

It is not the bare arrival `P^k1_E`, and a Fourier atom of `f` need not be
an atom of either factor in (5).

Because every terminal word lies in

```text
R_j=A_0 minus D_c,
```

the source and current satisfy

```text
support(1_E) subset D_c,
support(f) subset D_c^c                             (7)
```

up to null endpoints.

For every positive integer `a`, write

```text
P_a h(y)=(1/a)sum_(r=0)^(a-1)h((y+r)/a).            (8)
```

It preserves mass and obeys

```text
(P_a h)_hat(n)=h_hat(an).                           (9)
```

Normalize both ends of the handoff by the owner:

```text
A=P_c f,
B=P_c1_E.                                           (10)
```

If `A(y)>0`, some preimage `x` with `cx=y mod 1` belongs to
`support(f)`, so `y` lies outside `D_1`. The same argument on the source
puts `support(B)` inside `D_1`. Therefore

```text
support(A) subset D_1^c,
support(B) subset D_1,
AB=0 almost everywhere,                             (11)

integral A=rho_Q,
integral B=measure(E)>0.                             (12)
```

Both functions are nonzero nonnegative rational step functions.

## 2. Disjointness forces a bounded common owner multiple

Parseval applied to (11) gives the exact signed covariance

```text
sum_(n!=0) f_hat(cn)conjugate((1_E)_hat(cn))
 =-rho_Q measure(E)
 <0.                                                (13)
```

Indeed, (9) identifies the two Fourier factors with `A_hat(n)` and
`B_hat(n)`, while the omitted zero term is the positive product in (12).
The series is absolutely convergent by Cauchy--Schwarz.

For THM-2305's selected word, the magnitude of (13) has the exact branch
floors

```text
strict:
 rho_Q measure(E)
 >=586652368446484273/95291384904891722547000;

repeated-first:
 rho_Q measure(E)
 >=70913620890275833/95291384904891722547000.       (13a)
```

These are `alpha delta/3`, using the owner and return floors recorded in
THM-2296. Thus (13) is quantitatively nonzero before any individual
frequency is selected.

Consequently some nonzero integer `n` satisfies

```text
f_hat(cn)!=0,
(1_E)_hat(cn)!=0.                                   (14)
```

This already couples the original source and exact current word at the same
frequency. It can be landed uniformly.

Let `J_A,J_B` denote the numbers of nonzero jumps of `A,B`. The bilinear
endpoint-Prony lemma from THM-2296, applied to `A,B`, changes (14) to a
positive integer `q` with

```text
1<=q<=J_AJ_B-1,

f_hat(cq)!=0,
(1_E)_hat(cq)!=0.                                   (15)
```

For completeness, the endpoint product

```text
C_n=(2*pi*n)^2 A_hat(n)conjugate(B_hat(n))           (16)
```

is an exponential sum on at most `J_AJ_B` endpoint-difference nodes.
It is not the zero sequence by (13), while its formal value at zero is the
product of the two total jump sums and hence vanishes. A nonzero `L`-node
exponential sum cannot vanish at `L` consecutive indices, so one of
`1,...,L-1` is nonzero. This proves (15).

Use

```text
S=H+sum_i q_i+sum_h c_h.                            (17)
```

The full source and every canonical word have at most `2S` jumps.
Perron transport cannot create jump locations, and multiplication of two
step functions creates jumps only in the union of their jump sets. Hence

```text
J_B<=2S,

J_A
 <=#Jump(f)
 <=#Jump(P^k1_E)+#Jump(1_Q)
 <=4S.                                              (18)
```

Substitution into (15) gives

```text
1<=q<=8S^2-1.                                       (19)
```

This conclusion applies to a pure word and to a fork without orienting the
fork.

## 3. The first positive collision

The atom in (15) can be ramified by an unknown power of thirteen. The exact
coordinate that determines that power is geometric.

For `s>=0`, put

```text
A_s=P^s A,
B_s=P^s B,
I_s=integral A_s B_s.                               (20)
```

Every `I_s` is nonnegative and `I_0=0` by (11). Define

```text
r=min{s>=1:I_s>0}.                                  (21)
```

This integer exists. For any `L^2` function `h`,

```text
||P^s h-integral h||_2^2
 =sum_(n!=0)|h_hat(13^s n)|^2
 ->0                                                (22)
```

because the Fourier coefficients are square summable. Thus

```text
I_s->(integral A)(integral B)
     =rho_Q measure(E)>0.                           (23)
```

In support language, `r` is the first depth at which the two nonnegative
normalized image densities overlap on positive Haar mass. It retains more
than the assertion that two owner labels agree: it is an actual common
spatial scale.

Fourier expansion of (20) is absolutely convergent and gives

```text
I_s
 =sum_(n in Z)
    A_hat(13^s n)conjugate(B_hat(13^s n)).           (24)
```

At every scale, subtracting the multiples of thirteen from the full lattice
gives the complete shell ledger

```text
sum_(13 does not divide n)
 f_hat(c 13^s n)
 conjugate((1_E)_hat(c 13^s n))
 =I_s-I_(s+1).                                      (24a)
```

Minimality says `I_s=0` for `0<=s<=r-1` and `I_r>0`. Consequently every
aggregate shell before `r-1` vanishes:

```text
I_s-I_(s+1)=0,                         0<=s<r-1,     (24b)
```

whereas the shell at `r-1` is strictly negative:

```text
sum_(13 does not divide n)
 A_hat(13^(r-1)n)
 conjugate(B_hat(13^(r-1)n))
 =I_(r-1)-I_r
 =-I_r
 <0.                                                (25)
```

Using (9), this is equivalently

```text
sum_(13 does not divide n)
 f_hat(c 13^(r-1)n)
 conjugate((1_E)_hat(c 13^(r-1)n))
 =-I_r<0.                                           (26)
```

Thus the covariance is not merely nonzero somewhere in the owner subgroup.
It is nonzero on the one exact valuation shell selected by first spatial
collision. More precisely, `r-1` is the first nonzero **aggregate**
thirteen-adic covariance shell. Individual products can occur and cancel
inside an earlier shell; (24b) does not assert their termwise vanishing.

## 4. Bounded landing in a nonzero root residue

Equation (25) supplies some `n` prime to thirteen with both factors
nonzero. A residuewise endpoint-Prony argument makes it finite.

Put

```text
U=A_(r-1),
V=B_(r-1).                                          (27)
```

They have at most `J_A,J_B` jumps. For `a in {1,...,12}`, restrict (16) to

```text
q=a+13 ell.
```

After combining equal nodes,

```text
(2*pi*q)^2 U_hat(q)conjugate(V_hat(q))
```

is an exponential sum in `ell` on at most

```text
L<=J_AJ_B                                            (28)
```

nodes. Equation (25), and conjugation if necessary, says that at least one
positive residue progression is not the zero sequence. It cannot vanish at
all `L` consecutive indices

```text
ell=0,1,...,L-1.
```

Therefore there are positive integers `n,N` with

```text
13 does not divide n,
1<=n<=13J_AJ_B-1<=104S^2-1,

N=c 13^(r-1)n,                                      (29)

f_hat(N)!=0,
(1_E)_hat(N)!=0.                                    (30)
```

Since `c=13^lambda u` and both `u,n` are thirteen-units,

```text
nu_13(N)=lambda+r-1.                                (31)
```

The bound in (29) is on the outside unit `n`; the collision depth records
the complete ramification separately.

## 5. A quantitative collision-depth bound

The depth `r` is not just existential. There are two complementary
quantitative bounds.

First use literal support components. A nonzero step function with `J`
jumps has at most `J` positive support components; this deliberately
conservative convention also covers positive-to-positive jumps. Since

```text
0<=A,B<=1,
```

some positive component of `A` has length at least `rho_Q/J_A`, and some
positive component of `B` has length at least `measure(E)/J_B`. Once
thirteen-fold expansion makes either interval have length at least one, the
corresponding Perron image is positive almost everywhere on the circle and
must overlap the other nonzero density. Hence

```text
r<=min(r_A^comp,r_B^comp),                          (32)

r_A^comp=min{s>=1:13^s rho_Q>=J_A},
r_B^comp=min{s>=1:13^s measure(E)>=J_B}.
```

Using (18), a convenient sufficient condition is

```text
13^s rho_Q>=4S
or
13^s measure(E)>=2S.                                (33)
```

In particular, any `s>=1` satisfying both

```text
13^s delta>=12S,
13^s alpha>=2S                                      (34)
```

bounds `r` on the corresponding branch. Only one inequality is actually
needed when the sharper live mass is retained.

Second, for a periodic bounded-variation function,

```text
Var(P_a h)<=Var(h)/a,
||h-integral h||_infinity<=Var(h).                  (35)
```

Hence `P^s h` is strictly positive almost everywhere as soon as

```text
13^s>Var(h)/integral h.                             (36)
```

If either `A_s` or `B_s` is everywhere positive, its product with the other
nonzero nonnegative density has positive integral. Therefore

```text
r<=min(r_A,r_B),                                    (37)

r_A=min{s>=1:13^s>Var(A)/rho_Q},
r_B=min{s>=1:13^s>Var(B)/measure(E)}.
```

The LRC jump ledger gives the more explicit variation bounds

```text
Var(B)<=2S/c,

Var(A)
 <=[Var(P^k1_E)+Var(1_Q)]/c
 <=[2S/13^k+2S]/c.                                  (38)
```

Let `alpha` be the selected-owner mass floor from THM-2296 and `delta` its
return floor:

```text
strict:
 alpha=15041431/593783190,
 delta=39002430583/53493927587100;

repeated-first:
 alpha=5229541/593783190,
 delta=13560199813/53493927587100.                  (39)
```

For THM-2305's selected word, `measure(E)>=alpha` and
`rho_Q>=delta/3`. One fully explicit branch-uniform form of (37) is

```text
r<=min{s>=1:
 13^s>
 min(
   6S(1+13^(-k))/(c delta),
   2S/(c alpha)
 )
}.                                                  (40)
```

This keeps the source mass, return mass, owner size, and variation as
separate coordinates rather than hiding them in a universal spectral bank.

## 6. Exact stopping boundary

The theorem closes one genuine gap:

```text
source owner + prescribed clock + exact current word
  -> same-frequency source/current atom
  -> exact first-collision valuation shell.                         (41)
```

It does not close the following distinct gaps.

1. **Product arrival versus bare arrival.** The current coefficient is

   ```text
   f_hat(N)=(1_Q P^k1_E)_hat(N),
   ```

   a convolution in Fourier space. It need not make
   `(P^k1_E)_hat(N)` or `(1_Q)_hat(N)` nonzero.

2. **Owner subgroup versus pair subgroup.** Equation (29) lands in `c Z`.
   If THM-2276's prescribed pair carry is `K=m c`, its relation subgroup is
   `m c Z`. Nothing here forces `m` to divide `13^(r-1)n`, and nothing
   selects the base atom `K`.

3. **Collision versus orbit composition.** The spatial overlap in (21)
   occurs after pushing two owner-normalized densities by the same number
   of root steps. It is not an identification of the incoming and outgoing
   subsets in a THM-2305 owner cycle.

4. **All-clock use and the exact shell-grade gate.** Nothing in the
   disjoint-support proof uses `k=lambda+1`; it uses only

   ```text
   E subset D_c,
   Q subset D_c^c,
   integral_Q P^k1_E>0.                              (41a)
   ```

   It therefore applies verbatim to THM-2302's same-label positive-return
   arm at any fixed clock. At its common shell clock `k=b+1`, a THM-2306
   atom has grade

   ```text
   lambda_j+r-1,
   ```

   while every THM-2293 shell-graph vertex has grade `b`. The grades align
   exactly when

   ```text
   r=b-lambda_j+1:

   r=b for the depth-one owner j=1,
   r=1 for the depth-b owner j=2.                    (41b)
   ```

   Even on this exact collision-depth branch, the outside unit `n` in
   (29) need not have THM-2293's selected root character `kappa`. Thus the
   residual interface is now the explicit pair

   ```text
   collision-depth alignment + root-color alignment. (41c)
   ```

5. **The collision depth need not be one.** For any prescribed `R`, take
   two equal-width two-interval indicators on the normalized circle,

   ```text
   B_R: centers plus or minus 1/16,
   A_R: centers plus or minus 7/16,
   epsilon_R=1/(64*13^R).                            (42)
   ```

   Then `B_R` lies in `D_1`, `A_R` lies in `D_1^c`, and their images under
   `T^s` remain disjoint for every `0<=s<=R`. The centers stay at least
   `1/8` apart, while the two image half-widths sum to at most `1/32`.
   Both functions have four jumps. Thus no bound for `r` can depend only on
   the support sides or jump counts; the mass/variation ratios in
   (32) and (37) are essential.

6. **The multiplier-four witness survives.** THM-2299's exact local
   `c_1 -> c_2` carrier has

   ```text
   (1_E)_hat(4c_1)=0.
   ```

   It also gives a sharp nontrivial collision-depth control. In that
   witness, with `epsilon=10^(-12)`,

   ```text
   B=P_13 1_E
     =(1/13)1_F,
     F centers {1,15}/16, half-width epsilon;

   f=P^2 1_E
     =(1/169)1_R,
     R centers {3,13}/16, half-width 13 epsilon;

   A=P_13 f,
     support centers {7,9}/16, half-width 169 epsilon.        (43)
   ```

   Under each further root push, the cross-gap between the two center sets
   is `3/8` at even depths and `1/8` at odd depths. Through depth eight the
   expanded half-width sum is smaller than the corresponding cross-gap.
   At depth nine the image of either `A` component already covers the
   circle. Therefore the first positive collision is exactly

   ```text
   r=9.                                                (44)
   ```

   The shell in (29) consequently has the exact additional valuation
   `r-1=8`. This explicit local current realizes an
   owner-normalized first shell eight levels above the owner and still has
   the vanished base pair atom. The present theorem may choose another
   owner multiple or that deeper shell; it does not repair the base-phase
   cancellation. The witness is local and lies far below the quantitative
   return floor.

The first-collision shell is therefore a precise new sidecar, not a
relabeling of the old residue-energy statement and not a proof of LRC(14).

## 7. Exact verification

The companion uses only integer and `Fraction` arithmetic. It checks the
word-mass floors, owner-normalization support law on exact rational grids,
the jump-product bounds, the exact valuation-shell partition, the
residue-Prony bank, and the arbitrarily delayed two-interval collision
control in (42). Reproduce with

```bash
python3 04-computation/lrc14_owner_normalized_first_collision_thm2306.py
python3 -O 04-computation/lrc14_owner_normalized_first_collision_thm2306.py
```

Every load-bearing check raises explicitly, so optimized mode runs the same
audit. QED.
