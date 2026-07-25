---
id: THM-2299
title: "Rooted current-service energy and the base-phase no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. On the
  quantitative prescribed-return arm of THM-2296, the same selected owner,
  exact expiration clock, and current blocker-only residual support every
  nonzero C_13 root character. Moreover one named target blocker receives a
  directed edge of at least half the return mass; every character has energy
  strictly larger than one quarter of the squared return floor on that same
  named edge. A covariant endpoint-Prony argument gives, in each character,
  a signed named-service coefficient with gauge index at most 4S-1. This
  does not select the THM-2276 pair atom. An exact local profile-(1,3,5),
  multiplier-four witness retains a c_1-to-c_2 handoff at prescribed
  expiration, a one-sheet nonzero word in every root character, and a
  c_1-anchored minor equal to 8 modulo 13, while its exact pair atom
  vanishes. The witness is not a global scalar cover, no profile is
  excluded, and LRC(14) remains open.
source: codex-2026-07-25-rooted-service-base-phase
depends_on:
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
  - THM-2296-prescribed-expiration-return-or-bounded-ancestry-resonance
related:
  - THM-2276-shallow-owner-residue-aligned-crossing
  - THM-2278-two-shallow-proper-root-spectrum-and-gap-ancestry-activation
  - THM-2282-thirteen-adic-saturation-and-unit-anchored-minor
script: 04-computation/lrc14_rooted_current_service_phase_no_go_thm2299.py
output: 05-knowledge/results/lrc14_rooted_current_service_phase_no_go_thm2299.out
script_sha256: d25bf84532f5b673bf98596000cdbe2c9e18b13f212141018e3aec963eb942f7
output_sha256: a228c6c1a02e8ba91a85e9476b4d44529666175e98527a283909e0ef3e01b92d
hash_basis: working-tree bytes (LF)
---

# THM-2299 -- current service fires every root character, but not the base phase

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2296 retains the owner label and prescribed expiration clock, but its
scalar Perron density is the trivial root-character trace. The lawful
refinement is not to discard that density. Keep its thirteen predecessor
values as one rooted vector over the current blocker-only target.

On the positive-return arm, every nonzero character of that vector fires
pointwise and quantitatively. A gauge correction turns it into a periodic
current-service function and endpoint-Prony gives a bounded signed
coefficient in the chosen character. The smallest hard multiplier then
identifies the exact remaining loss:

```text
owner + exact clock + current target service + root character + anchored minor

still does not determine

the signed base-phase coefficient at the prescribed pair frequency.      (1)
```

## 1. The prescribed-return input

Use THM-2296's selected label `j`, source, current target, and clock:

```text
E_j
 =A_0 intersection D_(c_j)
       minus union_(ell!=j)D_(c_ell),

R_j=A_0 minus D_(c_j),

c_j=13^(lambda_j)u_j,
k_j=lambda_j+1.                                      (2)
```

Put

```text
G_j=P^(lambda_j)1_(E_j),                             (3)

P f(y)=(1/13)sum_(r=0)^12 f((y+r)/13).
```

Owner transport gives

```text
G_j>=0,
G_j takes values in 13^(-lambda_j) Z,
support(G_j) subset D_(u_j).                         (4)
```

The prescribed return is

```text
rho_j
 =measure(E_j intersection T^(-k_j)R_j)
 =integral_(R_j) P G_j(y)dy.                         (5)
```

On the quantitative return arm of THM-2296,

```text
strict:
 rho_j>=delta_s:=39002430583/53493927587100;

repeated-first:
 rho_j>=delta_r:=13560199813/53493927587100.         (6)
```

The target in (5) is still the literal current-service set `R_j`. Almost
every point of it lies in the guard-safe, five-unit-mask-free residual,
outside `D_(c_j)`. The global scalar cover therefore forces one of the other
two blockers at that same point.

Let those two labels be `a,b`, in a fixed order, and partition current
service by

```text
R_(j,a)=R_j intersection D_(c_a),

R_(j,b)=R_j intersection D_(c_b) minus D_(c_a).     (6a)
```

Up to the null exceptional set of the cover, these sets are disjoint and
cover `R_j`. Their ancestry return masses

```text
rho_(j,t)=integral_(R_(j,t))P G_j,       t in {a,b},
```

sum to `rho_j`. Choose a named label `t` with

```text
rho_(j,t)>=rho_j/2.                                  (6b)
```

This is an intrinsic time-oriented colored edge

```text
c_j -> c_t                                             (6c)
```

at the prescribed clock.

## 2. A rooted current-service word

Write

```text
zeta=exp(2*pi*i/13),
x_r(y)=(y+r)/13,                         0<=r<=12,

v_r(y)=G_j(x_r(y)),

M_a(y)=sum_(r=0)^12 v_r(y)zeta^(-ar),
                                      1<=a<=12.      (7)
```

If `v(y)` is nonzero, its support has at most two entries. Indeed, by (4)
every occupied root lies in `D_(u_j)`. Since `u_j` is a thirteen-unit, it
permutes the thirteen root labels, while an arc of length `1/7` meets at most
two of thirteen equally spaced roots.

For one positive entry, `M_a(y)` plainly does not vanish. For two positive
entries `A,B` at distinct roots, put

```text
omega=zeta^(a(r-s))!=1.                              (8)
```

The angle of `omega` is one of the twelve nonzero thirteenth-root angles.
The closest such angle to `pi` is at distance `pi/13`. Hence

```text
|A+B omega|
 >=(A+B)sin(pi/26)
 >(A+B)/13.                                         (9)
```

The first inequality is minimized when `A=B`; the second is the elementary
concavity bound

```text
sin(pi/26)>2(pi/26)/pi=1/13.                         (10)
```

Since

```text
sum_r v_r(y)=13 P G_j(y),                            (11)
```

equations (9)--(11) give, pointwise for every `a=1,...,12`,

```text
|M_a(y)|^2>(P G_j(y))^2
```

whenever `P G_j(y)>0`.

Integrate first on the full current target. Cauchy--Schwarz and
`measure(R_j)<=1` give

```text
integral_(R_j)|M_a(y)|^2dy
 >integral_(R_j)(P G_j(y))^2dy
 >=rho_j^2/measure(R_j)
 >=rho_j^2.                                         (12)
```

Thus all twelve rooted characters fire on current blocker service, with the
explicit branch floors

```text
strict:
 integral_(R_j)|M_a|^2
 >(39002430583/53493927587100)^2;

repeated-first:
 integral_(R_j)|M_a|^2
 >(13560199813/53493927587100)^2.                   (13)
```

This is stronger than unlocalized root energy: the factor `1_(R_j)` has not
been pushed, averaged, or replaced by an eventual target.

More sharply, integrate on the named set chosen in (6b):

```text
integral_(R_(j,t))|M_a(y)|^2dy
 >rho_(j,t)^2
 >=rho_j^2/4.                                       (13a)
```

Thus every character fires on the same named directed blocker edge. The
uniform named-edge floors are

```text
strict:
 (39002430583/53493927587100)^2/4;

repeated-first:
 (13560199813/53493927587100)^2/4.                  (13b)
```

## 3. The covariant gauge and a bounded signed lift

The raw root word is not periodic. Reindexing its roots gives

```text
M_a(y+1)=zeta^a M_a(y).                              (14)
```

The gauge-corrected function

```text
N_a(y)=exp(-2*pi*i*a*y/13)M_a(y)                    (15)
```

is periodic. Without a target restriction it obeys the exact affine
Fourier law

```text
(N_a)_hat(h)=13 (G_j)_hat(a+13h),          h in Z.  (16)
```

Retain the named current service by defining

```text
W_a=1_(R_(j,t))N_a.                                 (17)
```

Equation (13a) says `W_a` is nonzero for every nonzero character.

Let `J_G` be the number of jumps of `G_j` and let `K_t` be the number of
jumps of `1_(R_(j,t))`. As in THM-2296, Perron transport cannot create jump
locations. The named target is a Boolean combination of a subset of the same
nine scalar boundary banks used by `E_j`. Hence

```text
J_G<=2S,
K_t<=2S,                                            (18)

S=H+sum_i q_i+sum_ell c_ell.
```

Each jump of a root section `G_j((y+r)/13)` is the image under `T` of one
jump of `G_j`; summing the sections creates no new locations. Therefore the
step amplitude

```text
1_(R_(j,t))M_a
```

has at most

```text
B_a<=J_G+K_t<=4S                                    (19)
```

circle jumps.

On every constancy cell of that amplitude, `W_a` is a constant times
`exp(-2*pi*i*a*y/13)`. Distributional differentiation gives

```text
2*pi*i(h+a/13)(W_a)_hat(h)
 =sum_(x in Jump(W_a)) Delta_a(x)exp(-2*pi*i*h*x).  (20)
```

The right side is a nonzero exponential sum with at most `B_a` distinct
nodes. It is nonzero because otherwise the periodic function `W_a` would
solve

```text
W_a'=-(2*pi*i*a/13)W_a
```

globally, whose only periodic solution for `1<=a<=12` is zero. The
consecutive-zero Vandermonde lemma now supplies

```text
0<=h<=B_a-1<=4S-1
```

such that

```text
(W_a)_hat(h)
 =integral_(R_(j,t)) M_a(y)
    exp(-2*pi*i(h+a/13)y)dy
 !=0.                                                (21)
```

Equation (21) is a bounded, gauge-faithful, signed current-service lift in
the exact residue `a`. It is not a claim that
`(G_j)_hat(a+13h)` and `1_(R_j)_hat(h)` are separately nonzero. Restriction
to `R_j` is a convolution in Fourier space.

## 4. The smallest hard multiplier keeps every label and kills the pair atom

The exact stopping witness already occurs at multiplier four and the
smallest interior profile `(1,3,5)`. Take

```text
H=1,

(q_1,q_2,q_3,q_4,q_5)=(4,2,3,6,10),

(c_1,c_2,c_3)=(13,13^3,2*13^5).                    (22)
```

All six guard/unit coefficients are thirteen-units, `H` is odd, the five
`q_i` are distinct, and the blockers have the required strict valuations.
There is the primitive shallow pair

```text
p: 13q_1-4c_1=0,

K=A(p)=52,
kappa=K/13=4.                                       (23)
```

There is also an independent relation

```text
r:q_1-2q_2=0.                                       (24)
```

On the columns `(c_1,q_2)`, the rows `(p,r)` have minor

```text
det [[-4,0],[0,-2]]=8!=0 mod 13.                    (25)
```

Thus even the anchored mod-thirteen transversality sidecar is present.

Put

```text
epsilon=10^(-12),

F=
 (-1/16-epsilon,-1/16+epsilon)
 union
 ( 1/16-epsilon, 1/16+epsilon),                    (26)

E={ (8+y)/13 : y in F }.                            (27)
```

The root `8` is pair-aligned because

```text
4*8=6 mod 13.                                       (28)
```

At the two source centers,

```text
x=127/208, 129/208,
```

the guard is active, all five unit masks and both deep blockers are safe,
and `c_1` is strictly dangerous. At the prescribed terminal time,

```text
T^2x=3/16, 13/16,
```

the guard is active, all five unit masks, `c_1`, and `c_3` are safe, while
`c_2` is strictly dangerous.

The norm margin divided by the relevant slope in the interval coordinate is
at least

```text
3/540602608>epsilon.                                (29)
```

The minimum occurs for terminal `c_3`; hence the complete open intervals in
(26)--(27), not only their centers, retain every strict label.

Multiplication by thirteen maps `E` one-to-one onto `F` and maps `F`
one-to-one on each component onto the current `c_2`-only target

```text
R=T(F).                                              (30)
```

Consequently every point of `E` makes the prescribed

```text
c_1 -> c_2
```

handoff at time two. Moreover,

```text
G=P1_E=(1/13)1_F.                                   (31)
```

On every terminal root fibre over `R`, the ancestry word has exactly one
nonzero entry, of weight `1/13`. Thus every nonzero rooted character has
pointwise squared magnitude `1/169`, with the current target label retained.

Nevertheless the equal intervals in (26) have opposite quarter-turn phases
at frequency four:

```text
(1_F)_hat(4)=0.                                     (32)
```

Equations (23), (27), and (31) give

```text
(1_E)_hat(52)=0,
G_hat(4)=0.                                         (33)
```

For the pair character `a=4`, the support of the rooted word is exactly
`R`, so (16)--(17) make the same cancellation

```text
(W_4)_hat(0)=13G_hat(4)=0.                          (34)
```

Equations (25), (30), (31), and (34) are the sharp no-go: even one-sheet
root ancestry, exact current `c_2` service, the correct clock, and an
anchored unit minor do not prevent cancellation of the signed base phase
between two terminal components.

## 5. Connection and loss ledger

```text
source:
  THM-2296's quantitative prescribed-return arm;

map:
  stop one step before expiration, retain all thirteen Perron predecessor
  values, Fourier-transform in the root label, and apply the affine gauge;

preserved:
  selected owner j, exact depth lambda_j, exact clock lambda_j+1,
  ancestry multiplicity, current residual R_j, one named directed target
  blocker c_t, and every nonzero C_13 root character;

quantitative output:
  uncolored rooted energy >rho_j^2, named-edge rooted energy >rho_j^2/4,
  and a named-service signed gauge coefficient with 0<=h<=4S-1;

destroyed by taking only energy:
  the complex phase between different terminal base cells;

smallest hostile test:
  multiplier m=4, profile (1,3,5), two quarter-turn-separated cells,
  one rooted sheet, prescribed c_1-to-c_2 service, anchored minor 8;

needed sidecar:
  a genuinely global phase-coherence or signed cover identity coupling
  distinct terminal cells at the exact pair base index.                 (35)
```

The witness satisfies every displayed local scalar label and relation. It
is **not** asserted to satisfy the global cover away from its two source and
two target intervals. It therefore does not refute LRC, exclude a profile,
or decide which THM-2296 arm a genuine global counterexample occupies.

## 6. Exact verification

The companion uses only exact integer and `Fraction` arithmetic. It checks
the two return floors and squared energy floors, all twelve root characters,
the covariant jump ledger, the complete scalar row and valuation profile,
both relation rows and the anchored minor, every source and target label,
the exact minimum half-width margin, the one-sheet transport and masses, and
the quarter-turn cancellations (32)--(34). Reproduce with

```bash
python3 04-computation/lrc14_rooted_current_service_phase_no_go_thm2299.py
python3 -O 04-computation/lrc14_rooted_current_service_phase_no_go_thm2299.py
```

Every load-bearing check raises explicitly, so optimized mode executes the
same audit. QED.
