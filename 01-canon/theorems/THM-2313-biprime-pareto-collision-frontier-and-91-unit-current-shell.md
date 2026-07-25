---
id: THM-2313
title: "Biprime Pareto collision frontier and the 91-unit current shell"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For two
  nonzero nonnegative real step functions with disjoint supports, push both
  independently by the normalized Perron operators for thirteen and seven.
  Their positive-overlap cells form an upper set with a finite Pareto
  frontier. A frontier point on an axis isolates a negative exact
  one-prime shell. Every interior frontier point is the mixed Möbius corner
  and isolates a positive aggregate over frequencies coprime to 91.
  Endpoint-Prony lands a common atom with outside unit at most
  91 J_A J_B-1. For THM-2306's owner-normalized word current this is at most
  728 S^2-1 and records both exact valuations. The multiplier-four hostile
  carrier has an exact nine-point frontier, seven interior corners, and one
  grade-three corner, but it still does not select the prescribed root
  character, an incident pair edge, or cycle phase holonomy. No scalar row
  is excluded and LRC(14) remains open.
source: codex-2026-07-25-biprime-pareto-collision
depends_on:
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2306-owner-normalized-disjoint-support-and-first-collision-shell
related:
  - THM-2293-quadratic-root-energy-raises-the-ancestry-shell
  - THM-2302-same-label-expiration-dichotomy-and-pure-terminal-shell-no-go
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2310-quantitative-beta-shift-gluing-of-positive-handoffs
  - THM-2315-marked-target-gain-corolla-and-pairwise-composition-boundary
script: 04-computation/lrc14_biprime_pareto_collision_thm2313.py
output: 05-knowledge/results/lrc14_biprime_pareto_collision_thm2313.out
script_sha256: 97c16dcc67695151654108be224fbee7337afba84d72bcc3742d4e9be64d9fee
output_sha256: 1122f4df26c86d8bf3d9714dfedfb5c0b0ed5e6d04caa561eefb83d17ae39b5f
hash_basis: working-tree bytes (LF)
---

# THM-2313 -- biprime first collision is a Pareto antichain

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2306 pushed an owner-normalized source/current pair by thirteen until
its first positive spatial collision. That one-dimensional stopping time
isolated an exact thirteen-adic shell, but it left every other prime inside
an uncontrolled outside multiplier. Seven is not an arbitrary extra
direction: the desired shallow-pair edge color in THM-2293 and THM-2302
requires an outside multiplier coprime to `91=7*13`.

The correct two-prime object is not a single first time. It is the
**antichain of first collisions** in the product order on `N^2`:

```text
disjoint source/current pair
  -> upper set of positive (13-depth,7-depth) collisions
  -> finite Pareto frontier
  -> axes give one-prime shells
  -> interior corners give the square-free Möbius selector
       1-[13|n]-[7|n]+[91|n]=[gcd(n,91)=1].          (1)
```

Thus an interior corner records two exact valuation coordinates at once.
This is a genuine strengthening of the first-collision coordinate. It is
still a common source/current coefficient, not an incident edge of the
quadratic shell graph.

## 1. The abstract collision upper set

On `T=R/Z`, use

```text
P_a h(y)=(1/a)sum_(r=0)^(a-1) h((y+r)/a),           (2)

(P_a h)_hat(n)=h_hat(an),
P_a P_b=P_(ab)=P_b P_a.                             (3)
```

Let `A,B` be nonzero nonnegative real step functions satisfying

```text
AB=0 almost everywhere.                             (4)
```

Write `J_A,J_B` for their numbers of nonzero jumps. For `s,t>=0`, define

```text
A_(s,t)=P_13^s P_7^t A,
B_(s,t)=P_13^s P_7^t B,

I_(s,t)=integral_T A_(s,t) B_(s,t),                 (5)

Omega={(s,t) in N^2:I_(s,t)>0}.                     (6)
```

Here `s` is a thirteen-dynamical quotient depth. The coordinate `t` is an
auxiliary commuting seven-quotient scale introduced to resolve the
existing `91`-unit edge color; it is not elapsed time for the LRC
`T_13` orbit.

Then `I_(0,0)=0`, and `Omega` is an upper set in the coordinatewise order.
Indeed, suppose two nonnegative functions `U,V` have positive overlap and
put

```text
H={x:U(x)>0 and V(x)>0}.
```

The set `H` has positive measure. Its image under `x -> ax mod 1` has
positive measure, and every point of that image has a common preimage on
which both functions are positive. Hence `P_aU` and `P_aV` have positive
overlap. Apply this with `a=13` and `a=7`.

Both coordinate axes eventually enter `Omega`. For `p in {7,13}` and
`h in L^2(T)`,

```text
||P_p^r h-integral h||_2^2
 =sum_(n!=0)|h_hat(p^r n)|^2
 ->0.                                               (7)
```

Consequently

```text
I_(s,0)->(integral A)(integral B)>0,
I_(0,t)->(integral A)(integral B)>0.                 (8)
```

If `(s_0,0)` and `(0,t_0)` are positive, every Pareto-minimal point of
`Omega` lies in the finite rectangle

```text
[0,s_0] x [0,t_0].                                  (9)
```

The frontier is therefore a finite nonempty antichain. Every one of its
points is of exactly one geometric type:

```text
13-axis: (s,0), s>0;
7-axis:  (0,t), t>0;
interior: (s,t), s,t>0.                             (10)
```

This classifies each frontier point; it is not a claim that an entire
frontier has only one type. One frontier can contain both axis endpoints
and many interior points.

## 2. Fourier shells are the Möbius derivatives of collision

Use the Fourier convention

```text
h_hat(n)=integral_T h(x)exp(-2*pi*i*n*x)dx.          (11)
```

Parseval and (3) give an absolutely convergent bilinear series

```text
I_(s,t)
 =sum_(n in Z)
    A_hat(13^s 7^t n)
    conjugate(B_hat(13^s 7^t n)).                   (12)
```

Absolute convergence follows from Cauchy--Schwarz. Minimality now has three
different Fourier consequences.

### 2.1 A thirteen-axis point

If `(s,0)` is Pareto-minimal, then `I_(s-1,0)=0` and `I_(s,0)>0`.
Removing the sublattice of multiples of thirteen gives

```text
sum_(13 does not divide n)
 A_hat(13^(s-1)n)conjugate(B_hat(13^(s-1)n))

 =I_(s-1,0)-I_(s,0)
 =-I_(s,0)<0.                                       (13)
```

### 2.2 A seven-axis point

Likewise, a Pareto-minimal `(0,t)` gives

```text
sum_(7 does not divide n)
 A_hat(7^(t-1)n)conjugate(B_hat(7^(t-1)n))

 =I_(0,t-1)-I_(0,t)
 =-I_(0,t)<0.                                       (14)
```

### 2.3 An interior point

If `(s,t)` is interior and Pareto-minimal, its west, south, and southwest
predecessors all vanish. Inclusion--exclusion on the four divisor
sublattices gives

```text
sum_(gcd(n,91)=1)
 A_hat(13^(s-1)7^(t-1)n)
 conjugate(B_hat(13^(s-1)7^(t-1)n))

 =I_(s-1,t-1)-I_(s,t-1)-I_(s-1,t)+I_(s,t)
 =I_(s,t)>0.                                        (15)
```

Equation (15) is the mixed discrete derivative of the collision function.
Its selector is exactly

```text
1-[13|n]-[7|n]+[91|n]=[gcd(n,91)=1].                (16)
```

The sign change is structural: a first collision on one axis has the
negative first derivative `-I`, while a full two-prime corner has the
positive second derivative `+I`. Equations (13)--(15) are aggregate
statements. Earlier shells may contain individual products which cancel,
and a positive real aggregate need not contain an individually
positive-real summand. It does force at least one nonzero common Fourier
product.

This is the square-free divisor-lattice view of the object. More generally,
for distinct primes `p_1,...,p_d`, a fully interior Pareto-minimal collision
has full corner derivative

```text
sum_(gcd(n,product_i p_i)=1)
 A_hat((product_i p_i^(s_i-1))n)
 conjugate(B_hat((product_i p_i^(s_i-1))n))

 =(-1)^d I_s.                                       (17)
```

Indeed the divisor-lattice Möbius expansion has `2^d-1` zero lower
corners and one positive top corner. The present theorem uses only `d=2`;
(17) explains why adding a second prime is a change of representation, not
just a longer one-prime search.

## 3. Residuewise endpoint-Prony landing

The nonzero aggregate can be landed at bounded frequency without truncating
Fourier tails. At the relevant predecessor put

```text
U=P_13^(s-1)P_7^(t-1)A,
V=P_13^(s-1)P_7^(t-1)B                             (18)
```

for an interior corner, with the absent coordinate omitted on an axis.
Perron transport cannot increase the number of jump locations, so `U,V`
have at most `J_A,J_B` jumps.

For a positive integer `q`, integration by parts makes

```text
C_q=(2*pi*q)^2 U_hat(q)conjugate(V_hat(q))          (19)
```

an exponential sum on at most

```text
L<=J_AJ_B                                            (20)
```

distinct endpoint-difference nodes. Restrict to one residue progression

```text
q=a+m ell.
```

After equal nodes are combined, this is an `L`-node exponential sum in
`ell`. A nonzero such sequence cannot vanish at `L` consecutive integers,
by the Vandermonde determinant.

At least one positive residue progression occurring in (13), (14), or
(15) is not the zero sequence. Testing `ell=0,...,L-1` therefore gives the
following sharp bookkeeping bounds:

```text
13-axis:
  13 does not divide n,
  1<=n<=13J_AJ_B-1;

7-axis:
  7 does not divide n,
  1<=n<=7J_AJ_B-1;

interior:
  gcd(n,91)=1,
  1<=n<=91J_AJ_B-1.                                (21)
```

The constants come from the largest allowed first residue:
`12+13(L-1)=13L-1`, `6+7(L-1)=7L-1`, and
`90+91(L-1)=91L-1`.

Now instantiate THM-2306. Its positive word-restricted current and source
are

```text
f=1_Q P_13^k 1_E,
A=P_c f,
B=P_c 1_E,                                         (22)

support(A) subset D_1^c,
support(B) subset D_1,
J_A<=4S,
J_B<=2S.                                           (23)
```

Write the owner in both prime coordinates:

```text
c=13^lambda 7^gamma w,             gcd(w,91)=1.    (24)
```

Since `A_hat(q)=f_hat(cq)` and
`B_hat(q)=(1_E)_hat(cq)`, every frontier point supplies a positive
frequency `N` at which

```text
f_hat(N)!=0,
(1_E)_hat(N)!=0.                                   (25)
```

The three cases retain exactly the following coordinates:

```text
13-axis (s,0):
 N=c 13^(s-1)n,
 13 does not divide n,
 n<=104S^2-1,
 nu_13(N)=lambda+s-1;
 nu_7(n) is uncontrolled.

7-axis (0,t):
 N=c 7^(t-1)n,
 7 does not divide n,
 n<=56S^2-1,
 nu_7(N)=gamma+t-1;
 nu_13(n) is uncontrolled.

interior (s,t):
 N=c 13^(s-1)7^(t-1)n,
 gcd(n,91)=1,
 n<=728S^2-1,
 (nu_13(N),nu_7(N))
   =(lambda+s-1,gamma+t-1).                         (26)
```

Thus only an interior frontier point gives a genuinely `91`-unit residual
core after both displayed prime powers are removed. The bound in (26) is
on that residual core; the two collision coordinates retain all
ramification separately.

## 4. Exact nine-corner atlas of the multiplier-four carrier

The local hostile row of THM-2299 has profile `(1,3,5)`, source owner
`c_1=13`, and `epsilon=10^(-12)`. THM-2306 identifies its normalized
supports exactly:

```text
B=P_13 1_E
  =(1/13)1_F,
  F centers {1,15}/16, half-width epsilon;

f=P_13^2 1_E
  =(1/169)1_R,
  R centers {3,13}/16, half-width 13 epsilon;

A=P_13 f
  =(1/2197)1_C,
  C centers {7,9}/16, half-width 169 epsilon.       (27)
```

After `s` further thirteen-pushes and `t` seven-pushes, put

```text
M=13^s 7^t.
```

The two center sets are multiplied by `M` modulo sixteen. Their cross-gap
is independent of `t` and is

```text
g_s=3/8 if s is even,
g_s=1/8 if s is odd.                                (28)
```

The pushed supports have positive overlap exactly when

```text
170 M epsilon>g_s.                                  (29)
```

If neither image interval already covers the circle, this is precisely the
strict sum-of-radii test. If one does cover, (29) is automatic. With
`epsilon=10^(-12)`, (29) becomes the integer test

```text
s even: 1360*13^s*7^t>3*10^12;
s odd:  1360*13^s*7^t>  10^12.                     (30)
```

For `s=0,...,9`, the least positive `t` is

```text
(12,10,9,7,6,4,4,2,1,0).                           (31)
```

The cell `(6,4)` is dominated by `(5,4)`. The complete Pareto frontier,
with its exact scale `M`, is therefore:

| corner `(s,t)` | type | `13^s 7^t` | exact mixed valuation for `c=13` |
|---:|:---|---:|:---|
| `(0,12)` | 7-axis | `13841287201` | only `nu_7=11` is fixed |
| `(1,10)` | interior | `3672178237` | `(nu_13,nu_7)=(1,9)` |
| `(2,9)` | interior | `6819759583` | `(2,8)` |
| `(3,7)` | interior | `1809323971` | `(3,6)` |
| `(4,6)` | interior | `3360173089` | `(4,5)` |
| `(5,4)` | interior | `891474493` | `(5,3)` |
| `(7,2)` | interior | `3074677333` | `(7,1)` |
| `(8,1)` | interior | `5710115047` | `(8,0)` |
| `(9,0)` | 13-axis | `10604499373` | only `nu_13=9` is fixed |

There are nine minimal cells, seven of them interior. This is substantially
richer than the one-prime first depth `r=9`: the two-prime object sees a
whole antichain of incomparable first coalescences.

Both hostile normalized supports have four jumps. Its local mixed
endpoint-Prony bound is consequently

```text
91*4*4-1=1455,                                      (32)
```

rather than the global `728S^2-1`.

## 5. The grade-three corner still misses color and incidence

At THM-2302's common shell clock `b+1`, an owner with
`nu_13(c)=lambda` reaches shell grade `b` at a mixed corner exactly when

```text
lambda+s-1=b,
s=b-lambda+1.                                      (33)
```

For the hostile carrier, `lambda=1` and `b=3`. Exactly one point of its
frontier has the required thirteen-coordinate:

```text
(s,t)=(3,7).                                        (34)
```

The new theorem therefore does find a grade-three current/source atom in
this sharp local obstruction. But after division by `13^3`, its root
character is

```text
N/13^3
 congruent 7^6 n
 congruent -n                  (mod 13).             (35)
```

The hostile row's shallow-pair/root character is `kappa=4`. Matching it
would require

```text
n congruent 9 (mod 13),                             (36)
```

and the mixed aggregate does not select that residue. Moreover, because
`t=7`, the complete multiplier outside `13^3` contains `7^6`; it is not
coprime to `91` even though the **residual core after removing the
seven-coordinate** is.

Two further losses remain even if (36) happens accidentally:

- `f_hat(N)` is the coefficient of the word-restricted product
  `1_QP^k1_E`, not a second bare source vertex;
- a shared source/current frequency is not an incident
  `A-A'=m c_3` edge of THM-2293's shell graph.

There is also no hidden target-gain landing. THM-2315 identifies a fork
word with a twelve-point projective gain fibre. The present aggregate proves
that at least one unit residue progression fires, but does not canonically
choose that progression before endpoint-Prony and does not map its residue
to a point of the gain fibre. Conversely, an exact target gain does not
select this theorem's Fourier residue. These are independent address
coordinates until a new incidence theorem couples them.

The interior corner `(8,1)` avoids an added factor of seven, but it has
thirteen-grade eight rather than the hostile profile's grade three. The
frontier exposes the tradeoff exactly:

```text
grade alignment, root-character alignment, and 91-unit edge color

are three separate coordinates.                    (37)
```

## 6. Symmetric coalescence is not forward gluing

There is a useful but limited junction interpretation. For consecutive pure
transition carriers

```text
(F,k):X -> Y,
(G,l):Y -> Z,
```

let

```text
D=P_13^k1_F,             B=1_G.                     (38)
```

If `integral_G D>0`, a positive subset of `F` arrives directly in `G`, so
the two carriers compose with zero wait. If `DB=0`, the present theorem
applies to `D,B` and supplies a Pareto family of **symmetric future
coalescences** and their signed Fourier shells.

THM-2310 supplies a different statement: after a finite one-sided
thirteen-wait,

```text
Arr(F,k) intersection T_13^(-r)G
```

has positive measure. Its wait `r` is not either coordinate of the
symmetric Pareto frontier. The former compares a current set with the
future preimage of the next carrier; the latter pushes both functions into
a common future quotient.

Around a pure owner cycle one may therefore attach a finite forward wait
and a finite signed collision antichain at every junction. Different
junctions can choose different corners, frequencies, residues, and complex
phases. There is no telescoping phase product and no cycle holonomy.
THM-2305's fork remains a genuine hyperedge and is untouched.

## 7. Connection and loss ledger

The new representation can be summarized without conflating its outputs:

```text
source:
  THM-2306's disjoint owner-normalized source and exact word current;

object:
  the collision upper set Omega in the divisibility lattice N^2;

boundary invariant:
  its finite Pareto antichain, not one arbitrarily chosen first time;

map:
  take the square-free Möbius derivative across a minimal corner;

preserved:
  owner c, complete word-restricted current f, original source E,
  common frequency, exact 13- and 7-valuations at an interior corner,
  and an explicit residual-unit landing bound;

destroyed or unselected:
  prescribed root character, base component phase, pair-edge incidence,
  target-plane gain, one-sided orbit order, a common corner across
  different junctions, and cycle holonomy;

cheapest hostile probe:
  the exact multiplier-four two-interval carrier, whose nine-point
  frontier realizes both the new gain and every remaining loss.        (39)
```

The square-free Möbius derivative is reusable: whenever several commuting
Perron directions are genuinely intrinsic, a collision antichain can
separate their valuation coordinates. It should not be added cosmetically.
Here seven is justified by the existing `gcd(m,91)=1` edge color.

No scalar counterexample profile is eliminated, the local hostile carrier
is not promoted to a global scalar cover, and LRC(14) remains open.

## 8. Exact verification

The companion uses only integer and `Fraction` arithmetic. It exhausts the
upper-set/frontier rules on every `4 x 4` Boolean grid, verifies all three
divisibility selectors, all mod-13, mod-7, and mod-91 Prony banks through
128 endpoint nodes, the global jump bounds, and 123750 exact two-prime
valuation rows. It also computes the hostile cross-gap law, full collision
staircase, nine-point Pareto frontier, seven valuation pairs, local
`1455` bound, and the grade/color stopping residue.

Reproduce with

```bash
python3 04-computation/lrc14_biprime_pareto_collision_thm2313.py
python3 -O 04-computation/lrc14_biprime_pareto_collision_thm2313.py
```

Every load-bearing check raises explicitly, so optimized mode executes the
same audit. QED.
