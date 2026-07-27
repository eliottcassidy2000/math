---
id: THM-2462
title: "Mixed-radix root-phase orbit universality"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. The
  nearest-integer phase law converts every fixed
  oriented C_13 danger mask into one rational open phase interval.
  Any finite list of such open root-chart words can be realized on
  disjoint positive intervals along one fixed one-dimensional orbit
  by a mixed-radix divisibility chain with prescribed nonzero speed
  residues. Explicit speeds
  (14,644,27048,1190112,57125376,2856268800) realize all thirteen
  clean guard-danger charts of the THM-2458 atlas with the same
  literal gate {0,1}; the thirteen parent intervals are disjoint and
  have equal exact width 1/6664627200. A separate primitive affine
  speed tuple gives the same thirteen words. Thus unrestricted
  fixed-speed phase compatibility cannot exclude the static hostile.
  This is only a six-role local root-atlas realization, not a full
  fourteen-speed scalar cover, owner/valuation/word continuation, or
  canonical phase current. No scalar row is excluded.
source: codex-2026-07-26-mixed-radix-root-phase-orbit
depends_on: []
related:
  - THM-2401-common-filter-endpoint-or-first-death-certificate
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2456-two-root-replica-uniform-offset-boundary
  - THM-2458-clean-root-guard-danger-thirteen-chart-uniform-offset-hostile
script: 04-computation/lrc14_mixed_radix_root_phase_orbit_thm2462.py
output: 05-knowledge/results/lrc14_mixed_radix_root_phase_orbit_thm2462.out
script_sha256: e037ad210e9a152b8c145b5511b266422fa49da6627550626395db4de2e64263
output_sha256: 71ba61d6c0e5583b9ee57264c92bd6dab4986500ecf2a13f1c73a93701c3e465
hash_basis: working-tree bytes (LF)
---

# THM-2462 -- mixed-radix root-phase orbit universality

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The clean-root atlas reserved in THM-2458 is static: its thirteen
rows are individually valid finite root covers, but residue data alone
does not show that one fixed integer speed tuple visits every row.
This theorem supplies the missing phase coordinate and gives the
opposite of a phase obstruction:

```text
finite family of rational open root-chart words
  -> one fixed mixed-radix speed tuple
  -> one positive parent interval for every word.                 (1)
```

The construction is general.  In particular, one-dimensional
fixed-speed compatibility by itself cannot exclude a finite open
root atlas when the integer lifts are unrestricted.

## 1. Exact nearest-integer phase law

Use the common oriented physical section

```text
x=(y+r)/13,              r in C_13.                              (2)
```

Let `u>0` be a `13`-unit and write

```text
rho=u mod 13,

N_u(y)=floor(uy+1/2),

delta_u(y)=uy-N_u(y) in [-1/2,1/2).                              (3)
```

For every root `r`,

```text
u(y+r)/13
 =integer+(uy+rho r)/13.                                        (4)
```

Put

```text
s=-rho^(-1)N_u(y) mod 13.                                       (5)
```

Away from endpoints, the ordinary danger mask

```text
D_u(y)={r:||u(y+r)/13||<1/14}
```

obeys

```text
-1/2<delta_u(y)<-1/14
  -> D_u(y)={s,s+rho^(-1)};

 1/14<delta_u(y)< 1/2
  -> D_u(y)={s,s-rho^(-1)}.                                     (6)
```

Indeed, after multiplying the danger inequality by thirteen, the
integer offset from the nearest value is respectively `0,1` or
`0,-1`.

For the guard danger radius `1/7`, the same argument gives

```text
-1/2<delta_u(y)<-1/7
  -> G_u(y)
     ={g,g+rho^(-1),g+2rho^(-1),g+3rho^(-1)},

g=-rho^(-1)(N_u(y)+1);                                          (7)

 1/7<delta_u(y)<1/2
  -> G_u(y)
     ={g,g+rho^(-1),g+2rho^(-1),g+3rho^(-1)},

g=-rho^(-1)(N_u(y)+2).                                          (8)
```

The integer offsets in (7) are `-1,0,1,2`; those in (8) are
`-2,-1,0,1`.

Thus every fixed oriented ordinary pair or four-root guard is an open
interval condition on

```text
uy mod 13.                                                       (9)
```

This is the load-bearing representation for the construction below.

## 2. The thirteen target words

Use the residue and step signature

```text
(rho_K,rho_q,rho_H,rho_A,rho_B,rho_C)
 =(1,7,8,1,9,8),                                                (10)

(d_q,d_H,d_A,d_B,d_C)=(2,5,1,3,5).
```

Fix the literal gate

```text
K={0,1}.                                                        (11)
```

For a row `(g;q,a,b,c)`, the desired masks are

```text
G_g={g,g+5,g+10,g+15},

Q_q={q,q+2},

A_a={a,a+1},

B_b={b,b+3},

C_c={c,c+5}.                                                     (12)
```

The thirteen rows are

```text
(0; 9, 7,3,12),   (1;2, 9,5, 7),   (2;3,10,6, 3),
(3; 2,10,6, 7),   (4;10,7,2,11),   (5;6, 3,9,11),
(6; 2, 9,2, 7),   (7;8, 2,2, 6),   (8;2, 6,9,11),
(9; 8, 3,2, 7),  (10;9, 4,3, 3),  (11;2, 6,9, 5),
(12;8, 6,2,11).                                                (13)
```

Select the negative-delta branch in (6)--(7).  One real lift of the
six target phase arcs is

```text
J_K       =(-1/2,-1/14),

J_q(q)    =(-7q-1/2,   -7q-1/14),

J_H(g)    =(-1-8g-1/2, -1-8g-1/7),

J_A(a)    =(-a-1/2,    -a-1/14),

J_B(b)    =(-9b-1/2,   -9b-1/14),

J_C(c)    =(-8c-1/2,   -8c-1/14),                              (14)
```

read modulo thirteen.  Equations (6)--(8) turn membership in (14)
exactly into (11)--(12).

## 3. An explicit fixed-speed realization

Take successive multipliers

```text
(b_1,b_2,b_3,b_4,b_5)=(46,42,44,48,50)                        (15)
```

and fixed speeds

```text
u_0=14,

u_(i+1)=b_(i+1)u_i.
```

Explicitly,

```text
(u_K,u_q,u_H,u_A,u_B,u_C)
 =(14,644,27048,1190112,57125376,2856268800).                  (16)
```

Their residues are exactly (10), since the multiplier residues are

```text
(7,3,5,9,11)
 =(7/1,8/7,1/8,9/1,8/9) mod 13.                               (17)
```

For each row, start with `I_5=J_C(c)` and recursively choose an integer
`k_i` such that

```text
I_i
 =J_i intersect (I_(i+1)+13k_i)/b_(i+1)
```

is nonempty, for `i=4,3,...,0`.                                (18)

The exact choices used by the companion give

```text
|I_0|=1/476044800                                              (19)
```

for every one of the thirteen rows.  If `x_0` is the midpoint of
`I_0`, set

```text
y=(x_0+13)/14.                                                  (20)
```

Then `0<y<1`, and induction through (18) gives

```text
u_i y mod 13 in J_i                         for every i.         (21)
```

The thirteen exact midpoint values, in the order (13), are

```text
3034865843/3332313600,
9058389143/9996940800,
17899453493/19993881600,
18124187051/19993881600,
8996707229/9996940800,
501517651/555385600,
362322469/399877632,
3604728023/3998776320,
823790063/908812800,
360620951/399877632,
18217325957/19993881600,
9062455079/9996940800,
3003819443/3332313600.                                        (22)
```

Around each midpoint, the full open parent interval has width

```text
1/6664627200.                                                   (23)
```

The smallest distance between two displayed midpoints is

```text
3341/102009600,                                                 (24)
```

so the thirteen intervals are pairwise disjoint by a very large
margin.  On every point of each interval, not merely at its midpoint,
the six physical masks are exactly the corresponding row of
(11)--(13).

Consequently the normalized restriction of Haar measure to the union
of the thirteen intervals assigns equal chart weight `1/13`.  Each
root lies in four of the guard supports, so the resulting averaged
guard density is the positive constant

```text
4/13.                                                           (25)
```

Rescaling the common measure by `13/4` gives the canonical static
atlas weights `1/4`.  Positivity and uniformity, rather than this
inessential common scale, are the uniform-offset content.

The six speeds in (16) share a common factor only because it makes the
certificate especially transparent.  The companion also verifies the
primitive affine tuple

```text
(14,657,27607,1214721,58306621,2915331063),                     (26)
```

which has gcd one and the same residues.  Exact interval intersection
again realizes all thirteen words.  Thus a common dilation is not the
source of the phenomenon.

## 4. General finite open-cell realization

The construction is not special to (13).

> **Mixed-radix finite-word lemma.** Let `p` be a prime, let
>
> ```text
> rho_0,...,rho_m in F_p^*,
> ```
>
> and let `W` be a finite set of words.  For every `w in W` and every
> coordinate `i`, let
>
> ```text
> J_(w,i) subset R/pZ
> ```
>
> be a nonempty rational open interval.  Then there are fixed,
> pairwise distinct positive integers
>
> ```text
> u_i=rho_i mod p
> ```
>
> and pairwise disjoint positive rational intervals `Y_w subset T`
> such that
>
> ```text
> u_i y mod p in J_(w,i)
> ```
>
> for every `w`, every `i`, and every `y in Y_w`.                (27)

### Proof

Choose the speeds as a divisibility chain.  Work backwards.  Suppose
that for coordinate `i+1` and each word a nonempty rational interval

```text
K_(w,i+1) subset J_(w,i+1)
```

has been retained.  Choose one positive integer

```text
b_i=rho_(i+1)rho_i^(-1) mod p                                  (28)
```

so large that, simultaneously for every word, one component of

```text
(K_(w,i+1)+pk)/b_i,             k in Z,
```

lies compactly inside `J_(w,i)`.  This is possible because the
component starts have mesh `p/b_i` and their lengths tend to zero.
Quantitatively, if real lifts have lengths `L_w,L'_w`, it is enough
to take

```text
b_i>(p+L'_w)/L_w
```

for every word.  Retain such a component as `K_(w,i)`.

After reaching coordinate zero, choose

```text
u_0=rho_0 mod p
```

large enough to place distinct translates

```text
Y_w=(K_(w,0)+pn_w)/u_0
```

inside `(0,1)`.  Put

```text
u_(i+1)=b_i u_i.                                                (29)
```

For `y in Y_w`, induction gives

```text
u_i y
 =(b_(i-1)...b_0)(u_0y)
 mod p
 in K_(w,i)
 subset J_(w,i).                                                (30)
```

All choices can be rational, and the starting translates make the
`Y_w` disjoint.  Shrinking them to one common positive rational
length gives equal chart weights without changing (30). QED.

The same proof works for any finite union of open arcs by choosing
one component in each coordinate.  Primality is used here only to
write the prescribed nonzero residue ratios in (28); a unit sequence
works over any modulus.

### 4.1 Root-constant labelled bits can also be appended

The lemma also explains exactly how far finite owner/word conditioning
can go.  Append another unit speed `v` to the mixed-radix chain and
use the physical speed

```text
c=13v.                                                          (31)
```

On the common root chart,

```text
c(y+r)/13=v(y+r)=vy+integer.                                   (32)
```

Therefore the truth value of `D_c` is root-constant and is determined
only by `vy mod 1`.  A prescribed strict danger bit contains, for
example, the rational phase interval

```text
(-1/28,1/28),
```

while a strict safe bit contains

```text
(1/7,3/14).                                                     (33)
```

View either interval as one component in `R/13Z` and append it to
the target word.  Repeating this construction proves:

> **Labelled root-constant corollary.** In the finite-word lemma, one
> may additionally prescribe any finite list of independently
> supplied root-constant strict safe/danger bits on every target
> chart.  The fixed physical speeds of those factors may all be
> chosen divisible by thirteen, distinct, and positive.  The
> realizing parent intervals may still be made disjoint and
> equal-weight.

The companion checks a nontrivial two-label instance.  It appends unit
residues `1,2`, using multipliers

```text
1305,1302,
```

and prescribes, respectively, the parity word `g mod 2=0` and the
thirteen-bit word

```text
(0,0,1,1,0,1,1,1,1,0,1,1,0).                                  (34)
```

All twenty-six labelled truth values occur on exact equal-width
intervals.

This corollary treats every appended bit as an independent physical
factor.  It does not separate incompatible requirements on repeated
uses of the **same** speed.  A delayed semantic word in which one
speed is transported through several related endpoint maps must
first have a nonempty common phase cell for those linked uses.  Nor
does appending independently chosen blocker factors prove the global
scalar-cover condition.

## 5. What the realization does and does not settle

The source, target, map, and loss ledger are:

```text
source:
  a finite static atlas of open cyclic root masks with fixed residues;

map:
  replace each mask by its nearest-integer phase interval and nest
  preimages along a rapidly separated integer-speed chain;

preserved:
  one common oriented root chart, one fixed speed for each of the six
  roles, every literal mask, positive interval mass, and equal chart
  weights;

destroyed/not supplied:
  all scalar-cover roles outside the six displayed speeds, prescribed
  7- and 13-adic valuation ladders, source-owner and blocker typing,
  delayed semantic words, endpoint/cospan phase, and exact relation
  currents.                                                     (35)
```

In particular, (16) does **not** assert that the union of danger sets
of a full primitive fourteen-speed LRC row covers the circle.  It
does not install the two blocker truth bits, a THM-2305 pure/fork
word, a graft clock, or the coefficient phases needed by THM-2401.
The primitive control (26) removes only the cheap common-gcd objection.

The theorem instead establishes a stopping rule:

```text
finite open root-chart compatibility
```

cannot by itself exclude the uniform-offset locus.  A successful
closure must use a constraint which the mixed-radix construction does
not preserve -- most plausibly the full scalar cover, valuation/owner
ancestry, or one canonical endpoint phase.

## 6. Exact companion

Run

```text
python 04-computation/lrc14_mixed_radix_root_phase_orbit_thm2462.py
python -O 04-computation/lrc14_mixed_radix_root_phase_orbit_thm2462.py
```

The dependency-free `Fraction` companion:

- verifies both branches of the ordinary and guard nearest-integer
  laws for all `12*13` residue/start pairs;
- constructs all thirteen nested interval certificates for (16);
- directly enumerates the six physical root masks at every midpoint
  and checks every clean cover;
- proves the common width, disjointness, equal chart weights, and
  constant `4/13` guard density;
- verifies the full exact midpoint list (22);
- appends two independent divisible-by-thirteen factors and checks
  all twenty-six prescribed root-constant truth bits; and
- independently finds and checks all thirteen phase cells for the
  primitive tuple (26).

Both transcripts must reproduce

```text
05-knowledge/results/lrc14_mixed_radix_root_phase_orbit_thm2462.out
```

byte-for-byte after LF normalization.  Every truth-bearing executable
check uses explicit `require`, including under optimized Python.

## 7. Scope

This theorem proves a local fixed-speed physicalization of the
thirteen static root words and a general finite-open-cell realization
lemma.  It does not promote the provisional status of THM-2458, use
its static enumeration as a proved dependency, or prove that the
realized charts occur inside a valid global LRC counterexample.

No scalar profile is excluded, the ledger remains `165`, and LRC(14)
remains open.

QED.
