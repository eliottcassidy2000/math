---
id: THM-2322
title: "Local hostile c3-toothpick ladder"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In the exact local
  multiplier-four carrier of THM-2299 and THM-2318, every positive odd
  owner-normalized frequency is marked by both the local source and its
  exact c_2-only terminal current. Its grade-three, root-character-four
  slice splits into thirteen parallel bi-infinite arithmetic ladders whose
  rungs differ by c_3=2*13^5. One positive ray is
  N_m=13^3(17+338m), m>=0. Deleting one septimal congruence class from
  each rail leaves thirteen 91-unit local arithmetic paths whose edge
  multipliers are one or two. This is a theorem about the explicit local
  subcarrier E, not the full canonical exclusive set E_j: it does not prove
  incidence in THM-2293's canonical atom graph, exclude a scalar row, or
  prove LRC(14).
source: codex-2026-07-25-local-hostile-toothpick-ladder
depends_on:
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
  - THM-2306-owner-normalized-disjoint-support-and-first-collision-shell
related:
  - THM-2293-quadratic-root-energy-raises-the-ancestry-shell
  - THM-2313-biprime-pareto-collision-frontier-and-91-unit-current-shell
  - THM-2318-one-shot-three-prime-mobius-amplifier
script: 04-computation/lrc14_local_hostile_c3_toothpick_ladder_thm2322.py
output: 05-knowledge/results/lrc14_local_hostile_c3_toothpick_ladder_thm2322.out
script_sha256: 3c0f22b16adff63468909b00676b0cf48c069e1a7df2244359dccb0efb09224c
output_sha256: 1537fef871336b93852190a4917e8b9eabb460124d99b74fc5da6561b13d05e1
hash_basis: working-tree bytes (LF)
---

# THM-2322 -- the local hostile carrier contains an infinite marked toothpick ladder

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2318 uses a fully minimal collision cube to produce one common
source/current atom at grade three with a `91`-unit outside multiplier.
The cube does not prescribe its root character or place it on an incident
`c_3` edge.

For a general LRC source, those losses remain exact. The local hostile
carrier used as THM-2318's positive control has more rigid geometry,
however: both its source and current are equal-width pairs of symmetric
intervals. Their Fourier zero sets can be read exactly. The result is not
one atom but a self-similar arithmetic ladder:

```text
same local source and exact terminal current
  -> every rung has root character four and grade three
  -> one c_3 step advances to the next rung
  -> no rung is lost to component-phase cancellation.                (1)
```

This separates a failure of the general collision selector from a failure
of the underlying local object. The selector forgets residue and incidence;
the explicit local object actually has both at the local arithmetic level.

## 1. The exact local carrier

Use the multiplier-four profile from THM-2299:

```text
(c_1,c_2,c_3)=(13,13^3,2*13^5),
epsilon=10^(-12).                                  (2)
```

On the circle, let

```text
F=
 (-1/16-epsilon,-1/16+epsilon)
 union
 ( 1/16-epsilon, 1/16+epsilon),                    (3)

E={(8+y)/13:y in F}.                               (4)
```

This is the exact **local** source subcarrier from THM-2299. It lies in the
`c_1`-exclusive region, but it is not asserted to equal the full canonical
exclusive set `E_1`.

Let `P=P_13`. The two prescribed images of `F` are

```text
R=T(F):
  centers {3,13}/16, half-width 13epsilon,

C=T(R):
  centers {7,9}/16, half-width 169epsilon.          (5)
```

Every point of `R` has the same `c_2`-only terminal word. Consequently the
local exact word current is

```text
f=1_R P^2 1_E=P^2 1_E=(1/169)1_R.                  (6)
```

THM-2306 records the two owner-normalized functions

```text
B=P1_E=(1/13)1_F,
A=Pf=(1/2197)1_C.                                  (7)
```

Fix the Fourier convention

```text
h_hat(q)=integral_T h(x)exp(-2*pi*i*q*x)dx.         (8)
```

Perron transport gives, for every integer `q`,

```text
(1_E)_hat(13q)=B_hat(q)=(1/13)(1_F)_hat(q),
f_hat(13q)=A_hat(q)=(1/2197)(1_C)_hat(q).           (9)
```

Thus a common nonzero Fourier coefficient of `1_F` and `1_C` at `q`
becomes a common local-source/current coefficient at `13q`.

There is also a direct current formula. The centers of `R` are
`plus or minus 3/16`, so for every nonzero integer `N`,

```text
f_hat(N)
 =2 sin(26*pi*N*epsilon)/(169*pi*N)
    cos(3*pi*N/8).                                  (9a)
```

The normalization preserves mass exactly:

```text
measure(E)=4epsilon/13,
measure(R)=52epsilon,
integral f=52epsilon/169=4epsilon/13.               (9b)
```

## 2. Exact Fourier zero sets

For every nonzero integer `q`, direct integration of the two symmetric
interval pairs gives

```text
(1_F)_hat(q)
 =2 sin(2*pi*q*epsilon)/(pi*q) cos(pi*q/8),         (10)

(1_C)_hat(q)
 =2 sin(338*pi*q*epsilon)/(pi*q) cos(7*pi*q/8).    (11)
```

Both formulas are exact. If `q` is odd, neither cosine vanishes:

```text
cos(pi*q/8)=0       iff q=4 mod 8,
cos(7*pi*q/8)=0     iff 7q=4 mod 8
                    iff q=4 mod 8.                 (12)
```

The sine factors in (10)--(11) have the same integer zero condition. Since

```text
gcd(338,10^12)=2,
```

we have

```text
sin(2*pi*q*epsilon)=0
 iff 5*10^11 divides q
 iff sin(338*pi*q*epsilon)=0.                      (13)
```

The divisor in (13) is even. Therefore every nonzero odd integer `q`
satisfies

```text
(1_F)_hat(q)!=0,
(1_C)_hat(q)!=0.                                   (14)
```

Equations (9) and (14) prove the stronger **odd parity spine**

```text
(1_E)_hat(13q)!=0,
f_hat(13q)!=0                    for every nonzero odd q. (14a)
```

The current assertion can also be checked directly from (9a): `13q` is
odd, so its cosine cannot vanish, while a sine zero would require the odd
integer `13q` to be divisible by the even integer `5*10^11`.

## 3. The infinite c3 ladder

The grade-three, root-character-four frequencies with odd outside
multiplier are exactly

```text
N_k=13^3(17+26k),                         k in Z.   (15)
```

Indeed, an integer congruent to four modulo thirteen is odd exactly when
it has the form

```text
4+13(2k+1)=17+26k.
```

For (15), use (14a) with

```text
q=13^2(17+26k).
```

It follows that every member of (15) is simultaneously nonzero for the
local source and current. Split `k` modulo thirteen:

```text
r_(a,m)=17+26a+338m,
q_(a,m)=13^2 r_(a,m),
N_(a,m)=13^3 r_(a,m),

0<=a<=12,                              m in Z.      (16)
```

This partitions the whole odd family (15) into thirteen parallel
bi-infinite rails. On each rail,

```text
nu_13(N_(a,m))=3,
13^(-3)N_(a,m)=4 mod 13,                           (17)

(1_E)_hat(N_(a,m))!=0,
f_hat(N_(a,m))!=0,                                 (18)

N_(a,m+1)-N_(a,m)
 =13^3*338
 =2*13^5
 =c_3.                                             (19)
```

Thus the local common spectrum contains thirteen infinite arithmetic paths in the
exact grade-three, root-character-four shell, with every consecutive edge
equal to one deepest-blocker step. Every vertex is doubly marked:

```text
local source mark:          (1_E)_hat(N_(a,m))!=0,
exact c_2-current mark:     f_hat(N_(a,m))!=0.      (20)
```

For a concrete positive ray, take `a=0` and `m>=0`:

```text
r_m=17+338m,
q_m=13^2 r_m,
N_m=13q_m=13^3 r_m.                                (21)
```

The first two vertices are

```text
N_0=13^3*17=37349,
N_1=13^3*355=779935,
N_1-N_0=c_3=742586.                                (22)
```

Both outside multipliers `17` and `355` are coprime to `91`. In fact,

```text
r_m=0 mod 7 iff m=2 mod 7.                         (23)
```

Delete precisely those indices and order the survivors. Their successive
index gaps are one or two, so the corresponding frequency gaps are

```text
c_3 or 2c_3.                                       (24)
```

Both edge multipliers are coprime to `91`, and every surviving vertex has
a `91`-unit outside multiplier. This gives an infinite `91`-unit local
path, not just the first toothpick (22).

The same rule works on every rail: since `338=2 mod 7`, the deleted class
on rail `a` is exactly

```text
m=a+2 mod 7.                                       (24a)
```

The remaining local arithmetic edge multipliers are again one or two.

## 4. What this repairs, and what remains open

For the exact hostile positive control, THM-2318's three stated missing
coordinates refine as follows.

```text
root residue:
  repaired locally -- every N_m has prescribed character kappa=4;

c_3 incidence:
  repaired at the local arithmetic level -- consecutive rungs differ by
  c_3, and both endpoints retain both the local source and exact
  word-current marks;

target word:
  preserved locally -- f is the literal c_2-only terminal current;

canonical incidence:
  still open -- E is only a local subcarrier and is not identified with
  the full canonical E_j;

global scalar consequence:
  absent -- the carrier is not a global scalar cover and excludes no row.
                                                               (25)
```

In particular, (18)--(24) do **not** imply that any rung is a vertex of
THM-2293's canonical graph

```text
Gamma_(j,kappa)
 ={N:(1_(E_j))_hat(N)!=0,
      nu_13(N)=b,
      N/13^b=kappa mod 13}.                         (26)
```

Fourier support is not monotone under set inclusion: a coefficient of
`1_E` need not survive after the rest of `E_j` is added. Therefore this
theorem neither proves a marked canonical edge nor weakens the
same-index/incidence target left by THM-2319 and THM-2321.

The positive lesson is structural. The two symmetric interval pairs have
the same parity-protected Fourier zero set, while `c_3/13^3=338` is even.
Translation by one deepest-blocker step preserves oddness and the root
residue. This parity-plus-filtration mechanism is the exact
**toothpick self-similarity** behind the ladder.

No scalar profile is eliminated, and LRC(14) remains open.

## 5. Exact verification

The companion checks the two Perron image center/radius/mass ledgers, the
arithmetic derivation of (10)--(14a), `100000` initial odd-spine points and
ladder rungs, all thirteen parallel rails, the exact grade/root/step
identities, the first marked edge, the septimal deletion law, and the
one-or-two-step `91`-unit survivor path. Every load-bearing test raises
explicitly in ordinary and optimized Python.

Reproduce with

```bash
python3 04-computation/lrc14_local_hostile_c3_toothpick_ladder_thm2322.py
python3 -O 04-computation/lrc14_local_hostile_c3_toothpick_ladder_thm2322.py
```

The two transcripts must match

```text
05-knowledge/results/lrc14_local_hostile_c3_toothpick_ladder_thm2322.out
```

byte-for-byte after LF normalization.
