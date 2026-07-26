---
id: THM-2441
title: "Septimal ancestry event-period collapse"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED; INDEPENDENT IMPLEMENTATION
  PENDING. For a rational finite-interval source with common endpoint
  denominator D=13^K D_0, every pair of reduced clocks congruent modulo
  7D_0 has pointwise ancestry profiles differing by one explicit uniform
  vector. Consequently the complete centred endpoint-event word, every
  charged septimal mode, the nonflat target mass, and THM-2421's exact
  target-restricted Gamma are periodic in the delayed base-thirteen
  scale, with period dividing ord_(7D_0)(13). The fixed target
  denominator does not enter the modulus. The uncentred word is not
  periodic. This is a finite exact reduction, not a landing theorem or
  an exclusion of any LRC(14) row.
source: codex-2026-07-26-ancestry-event-period-collapse
depends_on:
  - THM-2421-all-clock-septimal-ancestry-endpoint-event-detector
related:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-2409-unfiltered-septimal-source-completion-and-word-phase-boundary
  - THM-2418-alternating-base-thirteen-septimal-carry-matrix-and-rank-one-boundary
  - THM-2426-compositional-thirteen-root-final-septimal-lane-exclusion
script: 04-computation/lrc14_ancestry_event_period_collapse_thm2441.py
output: 05-knowledge/results/lrc14_ancestry_event_period_collapse_thm2441.out
script_sha256: d79a9c16353fc709ef00c21a3fe6e6bec665e8a5cd6a3881091726ff00d92b44
output_sha256: fb97f34d5165ceb1859858456d55c592430813929760cc18069cb7c3e4ab72a8
hash_basis: working-tree bytes (LF)
---

# THM-2441 -- delayed ancestry has a finite exact event period

**PROVED + VERIFIED-EXACT + HOSTILE-AUDITED; independent
implementation pending.**

THM-2421 replaces enumeration of all inverse branches at a clock `R` by
a signed endpoint-event word. Its Section 7 nevertheless leaves an
apparently infinite delayed-clock search:

```text
R=13^k -> push every source endpoint -> test the actual terminal word.
```

For a fixed rational source, that search is finite. The load-bearing
fact is stronger than periodicity of one scalar: congruent reduced
clocks change the complete ancestry vector by a single known vector on
the constant line.

## 1. Setup and statement

Let

```text
E=disjoint_union_(i=1)^m [A_i/D,B_i/D)
```

be a finite union of half-open intervals in `[0,1)`, where all
`A_i,B_i` are integers and

```text
D=13^K D_0,                 gcd(D_0,13)=1.                  (1)
```

Put

```text
S=sum_i(B_i-A_i)=D measure(E).                              (2)
```

For a positive integer `N`, use the clock

```text
R=13^K N
```

and define THM-2421's septimal inverse-branch ancestry vector

```text
mathcal N_R(y)=(N_(R,0)(y),...,N_(R,6)(y)),

N_(R,c)(y)
 =#{0<=n<R:
      n=c mod 7 and (y+n)/R belongs to E}.                  (3)
```

The definitions are pointwise on the half-open terminal coordinate
`0<=y<1`; no almost-everywhere qualification is needed for the main
identity.

### The uniform-drift theorem

Suppose

```text
N'=N+7D_0 t,                         t integer,
R'=13^K N'.                                                 (4)
```

with `N,N'>0`. Then, for every `0<=y<1`,

```text
mathcal N_(R')(y)
 =mathcal N_R(y)+tS (1,1,1,1,1,1,1).                       (5)
```

In particular, both natural centred profiles are unchanged:

```text
mathcal N_R^glob(y)
 =mathcal N_R(y)-R measure(E)/7 (1,...,1),                  (6)

mathcal N_R^pt(y)
 =mathcal N_R(y)
   -(sum_c N_(R,c)(y))/7 (1,...,1).                         (7)
```

For every fixed `y`,

```text
mathcal N_(R')^glob(y)=mathcal N_R^glob(y),
mathcal N_(R')^pt(y)=mathcal N_R^pt(y).                     (8)
```

The global centring in (6) is generally rational. The pointwise
projection in (7) is the orthogonal projection away from the constant
line. Their equality across the two clocks, not integrality of the
centred coordinates, is the assertion.

## 2. Exact residue-count proof

Fix one component `[A/D,B/D)` of `E`. At terminal phase `y`, its
contributing prefixes are the integers in

```text
I_(A,B)(N,y)
 =[ceil(NA/D_0-y), ceil(NB/D_0-y)) intersect Z.             (9)
```

The apparent restrictions `0<=n<R` add nothing to (9). Indeed,
`0<=A<B<=D`; the lower endpoint in (9) is at least zero and the
upper endpoint is at most `R`, including at `y=0`.

Under (4), the two integer endpoints in (9) change by

```text
ceil(N'A/D_0-y)
 =ceil(NA/D_0-y)+7tA,

ceil(N'B/D_0-y)
 =ceil(NB/D_0-y)+7tB.                                      (10)
```

For integers `L,U,alpha,beta` and a residue `c mod 7`, exact
quotient-and-remainder counting gives

```text
#{L+7alpha<=n<U+7beta:n=c mod 7}

 =#{L<=n<U:n=c mod 7}+beta-alpha.                           (11)
```

Apply (11) to (10). Every residue count contributed by this component
increases by

```text
t(B-A).
```

Summing over the disjoint components proves (5).

Moreover,

```text
(R'-R) measure(E)/7
 =(7Dt)(S/D)/7=tS,                                         (12)
```

which proves the globally centred identity. Taking coordinate sums in
(5) proves the pointwise-centred identity.

This proof also covers negative `t`: swap the two clocks, or use the
floor-difference form of (11).

## 3. The full signed endpoint word is unchanged

Let `x=A/D` be a source endpoint, retaining its left/right sign. By
(4),

```text
R'x=Rx+7tA.                                                 (13)
```

Therefore the following data are identical at the two clocks:

```text
terminal phase          {Rx},
septimal label          floor(Rx) mod 7,
endpoint sign,
collision class with every other endpoint,
whether the event lies at the cut {Rx}=0.                  (14)
```

Thus the complete signed endpoint-event multiset is unchanged. Equation
(5) supplies the part which endpoint events alone do not retain: its
base chamber profile changes only by `tS(1,...,1)`. Consequently the
**centred endpoint-event word** -- signed labelled events together with
one centred base profile -- is exactly the same.

The word may be restricted to any fixed measurable target `Q` without
changing this conclusion. If `Q` is a rational finite-interval set,
inserting its fixed endpoints gives the same finite relative chamber
word at both clocks. The denominator of `Q` does not enter (4): after
inverse-branch disintegration, `Q` lives in the fixed terminal
coordinate `y`.

## 4. Exact power period

Now take the delayed base-thirteen clocks

```text
R_d=13^(K+d),                        d>=0.                  (15)
```

Because `13` is a unit modulo `7D_0`, define

```text
T=ord_(7D_0)(13).                                           (16)
```

Then

```text
13^(d+T)=13^d mod 7D_0,
```

so (5)--(14) give, for every `d>=0`,

```text
mathcal N_(R_(d+T))^glob
 =mathcal N_(R_d)^glob,

mathcal N_(R_(d+T))^pt
 =mathcal N_(R_d)^pt.                                      (17)
```

The exact period can be a proper divisor of `T`; (16) is a uniform
period bound.

THM-2418's physical-root relabelling

```text
ell=c+(-1)^(K+d) r mod 7                                  (18)
```

does not enlarge the bound. Reduction modulo seven shows that `T` is
even, since `13=-1 mod 7`. Hence `d` and `d+T` have the same reflection
in (18).

## 5. Every charged detector is periodic

For `e in F_7^*`, let

```text
Nhat_(R,e)(y)=sum_(c=0)^6 N_(R,c)(y) zeta_7^(-ec).
```

The added vector in (5) has no charged Fourier component. Therefore

```text
Nhat_(R',e)(y)=Nhat_(R,e)(y)                               (19)
```

pointwise, with its complex phase retained.

For any fixed measurable target `Q`, THM-2421 defines

```text
Gamma_Q(E;R)
 =integral_Q sum_c
    (N_(R,c+1)(y)-N_(R,c)(y))^2 dy                         (20)
```

and the nonflat target mass

```text
eta_Q(E;R)
 =measure{y in Q:mathcal N_R(y) is nonconstant}.            (21)
```

Equations (5), (19) immediately imply

```text
Gamma_Q(E;R')=Gamma_Q(E;R),
eta_Q(E;R')=eta_Q(E;R).                                    (22)
```

Thus (20)--(22), every charged pointwise mode, and the complete centred
event word all depend on `d` through only

```text
13^d mod 7D_0.                                              (23)
```

## 6. Centring is essential

The raw ancestry profile is not periodic. Take

```text
E=[0,1),                 Q=[1/3,2/3).
```

At every terminal phase,

```text
mathcal N_13
 =(2,2,2,2,2,2,1),

mathcal N_(13^3)
 =(314,314,314,314,314,314,313).                           (24)
```

The clocks are congruent modulo seven, but the two raw vectors differ by

```text
312(1,...,1).
```

Their centred vectors agree and

```text
Gamma_Q(E;13)=Gamma_Q(E;13^3)=2/3.                          (25)
```

This is the minimal information-loss boundary: the signed endpoint
events agree, but an uncentred base profile retains the irrelevant
number of complete seven-sheet blocks.

The period need not collapse to one. For

```text
E=[1/39,8/39) union [13/39,17/39),
Q=[1/5,4/5),
```

one has `K=1,D_0=3` and

```text
ord_21(13)=2.
```

The exact Gamma values at `R=13^(1+d)`, for `d=0,...,5`, are

```text
12/5, 16/15, 12/5, 16/15, 12/5, 16/15.                    (26)
```

Hence the alternating object can remain genuinely alternating even
though its infinite tail has collapsed to a finite exact classifier.

## 7. Consequence for the surviving LRC(14) graft

In the surviving branch of THM-2426/THM-2436,

```text
nu_7(c_3)<=M,
```

an ordinary `q_*` of depth `M` makes `(c_3,q_*)` an isolated
noncirculant target graft by THM-2367. THM-2365 says that a positive bare
owner drift survives sufficiently delayed words, but the missing
statement is owner-conditioned no-cancellation on the actual exclusive
source and terminal word.

For any fixed rational scalar-cover packet `(E_j,Q_(j,sigma))`, the
present theorem reduces its complete delayed-clock endpoint audit to

```text
at most ord_(7D_(0,j))(13) exact phase classes,             (27)
```

each computable from one base chamber and the signed source endpoints.
No enumeration of `13^k` branches and no unbounded scale search is
needed. The output must retain the centred signed word, the source
owner, the target graft, and THM-2418's parity/reflection; `Gamma` alone
still loses orientation and terminal landing.

This does **not** provide an explicit genuine global
`E_j,Q_(j,sigma)` packet, prove a linear current survives, remove a
scalar row, or prove LRC(14). It makes the cheapest decisive computation
finite once that packet is supplied.

## 8. Exact companion

The dependency-free companion performs:

1. an exhaustive audit of all `510` unions of grid cells through
   denominator eight, at `4,080` congruent clock pairs;
2. an exhaustive audit of all `1,222` single intervals on the grids
   `D=13,26,39`, at `2,444` congruent clock pairs;
3. a power-clock audit for `D_0=1,...,12` and `K=0,1`, visiting every
   residue class in all `24` exact cycles;
4. direct branch enumeration independently checking the residue-count
   evaluator on every deterministic modest-clock chamber;
5. the raw hostile (24)--(25), the genuine period-two positive control
   (26), a target whose denominator is absent from every source
   modulus, and `96` secondary random exact fixtures.

In total it checks `6,724` congruent-clock pairs, `92,306` endpoint and
chamber probes, `177,552` direct branch profiles, and the exact endpoint
word and target-restricted `Gamma` at every pair. Reproduce with

```bash
python3 04-computation/lrc14_ancestry_event_period_collapse_thm2441.py
python3 -O 04-computation/lrc14_ancestry_event_period_collapse_thm2441.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_ancestry_event_period_collapse_thm2441.out
```

after LF normalization. Every load-bearing check remains active under
optimized Python. QED.
