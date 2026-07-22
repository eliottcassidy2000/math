---
id: THM-2080
title: "The mixed 1/14--1/7 comb overlap is an exact two-fold law, and every rank-six terminal guard has a reduced resonance of product at most 36"
status: >
  PROVED. An exact Fourier--Bernoulli calculation gives the Haar overlap of
  one 1/14 danger comb with one 1/7 guard comb. Consequently, containment of
  a rank-six terminal safe set in one guard forces some reduced ratio q/h=a/b
  with ab<=36. This is a scale-free finite projective gate, not a proof of
  LRC(14) or a classification of the remaining resonance directions.
source: codex-2026-07-22-LRC-terminal-resonance
depends_on:
  - THM-2073
  - THM-2076
related:
  - THM-965
  - THM-1232
  - THM-1233
  - THM-2077
  - THM-2078
script: 04-computation/lrc14_mixed_guard_fold_referee_codex_20260722.py
output: 05-knowledge/results/lrc14_mixed_guard_fold_referee_codex_20260722.out
---

# THM-2080 -- mixed guard fold and the deepest terminal resonance gate

Put

```text
D_q={t in R/Z: ||qt||<1/14},
E_h={t in R/Z: ||ht||<1/7}.                              (1)
```

For an integer `r`, let `rbar` be its representative in `{0,...,13}` and

```text
F(r)=rbar(14-rbar).                                      (2)
```

Thus `0<=F<=49` and `F` is even under `r -> -r` modulo `14`.

## 1. Exact mixed-radius overlap

Let `d=gcd(q,h)`, `a=q/d`, and `b=h/d`. Then

```text
measure(D_q intersect E_h)
 = [8ab+F(b+2a)-F(b-2a)]/(196ab).                       (3)
```

In particular,

```text
measure(D_q intersect E_h)
 <=2/49+1/(4ab)
 = 2/49+gcd(q,h)^2/(4qh).                               (4)
```

### Proof

The substitution `u=dt` is Haar-measure preserving, so it is enough to use
the coprime frequencies `a,b`. The Fourier coefficients of the two centered
interval indicators are

```text
Dhat(n)=sin(pi n/7)/(pi n),
Ehat(n)=sin(2pi n/7)/(pi n)                             (5)
```

for nonzero `n`, with constant coefficients `1/7` and `2/7`. Orthogonality
leaves exactly the solutions `na+mb=0`, namely `(n,m)=(bk,-ak)`. Therefore

```text
measure(D_a intersect E_b)
 =2/49 + 2/(pi^2 ab) sum_(k>=1)
          sin(pi bk/7) sin(2pi ak/7)/k^2.               (6)
```

Use `2 sin x sin y=cos(x-y)-cos(x+y)` and

```text
sum_(k>=1) cos(2pi kx)/k^2=pi^2 B_2({x}),
B_2(x)=x^2-x+1/6.                                      (7)
```

For an integer `r`,

```text
B_2({r/14})=1/6-F(r)/196.                              (8)
```

Substituting (7)--(8) into (6) proves (3). The two folds lie in `[0,49]`,
so their signed difference is at most `49`; this proves (4). The open
endpoints have Haar measure zero and do not affect the calculation. QED.

This is the mixed-width companion to THM-965. The baseline `2/49` is the
independent product of the two comb densities; all correlation is carried by
the two residues `b-2a` and `b+2a` modulo `14`.

## 2. A necessary resonance invoice for any terminal guard

Let `Q` have `s` distinct positive speeds and suppose

```text
G_Q={t: ||qt||>=1/14 for every q in Q} subset E_h.       (9)
```

Then `E_h` is covered by the `s` open danger combs `D_q`: if a point of
`E_h` lay in none of them, it would belong to `G_Q`. The union bound and
(3) give the exact signed-fold necessity

```text
sum_(q in Q) [F(b_q+2a_q)-F(b_q-2a_q)]/(a_q b_q)
 >=8(7-s),                                             (10)
```

where `a_q=q/gcd(q,h)` and `b_q=h/gcd(q,h)`. Using (4) instead gives the
coarser but transparent resonance invoice

```text
sum_(q in Q) gcd(q,h)^2/(qh) >=8(7-s)/49.              (11)
```

These inequalities are informative exactly below rank seven. This matches
the Haar wall in THM-2076: scalar density alone reaches rank six and then
needs arithmetic correlation.

## 3. The depth-five terminal lane has a finite projective gate

In a depth-five THM-2073 tower, the terminal core `Q=Q_5` has six speeds and
the final odd guard `h=h_4` satisfies (9). Equation (11) becomes

```text
sum_(q in Q) gcd(q,h)^2/(qh) >=8/49.                   (12)
```

Hence one of the six summands is at least `4/147`. For that `q`, write the
reduced ratio as

```text
q/h=a/b,                 gcd(a,b)=1.                  (13)
```

Since `gcd(q,h)^2/(qh)=1/(ab)`, one has

```text
ab<=147/4,
and therefore                 ab<=36.                 (14)
```

Thus the unbounded deepest lane cannot be projectively generic: at least one
terminal speed lies on one of the finitely many coprime ratio directions
`a/b` with `ab<=36`. Because the tower guard is odd, its reduced denominator
`b` is odd. The conclusion is scale-free and survives arbitrarily large
terminal maximum.

The exact fold version (10) is stronger than (14) and should be retained in
any case split. It identifies the signed mod-14 beat/sum residue, not merely
the small product.

## 4. Finite audit and the remaining conductor

The companion independently computes both sides of (3) by exact midpoint
atom counts for all `1<=q<=40` and odd `1<=h<=39`. It also replays the
rank-six part of the THM-2078 box. Among the six hereditarily primitive,
divisor-complete cores through height `24` and all `144` guards allowed by
the THM-2077 height bound, the first-order exact overlap test leaves only
`28` pairs. THM-2078's full simultaneous bitset test removes all `28`.

The strict containment problem has therefore split into two genuinely
different layers:

```text
generic projective lane: closed by (14),
finite resonance directions ab<=36: need simultaneous endpoint/phase
                                     transport, not another scalar sum.    (15)
```

THM-1232/1233 show how component-span ledgers compactify a six-comb cover
once a located carrier gap and its ordered ratios are retained. They do not
apply verbatim here: `E_h` consists of guard-danger teeth, not complete
`h`-safe gaps, and terminal speeds may lie on either side of `h`. The next
valid target is a mixed-tooth analogue on each of the finite directions in
(14), keeping the tooth address and endpoint owners supplied by THM-2075.

## 5. Assumption challenge and Tournament Analysis

The challenged assumption is that an unbounded terminal maximum requires an
unbounded clock search. The correct first quotient is the reduced pair
direction `(a,b)` together with the signed residues `b+2a,b-2a mod 14`. It
preserves the exact pair correlation and destroys the common scale `d` and
the simultaneous six-comb endpoint word; that lost word is the mandatory
sidecar for the residual.

There is no faithful tournament on runners. Orienting the six speeds by the
size of their overlap with `E_h` gives a transitive priority order, but the
cover predicate is a union and can be decided by phase collisions invisible
to pairwise dominance. The useful objects are the six signed fold invoices,
the guard tooth address, and the simultaneous endpoint-cover hypergraph.
QED.
