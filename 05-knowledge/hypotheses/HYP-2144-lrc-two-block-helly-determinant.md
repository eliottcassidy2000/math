---
id: HYP-2144
status: OPEN proof-program; S599 gives bounded-window evidence and first small certificates
source: codex-2026-06-03-S599
related:
  - HYP-2142
  - HYP-2140
  - HYP-2107
  - HYP-2110
  - THM-398
  - THM-401
  - THM-404
---

# HYP-2144: the bounded two-block determinant automaton should have small Helly witnesses

## Claim

HYP-2142 identifies the large-owner Cprime residual as a bounded feasibility
problem in one multiplier `w`: each component of `G(S')` contributes a
two-block determinant language, and a cover exists only if all component
languages intersect below the dominance cutoff

```text
1 <= w <= floor((n-1) max(S') / n).
```

The next conjectural strengthening is that this bounded intersection usually,
and perhaps always in the true residual, has a small certificate of emptiness:
a singleton determinant wall, a pair of incompatible determinant rows, or a
bounded-size Helly witness.  The human proof object should be the minimal empty
subfamily, not the full CRT modulus.

## S599 Evidence

`04-computation/lrc_twoblock_helly_s599.py` enumerates only the dominance-bounded
window in `w`, avoiding the full CRT lcm.  For each component it unions endpoint
owner choices and records the allowed bounded multipliers.  When the global
intersection is empty, it extracts a minimal empty subset up to size four.

In the S581/S595 large-owner regime after Bprime and Lemma C (`BC_only`), the
deterministic sample found no bounded-live rows:

```text
n  rows  preempt  singleton  pair  triple  quad  high  live
6  1195      897        286    12       0     0     0     0
8  1200     1002        196     2       0     0     0     0
10 1200     1111         89     0       0     0     0     0
12 1200     1151         49     0       0     0     0     0
14 1200     1185         15     0       0     0     0     0
```

The first n=14 post-B/C determinant residual in the sample is already a singleton
wall:

```text
row=(2,3,6,7,8,11,13,14,17,19,20,22,23), w=1, w_bound=21
C0 length=3/644, owners=sL, allowed={}
```

The first genuine two-row Helly example appears at n=6:

```text
row=(1,12,13,14,15), w=2, w_bound=12
C0 allowed={1}
C1 allowed={2,4}
```

The full stack regime first removes the S598 total/one-sided caps and Lemmas
E/F.  It still finds no bounded-live rows in the sample; n=14 is completely
preempted by earlier gates, while smaller n retain singleton/pair determinant
witnesses:

```text
n  rows  preempt  singleton  pair  triple  quad  high  live
6  1194     1118         66    10       0     0     0     0
8  1200     1187         10     3       0     0     0     0
10 1200     1200          0     0       0     0     0     0
12 1200     1199          1     0       0     0     0     0
14 1200     1200          0     0       0     0     0     0
```

The first full-stack pair example is all-large-owner:

```text
row=(2,6,8,10,13), w=1, w_bound=10
C0 allowed={1}
C1 allowed={3,7,10}
```

## S601 Helly-Log Follow-Up

HYP-2152 and `04-computation/helly_scale_log_laws_s601.py` turn the informal
"Helly scale" phrase into a certificate entropy:

```text
Lambda_H(M) = log sum_{h<=H} binom(M,h),
```

where `M` is the number of live component languages and `H` is the Helly rank
cutoff.  For fixed `H`, `Lambda_H(M) = H log M + O_H(1)`, so log, loglog, and
logloglog factors come from how compressed the live component count already is.

The S601 deterministic sample (`1800` row attempts for each `n` and regime)
found `1113` extracted empty certificates, all rank one or two:

```text
h=1: 1084
h=2:   29
h>=3:   0
high_order_empty=0
bounded_live=0
```

This strengthens the proof-program interpretation: the sampled two-block
determinant residual is earning a bounded local Helly dividend before the global
CRT modulus appears.  In particular, the fresh `BC_only` n=14 sample had
`25` singleton walls and no pair/high/live rows; the full-stack n=14 sample was
fully preempted.

Rebase integration: while S599 was being pushed, monad-compute S598 added a
widened exhaustive Cprime box and an n=9 run, checking 6.24M configs with zero
tight cases.  That makes the computational Cprime backdrop stronger.  HYP-2144
is not mainly a volume check; it is an attempt to extract the human-scale
determinant subcertificate inside that already-stable residual.

## Why This Pushes HYP-2142

S581 saw zero live bounded states among exact-classified large-owner rows, but
some cases were hidden behind CRT modulus growth.  S599 asks a different
question: before building a global CRT object, can one already find a small
empty subfamily in the bounded window?

This shifts the proof target from

```text
all component languages have empty bounded intersection
```

to

```text
some small set of determinant rows has empty bounded intersection.
```

That is a Helly-style strengthening of the automaton.  It does not contradict
HYP-2042's warning that the full LRC clearance cascade has no fixed low Helly
number: the quotient here is not raw runner safe sets.  It is the post-Cprime,
post-owner, dominance-bounded determinant language.

## Proof Route

1. Prove a singleton wall criterion: characterize when one component has no
   allowed bounded `w`.
2. Prove a two-row incompatibility lemma: if two components force disjoint
   bounded arithmetic progressions in `w`, the row is loose.
3. Factor the pair lemma by prime powers.  This should explain when n=14
   collapses to singleton walls and when n=6/8 still require pairs.
4. Only after singleton/pair extraction fails should the proof build the full
   CRT/ZDD object proposed in HYP-2107.

## Tournament Analysis

Vertices are determinant certificate sizes:

```text
singleton_empty, pair_empty, triple_empty, quad_empty,
high_order_empty, bounded_live, preempted_gate
```

Pair observable:

```text
(certificate_rank, sampled_route_count, name)
```

Switch/gauge: smaller empty determinant subfamily beats larger/global residue
burden.  S599's fingerprint is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3_cycles=0
SCCs=7 singleton SCCs
Hamiltonian_paths=1
tie_HP=singleton_empty > pair_empty > triple_empty > quad_empty
       > high_order_empty > bounded_live > preempted_gate
```

## Assumption Challenge

Candidate vertices considered: runners, gaps, cap centers, safe components,
endpoint owners, owner pairs, residue classes, prime powers, and proof
obligations.  The chosen vertices are determinant component languages and
minimal empty subsets.

Preserved predicate:

```text
Can the bounded multiplier window contain a common w that covers every safe
component by v=nw?
```

Destroyed information: phase order beyond the bounded window, the full CRT
modulus, exact large multipliers already discharged by dominance, and runner
order not visible through endpoint-owner languages.  Challenged assumption: the
automaton emptiness proof need not start from a global modulus; small Helly
witnesses may be the correct proof objects.

## Honest Status

This is sampled evidence and a proof-program, not a theorem.  The script
enumerates bounded `w` directly, so it is deliberately not a scalable CRT
engine.  Its value is certificate extraction and concept-sharpening: it shows
where singleton and pair determinant walls appear, and it gives the next formal
lemmas to try before returning to a factored CRT/ZDD implementation.

**See:** `04-computation/lrc_twoblock_helly_s599.py`
(+ `05-knowledge/results/lrc_twoblock_helly_s599.out`),
`07-reflections/lrc-two-block-determinant-helly-s599.md`, HYP-2142,
HYP-2107, HYP-2140, THM-398, THM-401, THM-404.
