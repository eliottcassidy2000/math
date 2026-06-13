---
id: HYP-2107
status: algorithmic proof-program; sampled prototype, not a theorem
source: codex-2026-06-03-S581
related:
  - THM-398
  - HYP-2106
  - HYP-2105
  - HYP-2104
  - HYP-2103
  - HYP-2102
---

# HYP-2107: the large-owner Cprime residual should be a bounded CRT residue automaton in the multiplier w

## Claim

After THM-398, HYP-2104, and the S574 cover-to-congruence translator, the
remaining Cprime work should not return to prime-fiber enumeration.  Once the
long-interval criterion B' and Lemma C remove the easy cover failures, the
large-owner residual can be compiled into allowed residue classes for the
multiplier `w` in `v=nw`.

The proposed proof object is:

```text
safe component -> endpoint-owner congruence windows -> allowed residues for w
all components -> CRT intersection -> dominance-bounded residue set
```

If the resulting residue set has no positive representative below the dominance
cutoff

```text
w <= floor((n-1) max(S') / n),
```

then larger `w` are already loose by dominance and smaller `w` cannot align the
cover.  This gives a positive-measure certificate without exact prime fibers.

## Prototype Evidence

`04-computation/lrc_large_owner_residue_sieve_s581.py` samples primitive rows
with one multiple `n*w`, applies B' and Lemma C, then compiles the remaining
large-owner components into a CRT residue automaton.

Stored output:

```text
n total  B'   small  union  large-res  CRT-empty  bounded-empty  live  overflow
11  2500  92.2%  33.4%  94.4%        141   0.0%  75.2%   0.0%       35
12  2500  92.2%  35.5%  95.2%        119   0.0%  28.6%   0.0%       85
13  2500  93.3%  44.9%  95.8%        106   0.0%  20.8%   0.0%       84
14  2500  95.2%  57.9%  97.5%         63   0.0%   6.3%   0.0%       59
```

For exact-classified large-owner rows, the prototype found zero live bounded
`w` residue states for `n=11,12,13,14`.  The remaining cases are not live
counter-signals; they are modulus-overflow cases where the naive residue-set
implementation stopped before classifying the row.  Median exact moduli were
`23940, 42840, 45885, 88380`.

S599 follow-up: `lrc_twoblock_helly_s599.py` implements the promised minimal
contradiction extraction in the dominance-bounded window rather than the full
CRT modulus.  In deterministic samples after Bprime and Lemma C it finds zero
bounded-live rows; the remaining n=14 determinant residuals collapse to
singleton component walls, while n=6/n=8 include the first pair incompatibility
certificates.  This creates HYP-2144: prove the automaton emptiness by
singleton/pair/bounded-size Helly witnesses before resorting to a global CRT/ZDD.

## Why This Improves the Algorithm

S574 leaves an open residual: every component is short, and every component has
at least one large binding owner `u>=n`, so the endpoint congruence has slack.
The S581 extension observes that the slack is still periodic in `w`.  For one
endpoint owner `(u,k,eps)`, the relevant period divides

```text
u / gcd(u, k n + eps).
```

For one safe component, union the allowed residues over all left/right owner
pairs.  For many components, intersect these finite residue sets by CRT.  This
turns the residual from a geometric cover search into a finite automaton whose
states are residue classes and whose exits are:

```text
empty residue set,
no positive residue below dominance cutoff,
live bounded residue,
or modulus overflow.
```

The `w=0 mod M` class is a useful warning: pure all-`w` emptiness is too strong,
because some formal alignment classes survive at enormous multipliers.  The
dominance cutoff removes that ghost.  The correct proof target is bounded
residue emptiness, not global residue emptiness.

## Algorithmic Extensions

1. Replace the naive global residue set with a prime-factor CRT decision
   diagram.  Store constraints locally at each prime power and combine only
   when compatibility survives.
2. Extract certificates from empty intersections: a minimal incompatible subset
   of components should become a short proof, not a large table.
3. Split the overflow cases by modulus source.  A high modulus caused by many
   coprime large owners is evidence for even smaller density; a high prime-power
   modulus may need a local p-adic lemma.
4. Add dominance as a first-class automaton guard.  Residues whose first
   positive representative is above the cutoff should be discharged before any
   exact search.
5. Compile paper seed families into this automaton.  A failed prime-fiber seed
   should become a residue automaton clause, especially for total `n=12` where
   `C=23` leaves no nonunit shell residual.
6. Use meet-in-the-middle component intersection: intersect components into two
   residue summaries, then CRT-join the summaries.  This should replace the
   current overflow cap.
7. Look for a density theorem: if each component allows at most `alpha_i` of
   residues and the product is below the reciprocal dominance window, bounded
   emptiness follows without enumerating residue classes.

## Tournament Analysis

Vertices should not be runners.  Useful vertices for this stage are proof gates
and residue-automaton states:

```text
Bprime_long_interval, LemmaC_small_owner, bounded_CRT_empty,
live_CRT_state, modulus_overflow, exact_fallback
```

Pair observable:

```text
(unproved_fraction, modulus_growth, dependency_depth, -proof_strength, tie_order)
```

Switch: the harder residual burden beats the easier certified gate.  The S581
prototype remains transitive:

```text
score_hist: {0:1, 1:1, 2:1, 3:1, 4:1, 5:1}
directed_3_cycles: 0
sccs: 6 singleton SCCs
hamiltonian_path_count: 1
hardness_path: exact_fallback > modulus_overflow > live_CRT_state
               > bounded_CRT_empty > LemmaC_small_owner > Bprime_long_interval
```

Cycles should be sought only after component residue states are split by
prime-power factors; a cycle there would mean competing p-adic owner constraints,
not runner-level tournament complexity.

## Assumption Challenge

Candidate vertices considered: runners, gaps, safe components, endpoint owners,
owner pairs, residue classes of `w`, prime-power factors of the residue modulus,
and proof obligations.  The chosen quotient uses safe components and residue
automaton states.

Preserved predicate:

```text
Can the thin AP danger arcs of a multiple v=nw cover every component of G(S')?
```

Destroyed information: runner order, exact prime-fiber residue tuples, most
time-grid coordinates, and unneeded large values of `w` beyond the dominance
cutoff.  Challenged assumption: the large-owner residual is not a place to
return to exact `J(k,p)` enumeration; it is a bounded CRT feasibility problem in
one multiplier.

## Honest Status

This is a proof-program with a sampled prototype.  It does not prove Cprime or
LRC.  The next hard work is replacing the overflow cap with a factored CRT/ZDD
implementation and extracting minimal incompatible component certificates.

**See:** `04-computation/lrc_large_owner_residue_sieve_s581.py`
(+ `05-knowledge/results/lrc_large_owner_residue_sieve_s581.out`),
`07-reflections/lrc-large-owner-residue-automaton-s581.md`, THM-398,
HYP-2106, the two current HYP-2105 records, HYP-2104, HYP-2103, HYP-2102.
