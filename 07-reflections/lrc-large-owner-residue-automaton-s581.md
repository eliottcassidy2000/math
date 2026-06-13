---
source: codex-2026-06-03-S581
status: speculative algorithmic synthesis plus sampled prototype
tags:
  - lonely-runner
  - Cprime
  - large-owner
  - CRT
  - residue-automaton
  - endpoint-owner
  - algorithmic-math
---

# Large-owner residue automata after the cover-to-congruence translator

The new algorithm has a clearer shape after S574.  The old fallback was too
coarse:

```text
all-short residual -> maybe exact finite prime fibers
```

The better fallback is:

```text
all-short residual
-> endpoint-owner congruence windows
-> allowed w-residue automaton
-> dominance-bounded emptiness
-> only then exact fallback
```

The trick is that S574's large-owner slack still depends periodically on the
multiplier `w` in `v=nw`.  For endpoint owner `(u,k,eps)`, the residue period is
controlled by `u / gcd(u, k n + eps)`.  A component contributes a small periodic
language of allowed `w` classes.  A cover attempt is the intersection of these
languages over all components of `G(S')`.

## What the prototype says

`lrc_large_owner_residue_sieve_s581.py` applies:

```text
B' long component
Lemma C small-small owner
large-owner residue automaton
dominance cutoff
```

On sampled primitive rows:

```text
n=11: B' or Lemma C proves 94.4%; large residual 141; no live bounded CRT states among exact-classified rows
n=12: B' or Lemma C proves 95.2%; large residual 119; no live bounded CRT states among exact-classified rows
n=13: B' or Lemma C proves 95.8%; large residual 106; no live bounded CRT states among exact-classified rows
n=14: B' or Lemma C proves 97.5%; large residual 63;  no live bounded CRT states among exact-classified rows
```

The prototype still overflows on many large-owner rows because it materializes
residues directly.  That is an implementation defect, not a mathematical
counter-signal.  The exact-classified rows say the live residual is thinner than
"large owner" suggests.

The crucial correction was the dominance bound.  Pure w-free residue emptiness
is too strong: a residue class such as `0 mod M` can survive formally at huge
multipliers.  But if the first positive representative is above

```text
floor((n-1) max(S') / n),
```

then dominance already proves looseness.  This converts ghost alignments into
certificates.

## Creative upgrades

1. **Prime-power ZDD.**  Do not build residues modulo the full lcm.  Store each
   component's allowed residues over prime powers, then combine by a
   zero-suppressed decision diagram.  Most local choices are absent, so a ZDD
   should fit the proof object.
2. **Minimal contradiction extraction.**  When CRT intersection is empty or
   bounded-empty, extract a small incompatible set of components.  That set is a
   human-readable proof certificate.
3. **Density certificate.**  Instead of exact residues, bound the allowed density
   per component and prove that no allowed representative can occur before the
   dominance cutoff.  This is the analytic cousin of the automaton.
4. **Ghost-class theorem.**  Formalize why surviving classes first appear beyond
   the dominance cutoff in the sampled rows.  That may be the large-owner lemma:
   all feasible large-owner alignments require a multiplier too large to remain
   in the hard Cprime band.
5. **Paper-seed compiler.**  Translate `n=11,12,13` paper seed families into
   component residue clauses.  If a seed fails across primes, ask which `w`
   residue language it was trying to realize.
6. **Meet-in-the-middle CRT.**  Split components into two halves, summarize each
   half by residues and modulus factors, then join summaries.  This should make
   the current overflow cases inspectable.
7. **Owner-size descent.**  Large owner `u>=n` is itself a runner in `S'`.
   Attach its own endpoint-owner obligations and see whether high slack forces
   a private pivot or a pair-pinch witness.

## What to ignore

Do not return directly to `p^k` tuples, arbitrary ansatz grids, or runner
tournaments.  They remember the wrong thing.  The residual is now a
single-multiplier feasibility language.  Full exact fibers should appear only
if the residue automaton has live bounded states after factored CRT and
dominance pruning.

## Tournament Analysis

The vertices are proof gates or automaton states, not runners:

```text
Bprime_long_interval
LemmaC_small_owner
bounded_CRT_empty
live_CRT_state
modulus_overflow
exact_fallback
```

Observable:

```text
(unproved_fraction, modulus_growth, dependency_depth, -proof_strength, tie_order)
```

Switch toward harder residual burden.  The current fingerprint is transitive
with one Hamiltonian path.  That is good: it means the algorithmic queue is
ordered.  Nontrivial cycles should only arise after splitting `modulus_overflow`
by prime-power factors.

## Handoff

The next concrete script should be `lrc_large_owner_factored_crt_s582.py`:

```text
component -> allowed local residues modulo prime powers
local residues -> ZDD/BDD
ZDD -> bounded emptiness certificate or minimal live state
```

The success condition is not "more rows checked."  It is a short certificate
for each overflow row: either a minimal incompatible component set, or a live
bounded residue state with all owner pairs displayed.
