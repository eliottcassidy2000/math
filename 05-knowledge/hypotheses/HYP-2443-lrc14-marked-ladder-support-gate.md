# HYP-2443 - LRC14 marked ladder support closes the next proof gap

**Status:** OPEN proof-route refinement.  
**Source:** codex-2026-06-12.  
**Companions:** HYP-2438, THM-492, HYP-2256, HYP-2253, HYP-2241,
HYP-2167, HYP-2164.  
**Computation:** `04-computation/lrc14_marked_ladder_support_codex.py`;
stored output `05-knowledge/results/lrc14_marked_ladder_support_codex.out`.

## Claim

For LRC total runner count `n=14`, the HYP-2438 band-ladder route should be
proved through a marked support ledger.

For a speed set `S`, denominator `q`, and unit twist `a`, say runner `v`
blocks `a/q` when

```text
||a v / q|| <= 1/14.
```

If every unit twist at `q` is blocked, form the blocker hypergraph whose
vertices are runners and whose hyperedges are the sets of runners blocking each
unit twist.  Let

```text
tau_q(S) = minimum number of runners hitting all unit twists at q.
```

The sharpened proof route is:

1. finite ladder: find a witness in the fibered ladder
   `Q_K={d*m : d|14, m<=K}`;
2. support pressure: if the ladder is blocked, use the marked data
   `(tau_q, canonical blocker load, universal blockers)` to force either
   Bprime(any runner) or an owner-private deletion/apex-opening certificate.

In short, raw `q`-blocking is a scalar shadow.  The proof object is the marked
runner support that pays for blocking all twists.

## Evidence

The atlas audits `Q_27` and `Q_41`, where

```text
Q_K = {d*m : d in {1,2,7,14}, 1 <= m <= K}.
```

The known floor atoms remain walls:

```text
AP, Vstar, 2AP:
  no Q_27 or Q_41 witness,
  every denominator blocked,
  many high-cost cover-blocked q,
  tau_q reaches 13.
```

This is expected: they are the normalized wall atoms, not rows to prove loose
by the ladder alone.

The single-stranger family

```text
S(r) = 7*{1,...,12} union {r}
```

separates two notions that were easy to conflate:

- pure shell denominators `q <= 27` are blocked in exactly the sampled
  evader pattern `13|r` and `r mod 27 in {0,+/-10}`;
- the fibered ladder already catches many of these through addresses such as
  `q=91=7*13`, while the non-fibered rung-up rescue appears at `q=40` or `41`.

In the scan `r<1200`, pure-shell band-1 failures caught by `q<=41` are:

```text
r = 260,351,442,611,702,793,962,1053,1144.
```

The random structured search retained `13/600` high-pressure rows.  The
interesting rows are not those with many universal blockers; they are rows
with many cover-blocked denominators, where two or more marked runners share
the blocking cost.  That is exactly where owner-private deletion and
Bprime(any runner) have content.

## Interpretation

This is the LRC14 version of the recent support-gate pattern:

```text
72-code: scalar theta/design ledger is not enough; mark the length-16 support.
Erdos-Moser: scalar 2a+1 recurrence is not enough; mark the packet support.
LRC14: scalar blocked denominator is not enough; mark the runner support.
```

The missing coordinate is a set-cover/address coordinate.  It records whether
a denominator is blocked by one universal/apex runner, by a small shared cover,
or by a high-load runner that keeps reappearing across denominators.

## Next Lemmas

1. **Pure shell versus fibered ladder lemma.**  Formalize that the S(r) evaders
   are pure-shell band-1 failures but not fibered-ladder failures once `q=91`
   is admitted.
2. **Universal blocker lemma.**  If many blocked denominators are universal,
   the row has visible divisibility/apex debt and should reduce to HYP-2256.
3. **Repeated cover-load lemma.**  If no finite witness exists through a fixed
   rung and universal blockers do not dominate, some runner must appear in many
   canonical covers; this runner is the marked candidate for Bprime(any runner)
   or owner-private deletion.
4. **Wall exclusion.**  Separate AP, Vstar, and 2AP as normalized floor atoms
   before applying the support-pressure lemma to primitive loose targets.

## Tournament Analysis

Vertices are proof routes/support quotients, not runners.  Candidate vertex
sets considered include runners, denominators, unit twists, danger bands,
blocker hitting sets, apex speeds, stranger speeds, safe intervals, and proof
obligations.

The selected route tournament is transitive:

```text
marked_ladder_setcover
> owner_private_deletion
> apex_opening_modes
> Bprime_any_runner
> raw_Q_search
> raw_scalar_floor
```

The quotient preserves the proof-relevant predicate "how expensive is it to
block this rung?" while deliberately forgetting most raw time geometry.
