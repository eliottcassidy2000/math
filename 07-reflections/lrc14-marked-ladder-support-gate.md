# LRC14 marked ladder support gate

**codex-2026-06-12.**  A support-gate push on HYP-2438/THM-492, using the
recent marked-support lessons from the order-5 `[72,36,16]` code branch and
the Erdos-Moser tower branch.

## What moved

The old scalar question was:

```text
does some denominator q give a lonely tick?
```

The better question is:

```text
if q fails, which runners paid to block all unit twists?
```

For a fixed denominator `q`, every unit twist `a/q` is blocked by at least one
runner.  That is a hitting-set problem on the runners.  The minimum cover
`tau_q(S)` and the cross-denominator load of the covering runners are the
missing address coordinates.

This is the same shape as the new support-gate motif:

```text
72-code: mark fixed coordinates/tetrads, not just the theta coefficient.
Erdos-Moser: mark support packets, not just the 2a+1 scalar recurrence.
LRC14: mark blocker runners, not just the fact that q is blocked.
```

## The useful separation

The single-stranger rows

```text
S(r) = 7*{1,...,12} union {r}
```

were mentally doing two jobs at once.  The marked ladder atlas separates them.

They are failures of the **pure shell** band-1 horizon `q<=27`, with the clean
signature `13|r` and `r mod 27 in {0,+/-10}`.  But they are not failures of the
fibered ladder in the same sense: the address `q=91=7*13` already sees the
7-core as the proven `n=13` problem, while `q=40/41` is the non-fibered rung-up
rescue.

That is a real conceptual gain.  The phrase "band-1 failure" must say which
ladder it means: pure shell or fibered shell.

## The proof angle

For a denominator `q`, there are three qualitatively different ways to block:

1. a universal blocker, usually a divisibility/apex signal;
2. a small shared cover, where a few runners divide the unit twists;
3. a diffuse high-cost cover, where `tau_q` is large.

The known wall atoms AP, Vstar, and 2AP are diffuse all the way up: they block
every tested fibered denominator and reach `tau_q=13`.  Those should be
quarantined as normalized floor atoms.

For primitive loose targets, the pressure should become productive.  If
universal blockers dominate, go to apex debt and HYP-2256.  If cover-blocking
dominates without universal blockers, repeated canonical cover load names a
runner that should be deletable or Bprime-dodgeable.  The proof would become a
pincer:

```text
finite ladder witness
or marked blocker load creates an opening certificate.
```

## Concrete next work

1. Prove the pure-shell/fibered-shell separation for the S(r) family by hand:
   the `13|r`, `r mod 27 in {0,+/-10}` condition blocks pure shell `q<=27`;
   `q=91` or `q=40/41` supplies the escape.

2. Replace the script's canonical-cover heuristic by a canonical invariant:
   for each `q`, compute the full family of minimum blocker covers up to the
   row automorphism/stabilizer.  The load should be a set, not a tie-broken
   artifact.

3. Attach blocker type to HYP-2256's endpoint owners.  A universal blocker at
   many `q` should be an apex-debt row; a repeated non-universal blocker should
   be the owner-private deletion candidate.

4. Search not for counterexamples but for maximal support-pressure rows after
   excluding AP/Vstar/2AP.  Those rows should be the testbed for the repeated
   cover-load lemma.

## Honest scope

This does not prove LRC(14).  It improves the target: HYP-2438 should no
longer be attacked as a raw finite denominator search.  It should be attacked
as a support-realization theorem for the blocker hypergraph.
