---
source: codex-2026-06-02-S563
status: reflection + exact denominator-tier audit
tags: [LRC, sieve, lower-bound, recursive-sieve, Tao, Bedert, Tournament-Analysis]
---

# Recursive multi-sieve lower-bound ledger

The user asked for another long push toward improving the LRC lower bound
beyond Tao, now through recursive multi-sieves.  The honest answer is still:
not as a universal asymptotic theorem.  Tao's global surplus and Bedert's later
Riesz-product surplus are the right external benchmarks.

But the repo's case-wise route is not empty.  It says: stop asking for one
global inequality first.  Separate the speed sets into obligations.  The coarse
THM-369 sieve already gives the full conjectural bound `1/(k+1)` whenever some
`q <= k+1` divides no speed.  That removes the bulk.  The real problem is the
sieve-covered residual core.

S563 asks whether recursive denominator tiers can remove more of that core.

## The ledger

The tiers are deliberately simple:

```text
T0: q <= n
T1: n < q <= 2n
T2: prime powers with base prime <= n
T3: n-smooth denominators
T4: remaining q <= Q
```

For each deterministic speed set, the script scans closed rational witnesses
`t=a/q` and records the first tier that reaches `1/n` plus the best witness up
to `Q=220`.

The useful samples were chosen to be annoying:

```text
blind_lcm_2_14 = {1,...,12,lcm(2,...,14)}
blind_lcm_2_24 = {1,...,12,lcm(2,...,24)}
```

Both defeat every `q <= 14` denominator witness.  Still:

```text
t=2/27 gives margin 2/27 > 1/14.
```

That is the fine-window phenomenon in its most blunt form.  Cover all the
small denominators you like; the first spreading window can still open a full
lonely corridor.

The S562 residual packets also clear:

```text
n=14 packet:      t=6/23, margin=2/23
n=14 dyadic lift: t=3/23, margin=2/23
n=17 packet:      t=9/25, margin=2/25
n=18 packet:      t=2/29, margin=2/29
```

This links the two previous pictures.  HYP-2073 saw dyadic packet lifts export
visible gap into boundary debt.  S563 says some of those exported packets are
not hard yet; a nearby recursive rational denominator certifies the full LRC
bound.

## What this does and does not prove

It does not beat Tao or Bedert globally.  It gives a strategy for shrinking the
place where global asymptotic tools are needed.

The proof target becomes:

```text
Every residual obligation either has a local tier witness >= 1/n,
or exports positive product-frontier mass to child obligations.
```

That is a much more structured object than "all k-speed sets."  It is also the
right level for the cascade-product language: local clearance factors either
multiply to a certificate, or the failed factor is carried as explicit debt to
the next tier.

## Tournament Analysis

For this pass, vertices are speed-set proof obligations.  The observable is
the lower-bound certificate:

```text
best margin, first tier reaching 1/n, witness denominator.
```

The tie path is the displayed sample order.  The resulting fingerprint is
transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3_cycles=0
sccs=[1,1,1,1,1,1,1,1]
hamiltonian_paths=1
```

This is not surprising.  A lower-bound ledger orders obligations by certificate
strength.  It forgets exactly the data that could produce cycles: endpoint
owners, wall-crossing events, local cover-arc overlaps, and signed cross-prime
couplings.

The hidden transitivity fact the user emphasized should probably enter at that
next level.  If `(X,Y)` and `(Y,Z)` force `(X,Z)`, then the no-return condition
forbids the resonance triangle `Z -> X -> Y -> Z`.  In LRC terms, the next
Tournament Analysis quotient should remember which endpoint debt is being
carried forward, not just that a margin exists.

## Assumption challenge

I considered denominators, residues, CRT channels, residual packets, endpoint
owners, runners, gaps, circle sections, wall crossings, Fourier modes, matroid
circuits, and proof obligations as vertices.

I chose proof obligations for S563 because that quotient preserves the
predicate under test:

```text
does some recursive tier certify margin >= 1/(k+1)?
```

It destroys ownership and overlap data, so it cannot prove the no-return/cycle
mechanism by itself.  That loss is now a feature of the handoff: the next
script should refine each obligation into owner-labelled residual packets.

## Handoff

Build the next ledger as:

```text
ResidualObligation(
    speed_set,
    tier,
    witness_or_none,
    endpoint_owner_multiset,
    product_frontier_mass,
    child_obligations,
)
```

Then test the descent lemma:

```text
local witness
or positive frontier mass
or conserved export to a strictly lower/child obligation.
```

That would be a real path toward an improvement beyond Tao's operative regime:
not a bigger global surplus, but a proof that the Tao/Bedert hard core keeps
collapsing under recursive sieving.

**Artifacts:** `04-computation/lrc_recursive_multisieve_lower_bound_s563.py`,
`05-knowledge/results/lrc_recursive_multisieve_lower_bound_s563.out`,
HYP-2076.
