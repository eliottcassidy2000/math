---
source: codex-2026-06-03-S583
status: SUPPORTED proof-program refinement; Lemma A and Lemma B still open
tags: [LRC, fold-sieve, circuit-free, additive-relations, Phi, Tournament-Analysis, n14]
---

# Fold/sieve and randomness: 3-term folds are the boundary

**Prompt (user):** focus mostly on Lemma B, but go back and forth between
Lemma A and Lemma B many times, using one as noise for the other.

The useful split is now cleaner:

```text
Lemma A side: no low-rank circuit, especially no a+b=c
Lemma B side: a+b=c, so c t = a t + b t is a literal fold
```

I treated "no 3-term relation" as a circuit-free proxy and used it as a
control population while probing the fold side.  The control did its job.  It
made the fold signal stand out, and it made 4-term additive energy look less
like the right boundary.

## The fold is a real gate

For the local model, compare independent clearance

```text
x safe, y safe, z safe
```

with the folded clearance

```text
x safe, y safe, x+y safe.
```

The folded gate has a persistent density penalty from `n=6` to `n=14`; at
`n=14` it is still about `-0.01286` relative to the independent model.  This is
small, but it is not cosmetic.  A relation `c=a+b` really removes an independent
coordinate and replaces it by the diagonal equation `c t = a t + b t`.

That supports the user's "fold+sieve" language exactly: the fold is geometric;
the sieve still has to find the clock or prime-divisibility route that escapes
the folded shield.

## The 4-term control changed the target

The decisive stress test was AP versus shifted AP:

```text
AP(k)          = (1,2,...,k)
shifted_AP(k) = (k+2,k+3,...,2k+1)
```

They have identical ordered additive energy.  But AP has many 3-term folds and
is tight, while shifted AP has no `a+b=c` relation and has a large positive
margin.

At `n=14`, both have energy `1469`:

```text
AP:          M=1/14, 36 folds, safe measure 0
shifted AP: M=5/14, 0 folds, safe measure about 0.117
V*:          M=1/14, 31 folds, safe measure 0
```

So "4-term-rich" is not the hard row type by itself.  It becomes dangerous only
when it carries lower-rank fold structure.  This retires the old worry in a
more pointed form: 4-term structure can be noise, but 3-term structure is an
actual mechanism.

## What the deletion shadow says

For rows with a relation `a+b=c`, deleting `c` often exposes a good witness
that gets blocked when `c` is reattached.  In S583 this happened in `446/934`
relation rows.  At the same time, the exact witness denominator never equalled
the folded shield `c=a+b`.

That combination is a useful warning.  The naive pair-sum clock is often the
thing the folded runner is built to block.  A proof of Lemma B should not expect
the shielded denominator to win directly.  It should use the fold to define the
bad diagonal, then route through another modulus, owner pin, or prime-fibre
certificate.

This is where the paper-sieve machinery and HYP-2106 still matter: not as a
giant terminal enumeration, but as a way to route around a fold once the local
geometry has identified the shield.

## How it plugs into `Phi`

HYP-2112 says the multiple-of-`n` branch has an exact value functional

```text
G(v) = Phi(C)
```

and the worry set is `ker Phi`.  S583 does not prove `Phi>0`; it suggests which
routes should try to prove it:

```text
3-term-free row  -> discrepancy/equidistribution lower bound
3-term fold row  -> folded diagonal + sieve/owner/residue escape
```

This is compatible with HYP-2108 and HYP-2110.  Endpoint-cover positivity and
small-owner descent prove one component pokes out of the cover; `Phi` records
the exact uncovered measure; the fold/randomness split tells us which proof
language to use before we reach the component ledger.

## Tournament Analysis

Vertices were proof-state buckets, not runners:

```text
A_circuit_free_proxy
A_high_4term_no3
B_has_3term_fold
B_fold_shadow_blocked
```

The observable was `(cleanliness rank, coverage)` and the switch favored easier
proof states.  The tournament was transitive:

```text
score_hist={0:1, 1:1, 2:1, 3:1}
directed_3_cycles=0
path=A_circuit_free_proxy -> A_high_4term_no3
     -> B_has_3term_fold -> B_fold_shadow_blocked
```

This proof-bucket tournament is intentionally boring.  The next nontrivial
Tournament Analysis should use fold gates, residue states, endpoint-cover
components, or prime-fibre obligations as vertices.  The bucket version is a
route ledger; it should not be mistaken for the proof object.

## Assumption challenge

Candidate vertices considered: runners, 3-term relations, fold gates, additive
energy fibres, exact witnesses, endpoint-cover components, and proof buckets.
The chosen quotient preserves exact threshold clearance and plausible proof
route, but it destroys endpoint-owner circuits and prime-fibre residue
certificates.

The challenged assumption is now explicit: additive energy is not the right
hardness coordinate unless it descends to a low-rank circuit.  The AP and
shifted-AP control is the clean witness.

## Next proof targets

1. Prove a discrepancy lower bound for rows with no 3-term additive circuit.
2. Prove a folded gate lemma for `a+b=c`, phrased on `(a t,b t)` with
   `c t=a t+b t`.
3. Attach the fold gate to the existing endpoint-cover and small-owner
   machinery, so one route produces a positive `Phi` summand.
4. Replace the shielded pair-sum denominator by a routed modulus or
   prime-fibre certificate.

## Incoming convergence

After rebasing, two pieces of shared work landed that sharpen this reflection.
The S577 additive-circuit A/B probe verifies the fold-as-shield identity in its
most literal form: when `v_c=v_a+v_b`, runner `c` is at integer phase at the
`1/v_c` pinch clock.  Its circuit-free scan also keeps margins safely above
`delta=1/(k+1)` through `k=10`, matching the direction of S583's no3 buckets.

The follow-up S577 round-2/3 probes push this closer to the proof split:
near-tight rows are 3-term-rich, circuit-free rows vanish from the near-tight
bucket, and the remaining near-tight rows either have a multiple of `k+1` or
are discharged by the `j/(k+1)` delta-clock witness.  That is exactly the
handoff pattern HYP-2114 wants: folds explain local shielding, while Cprime and
`Phi` own the multiple branch.

Oracle S578o then supplied the missing conceptual name: 4-term energy is the
translation-invariant shadow of folds.  The S583 exact-maximin data is therefore
best read as support for that shared HYP-2114, with `Phi` and the deletion-shadow
diagnostics added to the proof route.

The new HYP-2113 tournament speedup stack is the right algorithmic host for this
lemma split.  Fold/shield exits should be labelled in the certificate DAG before
raw exact interval search or CRT is invoked.

**Artifacts:** `04-computation/lrc_fold_sieve_random_bridge_s583.py`
(+ `05-knowledge/results/lrc_fold_sieve_random_bridge_s583.out`).  Supports
shared **HYP-2114**.
