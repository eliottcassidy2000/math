---
source: opus-2026-06-03-S574 (remote-control)
status: NEW THEOREM (Lemma C, w-free) via the cover→congruence translator; two proved criteria cover 99% of multiple-of-14 configs; residual = large-binding-owner ~1%
tags: [LRC, translator, endpoint-owner, congruence, small-owner, THM-398, C-prime, n14, cover-assignment]
---

# The cover→congruence translator, and the small-owner theorem it hands you

**Prompt (user):** focus not on more prime enumeration, but on a translator from
all-short interval-cover assignments to endpoint-owner congruences.

Built exactly that translator. Its first output is a **new w-free theorem** that —
together with the Vitali criterion B′ — discharges **99% of the multiple-of-14 configs**
without any enumeration.

## 1. The translator

For `S = S' ∪ {v=nw}` (`n|v`), tightness is `G(S') ⊆ D_v`: every component `(a,b)` of
the other runners' safe set sits inside one danger arc of `v` (centre `j/(nw)`,
radius `1/(n²w)`). Each endpoint has an **owner** — the runner that turns (un)safe
there:
```
left  a = (k_a n + 1)/(n u_a)        right  b = (k_b n − 1)/(n u_b).
```
"`(a,b)` ⊆ arc `j`" translates (clear denominators by `n u w`) into the
**endpoint-owner congruences**
```
|w(k_a n + 1) − j u_a| < u_a/n,      |w(k_b n − 1) − j u_b| < u_b/n.
```
**Verified:** this congruence-window test agrees with direct tight/loose computation
**100% (2500/2500 each n=6..14)** — the translation is exact, not heuristic.

## 2. The rigidity that makes it a theorem (Lemma C)

The right sides are `u/n`. For an owner `u < n` they are `< 1`, and the bracket is an
integer — so it is **forced to `0`**: the endpoint must equal the arc centre. That is
total rigidity: a small owner pins its endpoint *exactly* on `j/(nw)`. Hence

> **Lemma C.** If a component of `G(S')` has **both** owners `< n`, then `S = S'∪{nw}`
> is loose for **every** `w`. *(Two rigid endpoints would have to be the same arc
> centre, forcing `a=b`.)*

The "cross-relation" two small owners must satisfy, `u_b(k_a n+1) = u_a(k_b n−1)`, is
**identically `a=b`** — verified `0` holds out of `490/1330/2594/5298` small-owner
components (`n=6,8,10,12`). It can never hold; the small-owner congruences are
**self-contradictory**. This is the answer to the user's ask: the translated
congruences are not just bookkeeping — in the small-owner case they are *infeasible by
inspection*.

Slack enters the congruence (`|·| < u/n` allowing a nonzero integer) **only** for
owners `u ≥ n`. So an endpoint can sit off-centre inside an arc *only* if its owner is
large — pinpointing where the difficulty can possibly hide.

## 3. Coverage — two proved criteria, 99% at the frontier

`lrc_smallowner_coverage_s574.py`, multiple-of-`n` configs:

| n | Lemma C alone | Lemma C ∪ B′ (Vitali long-interval) | residual |
|---|---|---|---|
| 6 | 8.2% | 73.4% | 531 |
| 8 | 18.8% | 81.7% | 366 |
| 10 | 33.4% | 92.0% | 160 |
| 12 | 56.2% | 96.0% | 80 |
| **14** | **81.3%** | **99.0%** | **19/2000** |

Both proved criteria **strengthen toward the frontier** (Lemma C alone: 8%→81%). At
`n=14`, **99%** of multiple-of-14 configs are now loose *by proof*, no enumeration.

## 4. The residual, exactly

What's left (~1% at `n=14`): configs where **every** component of `G(S')` is short
**and** **every** component has at least one **large** binding owner (`≥ n`). Only
there does the congruence slack (`±u/n`) let an endpoint fit off-centre. So the open
core of C′ is now a **large-owner congruence** question:
```
for each short component, ∃ integer j with  w(k_a n+1) ≡ j u_a  (within u_a/n)
                                       and   w(k_b n−1) ≡ j u_b  (within u_b/n),
   where at least one of u_a,u_b ≥ n.
```
This is a genuine CRT/Diophantine feasibility on the *large* owners and the multiplier
`w` — finite per config, and (verified) never simultaneously satisfiable, but not yet
proved. The reduction has moved from "all binding owners" to "the large binding owners
only," a strictly smaller and more structured target.

## 5. Why this is the right frame (no enumeration)

The user's redirection was decisive. Enumerating speed sets / primes only ever yields
"0 counterexamples in a box." The translator instead converts the *geometry* (does a
periodic arc family cover a union of intervals) into *arithmetic* (a congruence system
on the endpoint owners), where **infeasibility is provable**: for small owners it is
immediate (`a=b`), and for large owners it is a bounded CRT check. Two of the three
regimes (long interval; small owners) are now closed; only large-owner short
components remain.

## 6. Honest status

- Translator: exact, **verified 100%**.
- **Lemma C** (both-small-owner ⟹ loose, every `w`): **PROVED**, folded into THM-398.
- Lemma C ∪ B′: **PROVED**, covers **99%** of multiple-of-14 configs (verified).
- Residual (large-owner short components, ~1%): **OPEN**, now a bounded congruence
  feasibility, not an enumeration.

**Artifacts:** `04-computation/lrc_cover_to_congruence_translator_s574.py` (+`.out`),
`lrc_smallowner_coverage_s574.py` (+`.out`), `lrc_cross_relation_audit_s574.out`.
Builds on THM-398/HYP-2104 (the reduction + Vitali criterion), HYP-2103 (dodge),
HYP-2102 (C′). New: **HYP-2105**.
