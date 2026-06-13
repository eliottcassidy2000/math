---
source: opus-2026-06-08-S709 (CORRECTS S708 / MISTAKE-064; user caught the misread)
status: CORRECTION + honest treatment. Erdős Problem 64 = the Erdős–Gyárfás conjecture: min degree ≥3
  ⟹ a cycle of length a POWER OF TWO (2^k: 4,8,16,…) — OPEN, falsifiable. I had misread "2^k" as "2k"
  in S708 and proved/mislabeled the trivial even-cycle statement (THM-443, settled-true) as the open
  problem. Here: the real problem (THM-445), known results (cubic-planar Heckman–Krakovski; Bondy–Vince
  spectrum; Markström searches), computational verification (atlas n≤7; Petersen/Heawood/Pappus… girth
  ≥5 all caught by an 8- or 16-cycle; 48 random cubic n≤16 — 0 counterexamples), and the DYADIC framing:
  even (settled) vs power-of-two (open) = additive parity mod 2 vs the multiplicative 2-adic tower.
tags: [erdos-problem-64, erdos-gyarfas, power-of-two-cycle, dyadic, 2-adic, open-problem, correction,
  mistake-064, even-cycle, bondy-vince, heckman-krakovski, markstrom, cubic-graph, girth, cut-cycle]
---

# Erdős Problem 64: power-of-two cycles, and the dyadic target I missed

**Prompt (user, correcting me):** "the problem i asked you to work on was Erdős problem 64, Does every
finite graph with minimum degree at least 3 contain a cycle of length 2^k for some k≥2? it is
currently open and falsifiable, yet unproven."

The user is right and I was wrong. In S708 I read **`2^k` as `2k`**, proved the classical *even*-cycle
statement (min degree ≥3 ⟹ even cycle ≥4 — true, settled), and mislabeled it as the open problem.
Logged as **MISTAKE-064**. The tell I ignored: *a settled-true statement can't be an open falsifiable
problem.* This reflection records the real problem and does it honestly.

## The real problem (and how much harder it is)

> **Erdős–Gyárfás (1995):** every graph with `δ≥3` contains a cycle of length a **power of two**
> (`{4,8,16,32,…}`). OPEN.

The even-cycle lemma says the cycle spectrum meets `2ℤ` (every even number is a candidate). Erdős 64
asks whether it meets the **dyadic ladder** `{2^k}` — an infinitely thinner target. A counterexample
must be a min-degree-≥3 graph dodging `4, 8, 16, 32, …` **simultaneously** — high girth, sparse cycle
spectrum. That is why "even" is a one-line exercise and "power of two" has resisted for 30 years.

## What is known / what I verified

- **Proved for cubic planar graphs** (Heckman–Krakovski 2013); **no counterexample** in computational
  searches (Markström et al.). Min degree ≥3 forces a rich spectrum — **Bondy–Vince:** two cycles with
  lengths differing by 1 or 2 — but "rich" need not "hit a power of two."
- **Verified (S709):** all 173 `δ≥3` graphs on `≤7` vertices have a power-of-two cycle (a 4-cycle,
  forced since girth `≤4` there). In the substantive **girth-≥5** regime (no 4-cycle): Petersen
  (`{5,6,8,9}`), Heawood, Pappus, Desargues, Möbius–Kantor, Dodecahedral — **all** carry an 8- (and
  often 16-) cycle. 48 random cubic graphs `n≤16`: 0 counterexamples. Consistent with the conjecture;
  not a proof.

## The dyadic insight (the honest repo connection)

> **Between "even" and "power of two" is exactly the gap between additive parity (`mod 2`) and the
> multiplicative 2-adic tower (`{2^k}`).** THM-443 (even) is the additive face; Erdős 64 is the
> *multiplicative/dyadic* face — the same additive↔multiplicative split as THM-442 (the IE/`Δ³`
> cell-recursion vs the multiplicative-modular H-law). Power-of-two is the **pure 2-adic tower**,
> resonant with the repo's 2-adic depth (S701/S704 prime-power towers; the order-2 antipodal `σ`,
> S702) and the **Cayley–Dickson dimension tower `2^k`** (CLAUDE.md). This is a *lens/resonance*, not
> a reduction: Erdős 64 is a cycle-spectrum problem, not an LRC/tournament one. The S707 Pfaffian/
> even-dicycle bridge applies to **even** cycles (Pólya's `per=±det`), **not** to power-of-two — I
> must not conflate them (that conflation was the root of MISTAKE-064).

## What a real attack would need (honest next steps)

- The conjecture is about forcing a dyadic length. A promising angle (literature): bound the **girth**
  and the **cycle-length spectrum spread** so that the dyadic gaps (`[2^k, 2^{k+1}]` doubling) must be
  crossed — Bondy–Vince gives spread-1-or-2, but the dyadic gaps grow, so a single near-pair does not
  suffice; one needs many lengths in one dyadic window.
- Cubic case first (where it's proved planar): extend beyond planar via the structure of cubic graphs.
- Computational: push the no-counterexample search to higher girth / order (Markström's method).

## Honest status

- **OPEN.** Nothing here resolves it. The verification is consistent with the conjecture on all tested
  classes (atlas n≤7, girth-≥5 named graphs, random cubic n≤16) — no counterexample.
- **Correction:** supersedes the S708 framing; THM-443 rescoped to the even-cycle lemma; the real
  problem is **THM-445**; **MISTAKE-064** logged. HYP-2313's parity-covering lens stands, but its
  even-cycle leg is the *weaker* (settled) statement, not Erdős 64.
- **Lesson kept:** parse a stated open/falsifiable claim as a constraint on interpretation — if your
  proof settles it, you've misread it.

**Artifacts:** `04-computation/erdos_gyarfas_power_of_two_cycle_s709.py` (+`.out`). Theorem record
**THM-445**. Corrects THM-443, MISTAKE-064. Cites Erdős–Gyárfás, Heckman–Krakovski, Bondy–Vince,
Markström.
