---
source: opus-2026-06-02-S556 (remote-control)
status: rigorous proof-lite (impossibility for structured counterexamples) + a brainstorm of provable LRC statements with evidence
tags: [LRC, n14, n7, impossibility, proof-lite, even-fold, antipodal, sieve, tight-witness-lattice, S554, S553, S552o]
---

# LRC@14: a proof-lite impossibility, and adjacent statements we can prove

**Prompt (user):** attempt to prove it is impossible to create a counterexample
"like this" (a proof-lite); think of other LRC statements to make progress on,
using the repo's unique insights.

Notation: n=14, observer + 13 speeds; `||x||` = distance to nearest integer;
`M(S)=max_t min_i ||v_i t||`; LRC(14) ⟺ `M(S) ≥ 1/14` for all S; a *counterexample*
has `M(S) < 1/14`. `fold(S) = {v/2 : v∈S even}`.

## Part I — Two rigorous lemmas (the repo's fresh tools, cleanly stated)

**Lemma 1 (even-fold).** For an even speed `v=2u`, `||v t|| = ||u·(2t)||`. Hence
with `s=2t`, `min_{v even}||v t|| = g_fold(2t)`, and
`M(S) = max_t min(g_fold(2t), g_odd(t)) ≤ max_s g_fold(s) = M(fold(S))`. ∎

**Lemma 2 (antipodal).** For an *odd* speed `v`,
`||v(t+½)|| = ||v t + v/2|| = ||v t + ½||` (since `v` odd ⇒ `v/2 ≡ ½`). The two
values `||vt||` and `||vt+½||` cannot both be `< 1/14` (their arguments are `½`
apart, `> 2/14`). So an odd runner is dangerous at **at most one** of any
antipodal pair `{t, t+½}`. ∎  *Corollary:* any all-odd set is lonely at `t=½`.

## Part II — The proof-lite: structured counterexamples are impossible

**Theorem A (q=14 sieve; THM-369).** If `14 ∤ v_i` for all `i`, then `t=1/14` is a
witness (`||v_i/14|| = min(r,14−r)/14 ≥ 1/14`, `r=v_i mod 14 ≠ 0`). ∎

> **Corollary A′ (self-contained, the "like this" families).** *Every config
> obtained from the AP `{1,…,13}` by doubling any subset of its even runners is
> lonely at `t=1/14` — hence is not a counterexample.* Proof: doubling `2m`
> (`1≤m≤6`) gives `4m`; `14 | 4m ⟺ 7 | 2m ⟺ 7 | m`, impossible for `m≤6`. All
> other runners are original AP members in `{1,…,13}`, none `≡0 (mod 14)`. So no
> speed is a multiple of 14 and Theorem A applies. ∎ (Verified: all 23 distinct
> such configs lonely at 1/14; the only tight one is `V*={1..11,13,24}`.)

This rigorously rules out the entire S555 Part-A search space (and, via Theorem A,
every probe there lacking a multiple of 14). The honest content: the constructive
disproof routes fail because **a counterexample is forced to contain a multiple of
14** — the obstruction the sieve names.

**Theorem B (necessary conditions for any counterexample, rigorous).**
Any LRC(14) counterexample `S` satisfies:
1. `S` contains a **multiple of 14** (else Theorem A gives a witness);
2. that multiple lies in the **mod-7 singleton class** (multiples of 7) *and* is
   even — i.e. it is the meeting point of the mod-2 fold (Lemma 1) and oracle-
   S552o's mod-7 CRT singleton; the two factorisations of `14=2·7` obstruct on the
   *same* runner;
3. `S` has **≥ 7 speeds of one parity** (since `#even+#odd = 13` is odd, the
   minority parity has `≤ 6`, the majority `≥ 7`). ∎

These are modest but exact, and (2) is a genuinely new bridge between the two
prime-factor reductions developed this week.

## Part III — Adjacent statements, with evidence (the brainstorm)

Ranked by how close they are to a proof:

**C. Tight-witness lattice (CONJECTURE, strong evidence).** *Every tight n=14
config is lonely only at `t ∈ (1/14)ℤ`.* Within distance ≤2 of the AP (adds ≤46)
the only tight configs are **AP and V***, both witnessed exactly at
`{1,3,5,9,11,13}/14` (`lrc_n14_tight_lattice_B3_s556.out`); the n=5,6,8 sporadics
are likewise witnessed at `j/n` (S553). **Cleaner sufficient form, well-supported:**
*tight ⟹ no multiple of 14* (⟹ witnessed at 1/14 by Theorem A). Evidence: both
known tight configs have no multiple of 14, and a directed hunt of **2500
hill-climbs forced to contain a multiple of 14 reached tightness 0 times** — such
configs are robustly *loose* (`lrc_n14_tight_no_mult14_test_s556.out`).

> **A productive tension (heuristic for why LRC@14 holds).** Theorem B(1): a
> counterexample *must* contain a multiple of 14. Conjecture C′ (strongly
> supported): a multiple of 14 keeps a config *loose* (away from tight, let alone
> sub-tight). The two requirements pull in opposite directions — the very feature
> a counterexample is forced to have is the one that prevents non-loneliness.
> Turning this tension into a quantitative lower bound (`M(S) ≥ 1/14` whenever
> `14 | v` for some `v`) is a concrete, plausibly provable target.

**D. e≤6 fold reduction (CONJECTURE, 127/127).** By Lemma 1, for `e=#even ≤ 6` the
even danger occupies `≤ 6·(1/7) < 1` so evens are safe on a positive-measure set;
by Lemma 2 each odd is dangerous at ≤1 preimage of a doubled time. So
`LRC(14)[e≤6] ⟸ some even-good s has no odd-split` (S554). Open: prove no-odd-split.

**E. Spectral-gap completeness (OPEN, with my S553 correction).** Oracle's gap
reduction to antipodal transversals is *incomplete* when `2n−1` is composite
(`27=3³` at n=14): the non-unit-pair non-transversals (e.g. V*) escape it. A full
proof of the gap must cover **transversals ∪ non-unit-pair non-transversals**.

**F. Even-fold across the proven range (PROVABLE COROLLARY).** For every even
`n=2p` with `p ≤ 7`, Lemma 1 + LRC(p) give `M(fold) ≥ 1/p` whenever `#even ≤ p−1`,
so the even half of such configs is settled by a *proven* case. (n=14 uses
LRC(7); n=10 uses LRC(5); etc.) This is the general form of the "7 lever."

## Part IV — Honest status

The **proof-lite is real**: Theorems A, A′, B are rigorous and rule out explicit
infinite/finite families of would-be counterexamples (everything lacking a
multiple of 14, all AP-with-doubled-evens). They are, honestly, the q=14 sieve
plus elementary arithmetic — but that is exactly *why* the constructive disproof
attempts (S555) had to fail, and B(2) is a new unification.

Full LRC(14) is **not** proved. The two most promising next theorems are sharply
posed: **(C′)** *tight ⟹ no multiple of 14* (⇒ tight-witness lattice), and **(D)**
*no-odd-split for e≤6* (⇒ LRC(14) for the even-minority half, conditional only on
the proven LRC(7)).

**Artifacts:** `04-computation/lrc_n14_proof_lite_and_brainstorm_s556.py`,
`05-knowledge/results/lrc_n14_proof_lite_and_brainstorm_s556.out`,
`05-knowledge/results/lrc_n14_tight_lattice_B3_s556.out`. Builds on S554
(even-fold), S553 (V*, tight census), S552o (mod-7), THM-369. New: **HYP-2058**.
