# The Erdős-592 seam is char-2, the seam tracks the Schur arity, and strong Specker cannot come from valuations — it needs Larson's partial sums

**Source:** kind-pasteur-2026-06-16-S4. Dispatch: use the snippet (the (3,3)-Chang run, the
ternary sweep, the handoffs "climb to a t-uniform rung → constructive strong Specker", "explain
the v₂/v₃ asymmetry algebraically", the THM-460 m=3 enumerator) as inspiration, and **mine the
repo for esoteric things that fell through the cracks**. Four parallel Explore sweeps ran; the
highest-leverage orphan — surfacing in **three of the four** — was **T778** ("the seam tracks
the constraint arity," logged once, "speculative, no computation"). Resolving it gives the
algebraic backbone the snippet asked for, plus a barrier that reframes the strong-witness fork.

## The advance (THM-521): the seam is char-2, and tracks arity

The Erdős-592 triangle-free constraint becomes, on gaps, the **2-term Schur relation**
`g(x,z)=g(x,y)+g(y,z)` (gaps add along a chain). A translation-invariant witness grades gaps by
`v_p`; the level `L_0(p)` is the units `{n : p∤n}`. THM-469 found the seam (sum-free ⟺ p=2);
THM-521 gives the reason and the law (both proved, verified k≤12):

> **FULL:** `L_0(p)` is `k`-term-sum-free ⟺ **`p=2` and `k` even**.
> **DIAGONAL:** the degenerate `k·a=b` is killed by `v_p` ⟺ **`p∣k`**.

The crisp algebraic statement of the v₂/v₃ asymmetry (the snippet's explicit handoff):

> **The seam is `char(𝔽_p)=2`** — the unique field in which the all-ones 2-term Schur sum
> `1+1=0` is a non-unit. For `p=2`, odd+odd=even *always* (full sum-freeness); for `p` odd,
> `1+1=2` is a unit, so the minimal Schur triple `(1,1,2)` already sits in `L_0`. The triangle's
> arity is `k=2` (even, and `2∣2`), which is exactly why both the full and diagonal seams land on
> `p=2`.

**T778 resolved with a twist:** the seam *does* track arity — `p∣k` for the diagonal relation
(`p=3` for a 3-term constraint, etc.) — but *full* sum-freeness is **exclusive to char-2 with
even arity**; odd-arity relations have **no** fully-sum-free valuation grading. So "seam = Schur
arity" is true for the degenerate relation and refined (char-2 ∧ even) for the strong property.

## The barrier (HYP-2558): why strong Specker can't be invariant

The payoff for the *open* problem. Invariant (valuation-grading) witnesses for `ω^n↛(ω^n,3)²`
can use **only the `p=2` grading** — the unique fully-sum-free one for the even-arity-2 triangle.
But the `p=2` grading runs out at the **linear wall `t=2n+1`** (HYP-2396; the invariant algebra
"dies at `t=7`" for `n=3`, THM-470's timed-out master run). Hence:

> **A strong Specker witness (`Q(n,t)` SAT for all `t`) cannot come from any valuation /
> translation-invariant grading.** If it exists, it must use **non-invariant, value-dependent**
> structure — Larson's **partial-sum / interaction-scheme** features (the "join algebra," survey
> §5). This *explains* the `t=7` wall (the invariant algebra IS the capped `p=2` seam) and
> sharpens the fork: either no strong witness exists (HYP-2396, and constructive strong Specker
> is impossible), or it lives off the valuation gradings — which is precisely where the
> never-implemented **join algebra** (patterns of values ∪ partial sums) was supposed to look.

So the snippet's two handoffs are one statement: **the t-uniform rung that could give
constructive strong Specker must be the Larson partial-sum rung, not a valuation rung — because
the seam-arity law caps every valuation grading at the linear wall.** That redirects the search
decisively.

## The most promising things that fell through the cracks (mined this session)

Ranked, with the honest note that I verified only THM-521 here; the rest are leads:

1. **The join algebra (Larson partial sums) — the way past the barrier.** Proposed in the
   survey §5 ("patterns of values ∪ partial sums"), **never implemented**. THM-521 §D now says
   this is *the only* route to a strong witness. Highest-leverage open build.
2. **The shift graph `Sh([ω]²)` is the canonical ω² witness** (survey §3.1: the unique maximal
   triangle-free pattern set `{ADJACENT}` *is* the shift graph — the classical triangle-free
   infinitely-chromatic graph, appearing with zero design choice). Its role in chromatic-number
   problems (Hadwiger–Nelson) was never connected. The shift graph is the `k=1` "atom" of the
   gap-Schur ladder; worth a structural study.
3. **`{7,21}` is NOT 2-adic — a useful disambiguation** (T012, the killed `73=(111)_8`
   prediction). The OCF forbidden gaps `{7,21}` are *Helly/realizability* obstructions
   (THM-499 baby-Hodge), NOT a base/arity phenomenon. So the Erdős-592 "2-adic seam" and the
   OCF "`{7,21}` gaps" are **different 2's** ("the three twos": the seam is the *2 of pair-
   composition / Schur arity*; the gaps are the *2 of parity / objects*). Keeping them distinct
   prevents a tempting false bridge.
4. **Never-harvested computations** (Explore agent 1): the **first Chang number** (m=1 M=2
   tower game at (3,3) — timed out, never finished); the **t=7 master run** (4985 CEGAR clauses
   *not persisted* — add disk-logging and resume, or symmetry-break); the **m=3 (2,2) probe**
   (timeouts only — the open case's first concrete data point). These are concrete, bounded,
   restartable.
5. **Cluster-expansion-at-the-ordinal-limit / Borsuk–Ulam tournament parity** (agent 3): two
   genuinely unwritten speculative bridges (tournament partition functions → Hamkins infinite-
   game ordinal values; Rédei parity ← Ky Fan's odd-simplex lemma). Flagged as long-shots.

## One-line synthesis

> The Erdős-592 seam is `char(𝔽_p)=2` (the all-ones 2-term Schur sum vanishes only in 𝔽₂);
> "seam = Schur arity" holds for the diagonal (`p∣k`) and refines to char-2-∧-even for full
> sum-freeness; and because the triangle is even-arity-2, *every* valuation grading is the
> capped `p=2` one — so constructive strong Specker must climb the **Larson partial-sum** rung,
> not a valuation rung (THM-521 §D, HYP-2558).

## Status / honesty

PROVED & VERIFIED: THM-521 (both laws, elementary; computational k≤12). The char-2 reading and
T778 resolution are exact. The §D barrier is proved *given* invariant-witness = valuation-grading
(THM-469) + the HYP-2396 wall; whether a non-invariant strong witness exists is HYP-2558 (open).
The mined leads (join algebra, shift graph, never-harvested runs) are pointers for future
sessions, not results. Cross-links: THM-521, THM-469 (the seam), THM-453 (triangle = 2-term
Schur), THM-470 (the t=7 wall), HYP-2396 (the linear wall), HYP-2558 (the barrier), T778
(resolved), the-three-twos (the 2 of arity vs parity vs branching), the Erdős-592 survey §5
(the join algebra), the three-axes reflection (kind-pasteur-S3). Source: Erdős Problem 592;
Schur's theorem; Larson's interaction schemes; Chang 1972.
