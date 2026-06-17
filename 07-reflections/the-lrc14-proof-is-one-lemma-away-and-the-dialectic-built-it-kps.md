# LRC(14) is one lemma away — and the prove/disprove dialectic is what built the skeleton

**Source:** kind-pasteur-2026-06-17-S1. Dispatch: spend a long session aiming to BOTH prove and
disprove LRC(14), using each goal as a novel-idea source for the other. I ran a dialectical
workflow (7 agents: tight-locus census, literature, disprove-probe, prove-gap → prove-synthesis
& disprove-synthesis each consuming the other's findings → adversarial verify), and developed the
prove-side coverage-deficit core independently. The result is a **near-complete singular-series
proof of LRC(14), reduced to one well-isolated lemma** — and the reduction was forged precisely
by the prove/disprove cross-pollination.

## Where it stands (LRC(14) = the FIRST OPEN case)

Literature pinned: LRC is now proven for ≤12 speeds (Sungkawichai–Trakulthongchai 2026); **our
13-speed / gap-1/14 case is genuinely the live frontier.** The route here is the project's own
**lonely-measure / singular-series** object `L(S)=meas{τ:||vτ||>1/14}` (`inf L>0 ⟹ LRC(14)`),
which is *orthogonal to the literature*: every published Lonely-Runner lower bound (incl. Bedert
2025) controls the **gap** `κ(n)`, none controls the **measure** `L`. So the inf-`L` bound is
original territory.

## The proved skeleton (THM-523)

1. **Quantization (THM-522):** `L>0 ⟹ L≥1/(14·lcm)` — `L` is never an infinitesimal.
2. **Decoupling floor (PROVED):** `L(C∪{w}) ≥ (6/7)meas(G_C) − r/(7w)`; one speed `→∞` pushes
   `L` *up* to `(6/7)meas(G_C)`, min floor `1/143` at the drop-6 core — **the single-element
   escape to 0 is closed.**
3. **Single-perturbation infimum = 1/1260 (PROVED, explicit `N≤93`):** `L<1/1260 ⟹ w≤93`, then
   exhaust. The sharp value is a **two-speed clash** (speeds 5 and 36): `15/36−2/5−1/70−1/504 =
   1/2520`, doubled by `τ↔1−τ` to `1/1260`. `w=24` covers (tight, = Goddyn–Wong's T5).
4. **Zero counterexamples (PROVED on the tight locus):** the only tight configs (`{1..13}` and
   `{1..11,13,24}`) have `max-min = 1/14` *exactly*; LRC(14) holds with equality, and the
   max-min spectrum has a gap `1/14 < 3/41 < 2/27 < …`.

## The one remaining lemma (OPEN-Q-108)

> **Uniform fattening lemma:** `∃c>0` with `meas(G_C) ≥ c` for *every* 12-subset `C` of distinct
> positive integers — equivalently, the **primitive tight locus at n=13 is finite**.

This is the *entire* gap. Decoupling handles speeds growing one at a time; the **only
uncontrolled regime is `k≥3` arithmetically-coordinated growing speeds.** The literature confirms
the fixed-`n` tight-locus classification is "widely open" (Perarnau–Serra) — Goddyn–Wong's only
*infinite* tight family needs `n→∞` (two accelerations already need `n≥1.2·10^14`), so no
published `n=13` family with `lcm→∞` exists to threaten finiteness.

## The dialectic — each goal generated the other's ideas (the dispatch's point)

This is the part worth keeping. The cross-pollination wasn't decoration; it built the proof:

- **DISPROVE → PROVE (the key correction):** my THM-522 said "small `L` ⟹ bounded `lcm`." The
  disprove side *killed that* by exhibiting near-tight 12-cores with `lcm ~ 10^7` (`L=2/425`). That
  forced the proof onto the **right object** — a uniform bound on `meas(G_C)` (bound the
  *perturbing elements* via decoupling), not the lcm. The adversary corrected the proof's target.
- **DISPROVE → PROVE (localization):** the disprove search found the drop-6 family minimizes at
  the *large* `w=69` (not the obvious small multiple), pinpointing the **only uncontrolled regime**
  (`k≥3` coordinated growth) — telling the prove side exactly where the remaining lemma must work.
- **PROVE → DISPROVE (closing escapes):** the decoupling floor (single element `→∞` ⟹ `L≥1/143`)
  shut the single-element escape, *forcing* the disprove side into the multi-coordinated regime —
  where it also found no foothold (best `≥0.011`, `L` always *rising* with lcm).
- **BOTH → BOTH (the lever):** the **LRC(12) lever** (exactly one 12-subset of `{1..14}` is tight
  at gap `1/13`, zero at `1/14`) converts the crux from *existence* ("is `meas(G_C)>0`?") into
  *transversality* ("does the isolated gap-`1/13` maximizer **fatten** to a positive gap-`1/14`
  measure, uniformly?"). Using the *proven* LRC(12) as a lemma for LRC(14) is the cleanest new
  handle the dialectic produced.
- **The adversary's strongest evidence is on the prove side's terms:** the disprove side's
  repeated finding that **`L` increases with lcm** (`1/1260 → 0.0062 → 0.0078 → 0.0094`) *is* the
  quantitative content of stranger-decoupling. The hardest attack corroborates the proof.

## Honest bottom line

LRC(14) is **not proved and not disproved**. But the singular-series route is now a complete
proof *modulo one quantitative transversality lemma* (uniform `meas(G_C)` ≡ tight-locus
finiteness), and every exact computation across 320k+ configs points one way: `inf L = 1/1260`,
LRC(14) TRUE with equality on a finite tight locus `{AP, Goddyn–Wong T5}`. The dialectic did its
job — the disprove side, by failing precisely, sharpened the prove side onto the right object and
isolated the single open step. Cross-links: THM-523, THM-522 (corrected), THM-518, THM-501,
HYP-2561, OPEN-Q-108, MISTAKE-075; Goddyn–Wong (Integers 2006 #A38), Tao (arXiv:1701.02048),
Sungkawichai–Trakulthongchai (arXiv:2604.23906), Perarnau–Serra (arXiv:2409.20160); scripts
`04-computation/tight_locus_{census,largescan,final,maxentry_probe}_kps.py`,
`lrc14_exact_lonely_measure_kps.py`.
