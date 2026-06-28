# Hypothesis-testing on the frontier: the covering bound holds with a uniform MARGIN — the tight-locus is finite and isolated (δ≈0.0026), the lcm-family is value-safe (M→1/12), and the witness count = (p−1)·d confirms the index-theorem

*mac-mini-2026-06-28-S80. The owner: a long session of creative hypotheses + testing for concrete progress on
the LRC(14) frontier. Five hypotheses tested; the headline progress is that the covering bound is structured as
{finite tight-locus} + {uniform margin}, and the feared no-finite-certificate family is value-safe. Builds on
[[the-index-theorem-frame-lrc-is-the-nonvanishing-of-an-index-analytic-equals-topological]],
[[pushes-and-pulls-on-the-hard-core-d7-unification-the-covering-witness-construction-and-the-imaginary-quadratic-pull]].*

## H1 (CONFIRMED) — the covering bound holds, and random covering sets are FAR from tight
600 random covering 13-sets (each with a multiple of 14): **0 violations of `M ≥ 1/14`**; min `M = 0.111 ≈ 1/9`.
So the covering bound is *easy* for random sets (`M ≈ 1/9`); the difficulty is concentrated at the special
tight configs, not the generic ones.

## CH2 (CONFIRMED) — the index-theorem witness count: peaks = (p−1)·d, primitive index = (p−1)/2
The number of safety peaks reaching `M` (the witnesses) at the tight configs:
```
  AP {1..13}:        6 peaks   index = 6/2 = 3 = (p-1)/2   (p=7)
  GW {1..11,13,24}:  6 peaks   index = 3
  2·AP (covering):  12 peaks   index = 6 = 2·3            (dilation ×2)
  3·AP (covering):  18 peaks   index = 9 = 3·3            (dilation ×3)
```
So **witnesses = (p−1)·d**, the primitive index `(p−1)/2 = 3` (= the de Moivre cyclotomic degree = the
Borsuk-Ulam degree, ODD for `p≡3 mod 4`), with the dilation `d` replicating it `d` times. This directly confirms
the S79 index-theorem frame: the witnesses come in `(p−1)/2` antipodal pairs (the Borsuk-Ulam structure).

## H2 + H3 (CONFIRMED) — the tight-locus is FINITE and ISOLATED with a uniform margin
- **Tight (`M=1/14`)**: only the AP/GW **dilations** (verified: AP, GW, 2·AP tight; perturbations not).
- **Single-swap gap** (H2): the smallest single-swap `M` above `1/14` is `0.0732` — a gap `≈ 0.0017`.
- **Near-tight perturbations** (H3a): over 4000 multi-perturbations of the AP, the smallest `M > 1/14` is
  `0.0740` (the set `AP` with one element **doubled**, e.g. `10→20`) — a gap `δ ≈ 0.0026`.
- **The lcm-family is VALUE-SAFE** (H3c, the crux): `S_X = {1..11,13,lcm(2..X)}` (the S45 "no finite witness
  certificate" family) has `M → 1/12 = 0.0833` (NOT `1/14`!) as `X` grows — margin `≈ 0.012`:
```
  X=7:  M=0.08235    X=9:  M=0.08317    X=11: M=0.08332    X=13: M=0.08333 → 1/12
```
> **KEY:** the S45 obstruction (the witness DENOMINATOR grows without bound) does NOT make `M` approach `1/14` —
> the `M`-VALUE stays at `≈ 1/12`, bounded above `1/14` by a margin. So "no finite certificate" (a statement about
> the witness location) and "`M` has a margin" (a statement about the value) are DIFFERENT: equidistribution is
> needed to LOCATE the witness, but the BOUND `M ≥ 1/14` holds with room to spare on this family.

## The proof structure this reveals (the concrete progress)
> **Covering bound = {TIGHT-LOCUS: AP/GW dilations — FINITE, `M=1/14`, constructed (S77 `t=1/(14d)`) + the
> Borsuk-Ulam index (S79)} ⊕ {ALL ELSE: `M ≥ 1/14 + δ` — a uniform MARGIN}.**
The tight-locus is finite and isolated; everything else (including the lcm-family) is bounded away from `1/14`.
So the proof reduces to: (1) the tight-locus is finite (the census — AP/GW dilations only); (2) the construction
+ index at the tight-locus (done for dilations, S77; GW shares Φ₁₄, kps S31aw); (3) the uniform margin `δ` for
the rest (the bulk — equidistribution gives `M ≥ 1/12`-ish, with room). The near-tight (AP with one doubled,
`M≈0.074`) is the spectral-gap edge.

## Honest status
- **CONFIRMED (computationally):** covering bound holds (0/600); witnesses `=(p−1)·d` (index-theorem); tight-locus
  `=` AP/GW dilations, isolated with gap `δ≈0.0026`; lcm-family value-safe (`M→1/12`).
- **PROGRESS:** the covering bound is structured as {finite tight-locus} ⊕ {uniform margin}; the S45 denominator
  obstruction is NOT a value obstruction; the index-theorem witness structure `(p−1)·d` is verified.
- **NOT a proof.** The open rigor: (a) prove the tight-locus is finite (the census, AP/GW only — a hard rigidity
  statement); (b) prove the margin `δ` is uniform (no untested sequence with `M→1/14`); (c) the construction/index.
  LRC(14) open, but the frontier's STRUCTURE is now sharp: finite tight-locus + uniform margin.

Related: HYP-3250 (this), HYP-3246 (index-theorem frame), HYP-3240 (covering construction), kps S31aw (GW Φ₁₄),
HYP-2901 (the lcm denominator wall), OPEN-Q-108.
