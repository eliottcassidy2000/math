# The non-dilated-core rigidity residual is Tao's n=12 conjecture — the definitive frontier

*boxeph-2026-07-18-S114. Owner: prove the non-dilated-core rigidity residual. Honest, definitive outcome:
this residual **is** the open crux — Tao's n=12 optimistic conjecture, equivalent to LRC(14)'s covering
case (S94), height-unbounded (not a finite check), and provably beyond every tool the project has (S101–
S113). I did not prove it and will not fabricate a proof. This reflection consolidates the frontier: what
the residual is, why it resists, what surrounds it (all done), and what a proof would actually require.*

## What the residual is

> **Non-dilated-core rigidity:** `M(V) < 1/13` (covering) ⟹ the 12-core `V∖{v_max}` is a dilated AP
> `d·{1,…,12}`. Equivalently (S94): the offset-vanishing form; equivalently the n=12 inverse theorem
> `|C|=12, M(C)=1/13 ⟺ C = d·{1,…,12}` (HYP-4382), the equality/rigidity case of settled LRC(13).

It is **height-unbounded**: the tight cores are `d·{1,…,12}` for every `d` (verified: `M(d·{1..12})=1/13`
for `d=1,2,3,5,7,…`), so it is a genuine infinite rigidity statement, not a finite verification. It is
**sharp**: boundary families attain `M=1/13` exactly (S113). It is **equivalent to LRC(14)** (S94, proved),
so no reformulation is a strictly weaker sub-problem.

## Why it resists — the systematic exhaustion (S101–S113)

Across thirteen sessions the entire standard toolkit has been shown, each with its specific reason, to fall
short of this one residual:

| tool | reaches | provably stops at | ref |
|---|---|---|---|
| maximality / perturbation | active-runner pinning at `t*` | blind to interior small gaps | S101 |
| sieve completeness | `q'∣` some speed, `q'≤13` | sieve-complete families beaten at `q'>13` | S102 |
| continued-fraction descent | `lcm(13,14)∣v_max` | far-element divisibility only | S103 |
| BSG / PFR (additive) | structure *from* energy | needs energy input `M<1/13` doesn't supply | S104 |
| the `1/12` gap theorem | — | *stronger* than the crux, sharp, not a lever | S111 |
| descent recursion | large-`ρ` (single-killer) | loses factor ~2 at `ρ≈1` (compact) | S113 |

The unifying diagnosis (S104): the residual is a **Diophantine→additive-energy bridge** — "global rational
optimality ⟹ the core has maximal additive energy (is an AP)." The elementary tools never reach the
additive core; the additive tools presuppose it; the two sit on opposite sides of the one missing
implication. That implication is a **concentration/transference** statement, and Weyl provably cannot force
concentration (S95). No tool in the kit produces it.

## What surrounds it — all done or reduced

The residual is now *isolated*; everything else is settled or kernel-checked:

- **Non-covering** ⟹ sieve witness — kernel-pure Lean (`sieve_dispatch`, S106).
- **`≥2` outliers** ⟹ `M ≥ 1/13` — THM-726.
- **Single-killer** (val=14 / `ρ≥13`) ⟹ `M > 1/14` **unconditional** — THM-1007 (3-line balance lemma);
  the sharp `14/183`-uniqueness is THM-724 mod the *same* residual (S112).
- **Dilated-AP-core compact** ⟹ `M ≥ 1/13` — THM-1013 (dilated sieve), kernel-pure Lean.
- **The full reduction** `LRC(14) ⟸ LRC(≤13) + INVcov` — kernel-checked, down to the ledger's own target
  Prop (S105–S109; INVcov = this residual).
- **Sharpest concrete form of the residual:** the mod-25 **pair-blocking rigidity** (HYP-4622): a 12-set
  that blocks all ten unit `±`-pairs mod 25 and has `M < 2/25` is the AP — a concrete covering/linear
  constraint plus Freiman tightness, verified over `~150k` families but height-unbounded.

So LRC(14) rests on exactly one wall, approached from two sides (single-killer and compact, S112/S113):
the near-dilated-core rigidity.

## What a proof would actually require

Not another reformulation — S94 proved they are all equivalent. A proof needs a genuinely new mechanism
for the Diophantine→additive bridge: force **additive energy / AP-structure on the residues from global
rational optimality**. Candidate shapes, none in the current kit:

1. an `L⁴`/Fourier **concentration** bound: global optimality ⟹ a large Fourier coefficient of the residue
   set ⟹ `E(C)` maximal ⟹ (Freiman) AP — but concentration is exactly what Weyl cannot give (S95);
2. an **equality characterization** of settled LRC(13): if the Sungkawichai–Trakulthongchai proof pins the
   extremal `M(C)=1/13 ⟺ AP`, HYP-4382 follows and the residual collapses — this is the most likely route,
   and it lives inside a proof I do not have access to;
3. a **transference** turning "no better rational at any `q'`" into an additive-energy inequality — a new
   lemma, not standard additive combinatorics.

Each is a research-level result. This is Tao's n=12 optimistic conjecture; it is open in the literature.

## Net (honest, and direct)

I did not prove the residual, and no session-length effort will, because it **is** the open conjecture. The
project's genuine, durable output is everything *around* it — the complete reduction (kernel-checked to the
ledger target), the single-killer closure, the dilated-sieve and descent tools, and the precise map of why
each standard tool falls short. Those are done. The core is a named open problem awaiting a new idea, most
plausibly the equality case of settled LRC(13) (HYP-4382). Further "prove the crux via X" requests will
reach this same wall — the honest, useful next moves are elsewhere: sharpening HYP-4382 toward the LRC(13)
equality proof, engineering deliverables, or other-`n` LRC.

Cross-links:
[[the-compact-case-is-equivalent-to-LRC14-and-descent-is-the-wrong-tool-boxeph-S113]],
[[INV-val14-is-the-single-killer-case-essentially-done-not-the-open-compact-crux-boxeph-S112]],
[[bsg-pfr-attack-the-wrong-half-the-crux-is-the-diophantine-to-energy-bridge-boxeph-S104]],
[[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]],
HYP-4382 (n=12 tightness), HYP-4622 (mod-25 pair-blocking rigidity), THM-724/1007/1013.
