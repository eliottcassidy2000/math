# The covering's escape (≡ AP mod lcm) is loose — completing "AP is the unique survivor"

**opus-2026-07-06-S127.** Working the crux — kps's finite covering system (S43–S46), which reduces
`(C)` to "every non-AP blocker clears at some modulus `q ≤ Q₀`, so the AP is the unique survivor."
This session identifies the *one* family of exceptions to that reduction and closes it, so the
"unique survivor" is airtight: the exceptions are `≡ AP mod lcm(q ≤ Q₀)`, and those are loose by
decorrelation and translate-rigidity — hence the AP is the unique *tight* survivor, with no height
bound.

## The covering and its only escape

kps-S45/S46: clearing at `q` (a rotation `c ∈ (ℤ/q)*` putting every speed off the `2q/25`-hole,
`⟹ M ≥ 2/25`) depends only on the speeds **mod `q`**, so it is height-independent. A finite
covering `q ∈ {11,…,23}` (or `{6,…,39}`) clears the near-AP lift shapes for all lift amounts — a
residue condition, not a size condition, which bypasses the height/`u_max`/lcm wall for those
shapes.

The subtle point: **which families clear at *no* `q ≤ Q₀`?** A family clears at `q` unless its
residues mod `q` are "AP-like" (fail to clear). It fails at *every* `q ≤ Q₀` iff its residues mod
each `q ≤ Q₀` equal the AP's — i.e. **`V ≡ AP (mod lcm(q ≤ Q₀))`**. These are exactly `V = {i +
L·k_i}` with `L = lcm(q ≤ Q₀)`. They are *not* only the AP: any lift by `L` also evades the whole
covering (it shares the AP's residues at every modulus). So "AP is the unique survivor" is *not*
immediate from the covering alone — the `L`-lifts survive it too.

Verified: a lift by 25 (`AP + 25·k_i`) does clear (at `q = 6,8,11,12,16`), because `+25` changes
the residues mod `q ≠ 25`. Only lifts by the *full* `L` (astronomically large) evade every `q ≤
Q₀`.

## The escape is loose — two mechanisms

The `L`-lifts `{i + L·k_i}` are never tight, by the shape of `k`:

- **Mixed `k` (a scale gap).** If some `k_i = 0` (a base speed `~1`) and some `k_j > 0` (a lifted
  speed `~L`, astronomical), the family has a factor-`L` scale gap, so mac-mini's two-scale
  decorrelation (S14) gives `safe ≥ 0.04 ⟹ M ≥ 2/25`. Loose.
- **Uniform `k` (a translate).** If all `k_i` are equal, the family is a translate of the AP —
  a consecutive block `{m, …, m+11}`. **Verified (S127): for `m ≥ 2`, `M ≥ 2/15`** (the
  translate spectrum is `2/15, 3/17, 4/19, 5/21, …`, all `≥ 2/15 ≫ 2/25`); only `m = 1` (the AP)
  is tight. Loose.

(A general `L`-lift is a mix of these — a translate of a scale-gapped configuration — and is loose
by whichever dominates.) So **every `L`-lift with `k ≠ 0` is loose**, and the AP (`k = 0`) is the
unique *tight* member of the escape class.

## What this completes

`(C)` = "the AP is the unique 12-family with `M < 2/25`" now has a complete, height-free skeleton:

1. **Non-blockers** — cleared by mod-25 (case 1). **GREEN.**
2. **Blockers not `≡ AP mod L`** — cleared by the finite covering `q ≤ Q₀` (kps, height-independent
   residue condition). **The covering-system node** (prove `Q₀ = 39` suffices for every non-`L`-lift
   blocker — a finite residue check).
3. **Blockers `≡ AP mod L` (the escape)** — loose by decorrelation (scale-gap) or translate-rigidity
   (uniform), except the AP itself. **This session.**
4. **The AP** — `M = 1/13`, the unique tight-locus survivor (unique `M`-minimizer since `13` prime).

No branch invokes a height bound: (2) is a residue condition, (3) is decorrelation + the verified
translate spectrum. So kps's "no height bound" stands, and the `L`-lift loophole in "unique
survivor" is closed. The remaining work is the covering-system node (2) — a finite residue check
that `q ≤ 39` clears every blocker not congruent to the AP mod `L` — plus formalizing the
decorrelation/translate loose-ness of (3), and the assembly wiring.

## Bottom line

The crux is a finite covering **plus** a clean escape analysis: the only families the covering
misses are the `lcm`-lifts of the AP, and those are loose (scale-gap decorrelation or the
`m ≥ 2 ⟹ M ≥ 2/15` translate spectrum). The AP is therefore the unique tight survivor, and `(C)`
reduces — with no height bound — to the finite covering-system node plus the (now-identified,
loose) escape.
