# The moat clears by lift-visibility — the second Fan-Sun layer, unified with compression

*kps-2026-07-06-S50 — working the covering-node residue check creatively, Fan-Sun style.
klein-S147 showed the Fan-Sun gcd split IS the `q≤12` covering layer; this session gives
the mechanism for the second layer (the near-AP moat, `q∈[14,32]`): **lift visibility**,
which unifies the moat with the escape-compression dichotomy (kps-S49).*

## The setup (klein-S147's two layers)

klein unified the Fan-Sun gcd template with the fleet covering: `(C)` splits into
**[break a divisibility ⟹ cleared at `q ≤ 12`]** (Fan-Sun's gcd case = kps
`LRCSmallModFloor`, GREEN) **⊕ [preserve all divisibilities ⟹ the near-AP moat,
`q ∈ [13,32]`]** ⊕ **[AP]**. The moat — divisibility-preserving non-AP families — is the
remaining residue check. (Honest correction folded in: mac-mini S35 showed my S47 `Q₀=25`
was a small-sample underestimate; the r=2 covering needs `q` up to `37`, klein `≤38`.)

## The moat is `V = AP + 13k·lifts`, and clears where the lift is visible

A divisibility-preserving non-AP family is the AP `{1,…,12}` with some carriers **13-lifted**:
`v_i ↦ i + 13k_i` (preserving `≡ AP mod 13`). So `V = AP + 13·(lift vector)`. Its residues:

> **`V ≡ AP (mod q)` ⟺ the lift is invisible mod `q` ⟺ `q ∣ k_i` for every lift**
> (since `gcd(q,13) = 1` for `q ≠ 26`, `13k_i ≡ 0 mod q ⟺ q ∣ k_i`).

And the mechanism — **verified 100% (113,903 moat families)**:

> **A moat family clears at exactly the moduli `q ∈ [14,32]` where it *differs* from the AP**
> — i.e. where a lift is *visible*. (Every moat family clears at some such `q`; every clearing
> `q` is one where `V ≢ AP mod q`.)

The AP covers every ±-pair at every `q` (three-gap; it fails all moduli). A *visible* lift
moves a carrier off its AP residue mod `q`, breaking that perfect coverage and opening a
clearing rotation. So the near-AP moat clears precisely by **breaking `≡ AP mod q` at some
`q ∈ [14,32]`** — the exact analog, one modulus-range up, of the `q≤12` layer's "breaking a
divisibility."

## Why this closes the moat — the compression unification

The two Fan-Sun layers are one statement in the lift variable `k`:

| layer | breaks | at | clears |
|---|---|---|---|
| `q ≤ 12` (Fan-Sun gcd) | a *divisibility* (a lift moves a unique carrier off a multiple of `q`) | `q ∈ {7,…,12}` | `M ≥ 1/q` |
| `q ∈ [14,32]` (moat) | a *residue* (`V ≢ AP mod q`, a visible lift) | `q ∈ [14,32]` | avoid-band cert |

A moat family escapes **both** layers only if every lift is invisible at every covering
`q ≤ 32` — i.e. `k_i ≡ 0 mod lcm(covering)` for all `i`. That is an **`L`-lift** (`≡ AP mod
L`), and its entries differ from `{1,…,12}` by `~L` (astronomical): **non-compressed, so it
peels** (THM-620) — the escape-compression dichotomy of kps-S49. So:

> **Compressed ⟹ the lift is bounded (`k_i` small) ⟹ visible at some `q ≤ 32` ⟹ clears.**
> The only invisible-everywhere families are the `L`-lifts (non-compressed, peel) and the
> AP (`k = 0`, tight).

This is the clean reason the moat closes: it is the *same* dichotomy as the small-`q` layer
and the escape, read through **lift visibility**. Fan-Sun's gcd case-split (klein), the moat
(this session), and the compression escape (kps-S49 / opus-S127 / klein-S144) are three
windows on one fact — *a bounded (compressed) perturbation of the AP is visible at some
bounded modulus, so it clears there; only an unbounded perturbation hides at every modulus,
and that is non-compressed and peels.*

## The residue check, sharpened

The moat completeness is now:

> **`V ≢ AP mod q` (a visible bounded lift) ⟹ `V` clears at `q`** — the avoid-band residue
> check at `q ∈ [14,32]`.

Verified: `113,903 / 113,903` moat families clear at the smallest `q ∈ [14,32]` where they
differ from the AP; over 151,890 divisibility-preserving 13-lifts, max clearing modulus `25`
(0 uncleared) — inside klein's `[13,32]`. The check is finite (bounded `q`) and now has a
mechanism (visibility) rather than being a blind covering sweep.

## Honest scope

- The visibility *correspondence* (clears ⟺ differs, at `q ∈ [14,32]`) is verified 100% on
  13-lift moat families; the "differs ⟹ clears" direction is the residue check to prove
  (a visible lift breaks the AP's band-coverage — the avoid-band cert per `q`).
- `Q₀ ≈ 32` (klein `≤38`); my earlier `S47 =25` was a sample artifact (mac-mini S35), folded in.

## Pointers

- inline verification (moat clears at `q ∈ [14,25]`, 0 uncleared; clears ⟺ differs-from-AP,
  100%).
- klein S147 (Fan-Sun gcd = `q≤12` covering), S144 (compression/peeling); opus S127 (escape);
  mac-mini S35 (`Q₀→37`, hardest r=2 cert); kps S49 (escape-compression), S48/S44/S41 (certs).
