# The "free sieve window" lever is refuted: two small gaps force a better time at a *medium* modulus (13 < q′ < q), not a sieve window

*boxeph-2026-07-18-S102. Owner: prove the S101 lever "two interior small gaps force a free sieve
window." Outcome: **the lever is false as stated** — I exhibit a 2-small-gap configuration that is
**sieve-complete** (no witness at any `q′≤13`) yet is beaten at `q=24`. So the crux's forcing operates at
**medium moduli `13<q′<q`**, not at the sieve moduli `q′≤13`. This corrects my own S101 naming and
localizes the inverse theorem to the medium-modulus band — where the additive/CF-depth content lives, not
the elementary sieve. The crux (= open LRC(14), S94) is **not** proved. Verified S102 computation.*

## The lever, and its refutation

S101 named the finishing lever: *two interior small gaps (at the maximizer modulus `q`) force a free
sieve window at some `q′≤13`, hence a witness `M≥1/13`, contradiction.* Testing it directly:

Take `V = {14,15,28,42,56,70,84,98,112,126,140,154,169}` — the deep-well residue AP with a twin inserted
(`+15`, `−168`), giving **two small gaps** (`1` and `13<val=14`) at `q=183`. Then:

| search band | best `M` | at |
|---|---|---|
| `q′ ≤ 13` (**sieve windows**) | **`0`** | — (every sieve modulus killed: `13∣169`, `11∣154`, …) |
| `14 ≤ q′ ≤ 26` | `1/12` | `q=24` |
| all `q′` | `1/12` | `q=24` |

So `V` is **sieve-complete** — it has *no* witness at any `q′≤13` — yet `M(V)=1/12 > 1/13`, realized at
`q=24`. **Two small gaps did not force a sieve window.** The lever, as I named it in S101, is false: the
mechanism that produces `M≥1/13` lives at a *medium* modulus `q′>13`, decoupled from the sieve.

(Other 2-gap configs — e.g. `twin lo+hi` — *are* beaten at `q=11≤13`; but the point is that some are not,
so "2 gaps ⟹ sieve window" is not a theorem. The beating modulus is `11, 24, 26, …` with no uniform
`q′≤13` rule.)

## What this actually shows: the crux is a *medium-modulus* phenomenon

The three modulus regimes are now cleanly separated:

- **Sieve moduli `q′≤13`** (elementary, S101 lemma): `M<1/13 ⟹ q′∣` some speed. This is *sieve-completeness*
  and it is **not** enough — a sieve-complete family can still have `M>1/13` beaten at `q′>13`.
- **The maximizer modulus `q≈13val`** (large): where the residue band and the gap structure live.
- **Medium moduli `13<q′<q`**: **this is where the crux forces.** Two small gaps at `q` are undone by a
  better time at a medium `q′`. Proving "≥2 gaps at the maximizer ⟹ a better time at some medium `q′`" is
  exactly the deep-well isolation (`non-AP ⟹ shallow CF ⟹ M>1/13`, the definitions' `q*≤50` shallow
  binding), and it is the open inverse theorem.

This sharpens S101's obstruction: difference-closure is not only **non-variational** (invisible to
perturbing `t*`, S101) but also **non-sieve** (invisible to the `q′≤13` completeness). The content sits
in the medium band `13<q′<q` — precisely the "any power-saving / additive-dimension" regime the project
has repeatedly localized to (klein-S279 `≥6`-linear; S92 additive dimension 2), and it is untouched by
either elementary tool.

## Why the medium band is the hard part (structural)

A witness at medium `q′` needs `min_v ‖a′v‖_{q′} ≥ q′/13`, i.e. all 13 residues `a′v mod q′` avoid the
band of half-width `q′/13`. As `t` moves from the maximizer `a/q` to a medium `a′/q′`, the residues do
**not** translate rigidly — each `v_i` moves at its own rate `v_i` — so a "hole" in the residues at `q`
does not become a "hole around 0" at `q′`. The map between the `q`-picture (gaps) and the `q′`-picture
(band-avoidance) is the arithmetic of continued-fraction descent, not a geometric rotation. That is why
no local/geometric construction (rotation, translation, single-window) reaches it, and why the honest
statement remains an **additive inverse theorem** on the residue set, not a sieve or packing argument.

## Net (honest)

- **Refuted:** the S101 "free sieve window" lever. A sieve-complete 2-gap family (`M=1/12` at `q=24`)
  shows two small gaps do not force a `q′≤13` window.
- **Localized:** the crux's forcing lives at **medium moduli `13<q′<q`** — not the sieve, not the
  maximizer. This is the additive/CF-depth regime (deep-well isolation), and it is the open inverse
  theorem (= LRC(14), S94).
- **Not proved:** the crux. I will not fabricate a proof of an open problem; this session removes a false
  lead (sieve) and pins the true content to the medium band.

The reduction chain's honest terminus: LRC(14) covering case ⟺ *at the maximizer, ≥2 small gaps force a
better time at a medium modulus* ⟺ the deep-well isolation ⟺ Tao's `n=12` inverse theorem — open, and
provably not reachable by the sieve (this session) or by maximality (S101).

Cross-links:
[[the-route-B-crux-is-the-open-inverse-theorem-what-covering-gives-and-why-maximality-cannot-finish-boxeph-S101]],
[[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]],
[[the-abandoned-attempts-decoded-the-crux-is-additive-dimension-two-not-any-scalar-boxeph-S92]],
[[both-lrc14-routes-bottom-on-the-same-multilinear-cancellation-not-one-mollification-klein-S279]].
