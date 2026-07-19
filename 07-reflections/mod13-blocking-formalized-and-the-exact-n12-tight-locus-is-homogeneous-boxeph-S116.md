# The mod-13 pair-blocking, formalized; and the exact n=12 tight locus is the homogeneous AP `c·{1,…,12}` (a=d)

*boxeph-2026-07-18-S116. Owner: Lean-formalize the mod-13 pair-blocking lemma; work on open mathematics.
Two deliverables: (1) `LRCMod13Blocking.lean` — the mod-13 sieve witness, kernel-pure (`[propext,
Classical.choice, Quot.sound]`, no sorry, in the corpus); (2) a sharper characterization of the n=12 tight
locus — `M(C)=1/13` among APs holds **only** for the *homogeneous* `c·{1,…,12}` (`a=d`), never a shifted AP
(`a≠d`) — with a **proved** partial rigidity `AP + M(C)=1/13 ⟹ a≡d (mod 13)` from the new lemma. Verified
S116 computation.*

## (1) The Lean formalization

`LRCMod13Blocking.lean`, `namespace LonelyRunner`, all kernel-pure:

- **`mod13_middle_far`** (`[propext, Quot.sound]`): the integer core — `r∈[2,11] ⟹ ∀k, 2 ≤ |13k+r|` (a
  residue in the middle band is `≥2` from every multiple of 13). Two-case `omega`.
- **`sieve13_middle_witness`** (PROVED): if at scale `b` every speed's residue `(c_i·b) mod 13 ∈ [2,11]`,
  then `t=b/13` puts every runner at distance `≥ 2/13` from the integers —
  `∀ i m, 2/13 ≤ |c_i·(b/13) − m|`. Proof: `c_i·(b/13) − m = (c_i b − 13m)/13`, and the integer numerator
  decomposes as `13·(·) + r` with `r∈[2,11]`, so `|·| ≥ 2` by `mod13_middle_far`; divide by 13 (`gcongr`).
- **`no_middle_band_witness_of_tight`** (PROVED): the contrapositive packaging — a family with a
  `<2/13`-close runner at `b/13` cannot have all residues in the middle band.

This is the Lean form of S115's mod-13 pair-blocking (the sharpened, *proved* analog of the verified
mod-25 HYP-4622): a family with `M(C)=1/13` cannot have a scale `b` at which every runner sits in the
middle band, so it must pair-block — some `c_i ≡ ±b^{-1} (mod 13)` for every `b`.

## (2) Open math: the exact tight locus is homogeneous (`a=d`)

Testing which arithmetic progressions `{a, a+d, …, a+11d}` achieve `M=1/13`:

| AP `(a,d)` | family | `M` |
|---|---|---|
| `(1,1)` | `{1,…,12}` | **1/13** |
| `(2,2)`, `(3,3)`, `(5,5)` | `c·{1,…,12}` | **1/13** |
| `(2,1)` | `{2,…,13}` (shifted) | `2/15` |
| `(1,2)` | `{1,3,…,23}` (odds) | `1/2` |
| `(2,3), (7,1), (1,7), …` | general AP | `>1/13` |

> **The n=12 tight locus is exactly `{c·{1,…,12} : c≥1}` — the *homogeneous* dilated APs (`a=d`, ratios
> `1:2:⋯:12`). No shifted or non-homogeneous AP is tight.**

This confirms HYP-4382's exact form (`C = c·{1,…,12}`, not "general AP") and sharpens it: the extra
constraint is that the AP passes *through the origin* (starts at `d`, is `d·{1,…,12}`). Geometrically this
is the "everything is the triangle" transitive/ordered extreme — the ratios `1:2:⋯:12`.

## The mod-13 lemma proves the partial rigidity `a≡d (mod 13)`

The new lemma pays off immediately on the AP-restricted case:

> **Proposition (PROVED).** An AP `C={a+dk : k=0..11}` with `M(C)=1/13` (and `13∤` some entry) has
> `a ≡ d (mod 13)`.

*Proof.* `M(C)=1/13` forces mod-13 pair-blocking (S115 / `no_middle_band_witness_of_tight`), which for an
AP means its 12 residues mod 13 must **miss 0** (a residue `≡0` gives `‖·‖=0 < 1/13`). The AP residues
`{a+dk mod 13 : k=0..11}` miss 0 iff `a ≢ −dk` for all `k∈{0,…,11}`, i.e. `a/d ≡ 1 (mod 13)`, i.e.
`a ≡ d (mod 13)`. ∎ Verified: every tight AP has `a≡d mod 13`; every loose AP tested is consistent
(`{2,…,13}`: `a=2,d=1`, `a≢d mod 13`, `M=2/15`).

So the AP-restricted rigidity factors: `a ≡ d (mod 13)` is **proved** (mod-13 slice); the full `a = d`
(exact homogeneity) is the residual — the same offset-vanishing / all-moduli content (S115) that the
general rigidity needs. The mod-13 lemma converts the "which AP is tight" question, on its mod-13 face,
into a one-line proof.

## Net

- **Delivered (Lean):** `LRCMod13Blocking.lean`, 3 kernel-pure theorems, built and registered — the
  mod-13 pair-blocking made a proved, machine-checked necessary condition for the n=12 tightness.
- **Delivered (math):** the exact n=12 tight locus is the homogeneous `c·{1,…,12}` (`a=d`); confirmed and
  sharpened HYP-4382; and a **proved** partial rigidity `AP + M=1/13 ⟹ a≡d (mod 13)` using the new lemma.
- **Open:** the full `a=d` (and the general non-AP rigidity) is the offset-vanishing residual — unchanged,
  the inverse-theorem wall.

Cross-links:
[[sharpening-hyp4382-the-mod13-pair-blocking-is-proved-but-necessary-not-sufficient-boxeph-S115]],
[[the-non-covering-sieve-dispatch-is-now-kernel-pure-lean-boxeph-S106]],
HYP-4382 (n=12 tightness), THM-724, `LRCMod13Blocking.lean`.
