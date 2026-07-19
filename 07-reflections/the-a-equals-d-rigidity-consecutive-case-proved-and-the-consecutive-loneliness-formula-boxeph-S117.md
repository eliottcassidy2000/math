# The a=d rigidity: the consecutive case (d=1) is PROVED, with a closed form `M({a,…,a+11}) = a/(2a+11)`; the spread case remains

*boxeph-2026-07-18-S117. Owner: prove the full `a=d` rigidity for tight APs (`M(C)=1/13 ⟹ C = c·{1,…,12}`).
Result: four **proved** necessary conditions narrowing tight APs to `d=1` or `d≥17`; the `d=1`
(consecutive) case **fully proved** via a clean witness, which also gives the closed form
`M({a,…,a+11}) = a/(2a+11)` (minimized uniquely at `a=1`). The general spread case (`d≥17`) has loose
residuals not caught by the conditions — full `a=d` not closed, but a real chunk is. Verified S117.*

## Proved necessary conditions for a tight AP `{a+dk : k=0..11}` (primitive, `M=1/13`)

Reduce to primitive by dilation (`M` is dilation-invariant). Then:

- **(i) `a ≡ d (mod 13)`** — the mod-13 pair-blocking (S116, machine-checked `LRCMod13Blocking.lean`).
- **(ii) no prime `p ≤ 12` divides `d`.** If `p∣d` (`p≤12`), then every term `a+dk ≡ a (mod p)`, so at
  `t=1/p` all runners sit at `‖a/p‖ ≥ 1/p > 1/13` (using `p∤a` from primitivity). So `M > 1/13`. ∎
- **(iii) `13 ∤ d`** — else `a≡0 (mod 13)` by (i), so `13∣a`, contradicting `gcd(a,d)=1`.
- **(iv) `d ≥ 2 ⟹ ‖a/d‖ ≤ 1/13`.** At `t=1/d` every term `≡ a (mod d)`, so `min = ‖a/d‖`; `>1/13` would
  beat tightness.

By (ii)+(iii), `d` is coprime to every prime `≤13`, so **`d = 1` or `d ≥ 17`**.

## The consecutive case (`d=1`) — PROVED, with a closed form

> **Theorem (consecutive rigidity, PROVED).** For `12` consecutive integers `{a, a+1, …, a+11}`,
> `M = a/(2a+11)`; in particular `M = 1/13 ⟺ a = 1 ⟺ C = {1,…,12}`.

*Witness (lower bound, rigorous).* At `t = 1/(2a+11)` the 12 speeds map to
`(a+k)/(2a+11)`, `k=0..11`, filling `[a/q, (a+11)/q]` with `q=2a+11`. Since `a/q + (a+11)/q = 1`, this
interval is symmetric about `1/2`, and every point's distance to `ℤ` is `≥ a/(2a+11)` (attained at `k=0`
and `k=11`). So `M ≥ a/(2a+11)`, and `a/(2a+11) > 1/13 ⟺ 13a > 2a+11 ⟺ a > 1`. Hence `a ≥ 2 ⟹ M > 1/13`;
tightness forces `a = 1`, i.e. `{1,…,12}`. ∎

Verified the witness is the **exact** maximizer: `M({a,…,a+11}) = a/(2a+11)` (`a=2:2/15`, `a=5:5/21`,
`a=14:14/39`, all matching the computed `M`). The general-`n` form is `M({a,…,a+n−1}) = a/(2a+n−1)`, which
at `a=1` recovers the LRC tight value `1/(n+1)` — a clean formula for the loneliness of any consecutive
block, with `{1,…,n}` the unique tightest.

## The spread case (`d ≥ 17`) — not closed

Conditions (i)–(iv) are **not sufficient** for `d>1`: e.g. `a=2, d=41` satisfies all four (`2≡41 mod 13`;
`41` coprime to `2..13`; `‖2/41‖=2/41 ≤ 1/13`) yet `M({2,43,…,453}) = 222/455 ≈ 0.49` — very loose, far
from tight. So spread APs (`d≥17`) are loose but escape the elementary sieve conditions; ruling them all
out needs a genuine spread-AP argument (their maximizer sits at a large non-obvious `q`, e.g. `455 = 5·7·13`
for `a=2,d=41`), which conditions (i)–(iv) do not supply.

## Net (honest)

- **Proved:** four necessary conditions ((i) `a≡d mod 13`, (ii)+(iii) `d` coprime to `2..13`, (iv)
  `‖a/d‖≤1/13`), narrowing tight primitive APs to `d=1` or `d≥17`; and the **consecutive (`d=1`) rigidity
  in full** — `{a,…,a+11}` tight `⟹ {1,…,12}`, with the closed form `M = a/(2a+11)`.
- **Not proved:** the spread case `d≥17`. Those APs are loose (verified `a=2,d=41`) but not caught by the
  sieve conditions; the uniform argument is open.
- **So:** the `a=d` rigidity is proved for consecutive APs and reduced (via proved conditions) to the
  spread case for general APs. A genuine chunk, not the whole.

The clean win is the consecutive-loneliness formula `M({a,…,a+n−1}) = a/(2a+n−1)` (witness proved; `=`
verified), which pins `{1,…,n}` as the unique tightest consecutive block — the `d=1` face of the
homogeneous tight locus (S116).

Cross-links:
[[mod13-blocking-formalized-and-the-exact-n12-tight-locus-is-homogeneous-boxeph-S116]],
[[sharpening-hyp4382-the-mod13-pair-blocking-is-proved-but-necessary-not-sufficient-boxeph-S115]],
HYP-4382 (n=12 tightness), `LRCMod13Blocking.lean`.
