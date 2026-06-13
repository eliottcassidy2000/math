# LRC(14): the evader condition "13|r" is just the cheapest clock failing

**Source:** kind-pasteur-2026-06-13-S2 (THM-497 Cor B′; follow-up to THM-497).

codex's "five evaders" — the hardest known multiple-of-14 configs, of the form
`7·{1,…,12} ∪ {r}` — all satisfy `13 | r`, an empirical necessary condition the
Pisano-quotient reflection characterized as a product obstruction ("missing
Pisano class at 27 × killed 13-clock"). This session found the cheap explanation
of the `13|r` half: it is *exactly* the band-0 lemma (THM-497 B) at shell
`q = 13`.

The dilated core `7·{1,…,12}` never has `13 | 7k` for `k ≤ 12`, so by the band-0
lemma the rational `t = a/13` is a strict witness for the whole config **unless
`13 | r`**. In other words, the only way for the stranger to survive the cheapest
nontrivial clock — the prime `q = 13`, where the band is just `{0}` and the
witness test is pure divisibility — is to be a multiple of 13. The famous
necessary condition is not deep arithmetic; it is the 13-clock refusing to be
dodged by anything that isn't sitting on it. (Verified: `q=13`-witness `⟺ 13∤r`
on 2247 configs, 0 mismatches.)

This reframes the hard core cleanly. A config of the dilated-core shape is loose
the instant it fails to occupy a small clock; to be hard it must *pay* the 13|r
clock (and then the mod-27 resonance for the next band). The hardness is the
*resonant tail* of the cheapest clock, not a separate phenomenon. And it is not
special to the prime 7: the `d=3` core `3·{1,…,12} ∪ {r}` produces its own
climbers (`h = 39` at `r = 182, 364`, `13|r`, with a *different* mod-27 signature
`{13, 20}` rather than the `d=7` family's `{0, ±10}`). So "dilated `{1..12}` core
+ 13|r + mod-27-resonant stranger" is a `d`-parameterized family of hard-but-loose
configs — the evaders are one fiber of it.

Two honest limits. First, this only explains the *necessary* `13|r` condition (the
band-0 floor); the *sufficient* hardness still lives in the band-2 mod-27
resonance, which is the over-correlated multiplicative-character regime THM-497 D
flagged as the open core. Second, the resonance/compressibility exploration that
opened this session mostly *re-confirmed* the known structure (hard core =
structured block + resonant stranger) rather than cracking it — the genuinely new
leverage remains THM-497's cardinality dichotomy (counting can never prove it) and
tool-domain boundary. The value added here is the bridge: the band-0 lemma, a
five-line proof, accounts for half of codex's empirical evader condition, and the
phenomenon parameterizes over the dilation `d`.

Cross-links: THM-497 (band-0 lemma + Cor B′), THM-492 (the five evaders),
HYP-2482/2483/2484, [[lrc14-covering-cardinality-permits-structure-forbids-kps1]]
(the cardinality dichotomy this builds on).
