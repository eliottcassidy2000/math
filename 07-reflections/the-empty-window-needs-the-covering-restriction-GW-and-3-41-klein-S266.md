# The empty window needs the covering restriction — Goddyn–Wong and 3/41 are the witnesses

*klein-2026-07-12-S266. Owner directive: another mining/verification session. I set out to verify
the freshly-consolidated frontier — opus-S246's reframe that LRC(14) ⟺ the single rigidity
`[M(S) < 2/27 ⟹ S is a dilated interval {1..13}]` (HYP-6155, marked CONFIRMED yesterday). Verifying
it broke it: the statement as written is false, in two independent and already-canonical ways, and
the fix is a covering restriction the reframe dropped.*

*Concurrency note: opus self-corrected in S247 (same day) on the first point — "the window (1/14,2/27)
is non-empty, 3/41 realized" — reframing the crux as "AP minimizes M." My verification lands
independently on the same 3/41 and sharpens the reframe: the minimizers are TWO families `{AP, GW}`
(GW also achieves 1/14, so "AP minimizes M" is non-unique), and the operative restriction is
covering — both minimizers, and the whole inhabited window, are non-covering. The precise correct
crux is boxeph's `DC ⟹ M ≥ 1/13`, not "AP minimizes M."*

---

## What I checked, and what it returned

I computed the exact loneliness spectrum `M(S) = max_t min_i ‖v_i t‖` (via the THM-668 pair-sum
ruler, cross-checked on a 2·10⁵ grid) over the small-M region — the only place families near the
tight value live. Two families settle it:

- **The window is not empty.** `{1,2,3,4,5,6,7,8,9,10,11,13,36}` has `M = 3/41 ≈ 0.07317`, which lies
  strictly inside `(1/14, 2/27) = (0.07143, 0.07407)`. Exact (ruler `q = 41 = 5+36`, witnesses
  `t = 17/41, 24/41`) and grid-confirmed. This is the **third-mediant** `3/(3k+2)` at `k=13` — already
  in canon as HYP-2621/opus-S118, which explicitly titled it *"first-gap emptiness is non-monotonic
  in N — arithmetic, not metric."*
- **The tight locus is two families, not one.** The **Goddyn–Wong set** `{1,…,11,13,24}` has
  `M = 1/14` **exactly** (ruler `q=14`, grid 0.071429), is primitive (gcd 1), and is **not** a dilated
  AP (its gaps are `1,1,…,1,2,11`). This is the published Goddyn–Wong 2006 tight family, already canon
  as THM-612/652/658.

So `[M < 2/27 ⟹ dilated interval]` fails twice: `M(GW) = 1/14 < 2/27` with GW non-AP, and
`M(\{…,13,36\}) = 3/41 ∈ (1/14, 2/27)` with that family non-AP. The reframe (HYP-6155) contradicts
its own canon (THM-612, HYP-2621).

## The one thing both witnesses share: they are non-covering

`GW` has no multiple of 14. `{1,…,11,13,36}` has no multiple of 14. Both are **non-covering**, hence
both are dispatched by the elementary `t = 1/14` sieve (THM-366/523) — neither threatens LRC(14),
and both sit at or above `1/14`. That is not a coincidence; it is the whole content. The empty-window
rigidity is true, but **only after the covering restriction**:

> **Corrected statement.** Among *compressed* (max ≤ 13·min) covering families, `M ≥ 1/13`
> (equality on the compressed near-dilate `2·{1..12}∪{13}`). Among *all* primitive families the tight
> locus is `{dilated AP} ∪ {Goddyn–Wong}` (both non-covering), and the sub-tight window `(1/14, 2/27)`
> is inhabited (by non-covering families like the `3/41` third-mediant).

**A second false floor, same genus.** While verifying I also re-broke *boxeph's* replacement crux
`DC ⟹ M ≥ 1/13` (HYP-6150, likewise fresh-CONFIRMED). It is false: `{1,…,12, 182}` is primitive and
divisor-complete (`182 = 14·13` supplies both the missing 13 and 14) with `M = 14/183 ≈ 0.0765 < 1/13`
— the long-known **deep well / covering-min** (MISTAKE-097, opus-S52). boxeph's exhaustive census
missed it because the deep well first appears at `Vmax = 182`, far outside the `Vmax ≤ 30` box — the
exact box artifact MISTAKE-141 warned against, one layer down. The number line, all verified:
`1/14 < 2/27 < 14/183 (deep-well DC) < 1/13 (compressed floor) < 3/37`. So the DC floor is `1/13`
only on the **compressed** subclass (`max ≤ 13·min`); non-compressed DC carries a far element and
**peels** (drop `182` → `{1..12}`, `M = 1/13` by LRC(13)), its own `M = 14/183` claimed only `> 1/14`.
Note `14/183 > 2/27`, so it does *not* break opus's weaker `DC ⟹ M ≥ 2/27` — only the sharp `1/13`.

The empty window that IS clean lives one level down: the **peeled 12-runner base** `(1/13, 2/25)`
(HYP-4151, proved at `r=1` over the mod-13 *field*). The literal 13-runner window is scarred by
`3/41`, which appears precisely because `k−1 = 12` is composite (`6 | 12` fires THM-539's
primorial-gate `a=3` dip), and the doubly-prime peel to the 12-base (`k−1=11` prime, `k+1=13` prime)
is exactly what removes it.

## Why the composite modulus forces the restriction

This is the same apex-7 zero-divisor that makes LRC(14) the first open case. At the prime `k=12`
level, modulus 13 is a field: the tight residues `{v_i·a mod 13}` are 12 distinct nonzero residues,
which *are* the full unit group `(ℤ/13)*`, so tightness forces the AP uniquely — no second family.
At `k=13`, modulus `14 = 2·7` is not a field: only 6 units `{1,3,5,9,11,13}`, with `7` the
zero-divisor. GW inhabits the **non-unit hole** — mod 27 its residues miss the non-unit antipodal
pair `{12,15}` and double the gcd-3 shell `{3,24}` — so it is a *non-transversal* tight family
unreachable by any unit-witness, the exact composite analogue of the `n=6` second tight set
`{1,3,4,5,9}`. The clean "residues form the AP" conclusion of the field case simply does not survive
verbatim to the composite; the surviving object is `M=1/14 ⟹ S ∈ {AP, GW}`, both non-covering.

## The correction, and what it points at

The reframe's *value* stands: `E₃` (Schur triples) is the right translation-sensitive lever, DC is
loose, and the four faces (three-gap, pigeonhole, E₃, Farey-ladder) do unify. What must be struck is
the **unrestricted** rigidity and the "empty `(1/14, 2/27)`" — both need the covering hypothesis, and
opus's own `DC ⟹ M ≥ 2/27` is superseded by boxeph's sharper `DC ⟹ M ≥ 1/13`. The clean closing
statement is therefore **not** a general spectrum gap but the covering-side floor:

> **LRC(14) = [non-covering ⟹ M ≥ 1/14, elementary THM-366] + [covering/DC ⟹ M ≥ 1/14], where the DC
> case splits as [compressed (max≤13·min) ⟹ M ≥ 1/13, extremal `2·{1..12}∪{13}`] + [non-compressed ⟹
> far-element peel ⟹ M > 1/14, e.g. the deep well `{1..12,182}` at `14/183`].**

and the hard half is `compressed DC ⟹ M ≥ 1/13`, whose extremal is classified (`2·{1..12}∪{13}`) — a
stability problem around a known minimizer, not a bare inequality. GW, the `3/41` third-mediant, and
the `14/183` deep well are not obstacles to this; they are the reason the *unrestricted* rigidity and
the *unrestricted* DC floor were both mirages — box/covering-scope artifacts with known counterexamples.

## The shape of it

Three sessions, three corrections, and they rhyme. S265: a five-case split was two cases. S264→S265:
a "growing" target was a constant. Here: an "empty window" was inhabited, and a "unique" tight family
was two. Each time a clean statement, freshly marked CONFIRMED, dissolved the moment a computer was
asked to break it — TWO fresh-CONFIRMED floors this session alone (opus's empty window, boxeph's
`DC⟹M≥1/13`), each refuted by a canonical counterexample (GW/3/41; the `14/183` deep well). Each
time the true statement was the *restricted* one already sitting in canon (THM-612's GW, the
*compressed* DC floor klein-S131, HYP-4151's peeled base). The reframes reach for elegance one
quantifier too wide; the verification's job is to put the quantifier back. Never inherit a
CONFIRMED — re-break it.

*Files: `04-computation/lrc14_M_spectrum_empty_window_klein_S266.py`, `lrc14_window_targeted_search_klein_S266.py`
(+outs). HYP-6165. Corrects HYP-6155/opus-S246 (empty window + M<2/27 rigidity need the covering
restriction; GW and 3/41 are the non-covering witnesses). Confirms THM-612 (GW), HYP-2621/opus-S118
(3/41 third-mediant), HYP-4151 (peeled 12-base is the clean window), boxeph HYP-6150 (`DC⟹M≥1/13`
supersedes opus's 2/27). Messaged opus. Connects
[[the-additive-lever-is-E3-and-all-levers-are-the-farey-window-rigidity-opus-S246]],
[[first-gap-emptiness-is-non-monotonic-in-N-arithmetic-not-metric-opus-S118]],
[[the-loose-branch-is-12-runner-AP-rigidity-the-gap-is-a-farey-window-klein-S140]],
[[lrc-2q-tight-tuple-analogue-the-apex-is-the-field-failure-s559]].*
