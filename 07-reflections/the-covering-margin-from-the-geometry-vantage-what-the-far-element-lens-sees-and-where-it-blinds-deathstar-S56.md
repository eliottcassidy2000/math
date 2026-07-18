# The covering margin from the geometry vantage: what the far-element lens sees, and where it blinds

**death-star-2026-07-17-S56.** A long-session exploration of the LRC(14) frontier from the
covering-geometry vantage I built this session (Lemma A bounded denominators, Lemma B good-component
widths, Theorem T / THM-1002 the quantitative top-gap margin). Purpose: push the geometric lens as far
as it goes on the *margin floor* half of the nucleus, map exactly what it resolves and where it blinds,
and record the equivalences that pin the residual. This is exploration + one new theorem (THM-1002),
not a claim on the hard core.

---

## 1. The frontier, as it now stands (post THM-995 IX/X)

The reduction bottoms on one object, and it has been sharpened to a clean testable form:

- **Sieve-margin lemma (PROVED, THM-995 IX):** if any `q ∈ {2..13}` divides no speed, `M ≥ 1/q > 1/14`
  at `t = 1/q`. So a tight family covers `2..13` and must **miss 14**.
- **The residual:** *does any primitive family covering all of `2..14` have `M = 1/14`?* Conjectured no
  (empirical covering-floor, this session's strong descent: `M ≥ 1/11`, at low-denominator resonances).

**The residual is exactly the census "no mult of 14" question.** Covering `14` means a multiple of 14;
a non-covering mult-14 family already has a sieve margin; so

> primitive fully-covering ⟹ `M > 1/14`  ⟺  primitive tight ⟹ **no multiple of 14**.

And this is implied by the THM-997 no-ghost residual: a multiple of 14 kills every base resonance
`a/14` (it sits at `0` there), so a tight mult-14 family would be lonely only at denominator `≥ 28`
(Lemma A) — which "primitive ⟹ no denominator-`>14` loneliness" forbids. The chain
**covering residual ⟸ no-mult-14 ⟸ THM-997 no-ghost residual** ties the whole session's threads to one
wall. `2·AP` (tight, covering, has 14) shows the escape hatch is exactly **non-primitivity**.

## 2. What the geometry lens resolves: the far-element regime, quantitatively

**THM-1002:** `M(V) ≥ Vmax · M(V∖{Vmax}) / (Vmax + v₂)`. This is Lemma B (a good-component of the
sub-family, wider than one arc of the removed speed, forces loneliness) plus a Lipschitz lower bound on
the sub-family's optimal component. It is *sharp* where a top element dominates:

- as `Vmax → ∞` the bound `→ M(V∖{Vmax})`, and since `M(V) ≤ M(V∖{Vmax})` always, **`M(V) = M(V∖{Vmax})`**
  — a large speed is invisible to `M` (the dominance-dodge, now with an explicit threshold);
- `{1..12,500}` sits at `M = 1/13` with the bound within 2%; the deep well likewise.

So the *entire* far-element / dominant-outlier stratum of the covering problem is closed by an
elementary, explicit inequality — no enumeration, no harmonic analysis. This is the geometric form of
THM-721/755/758 and of the covering-min rigidity's "far element can't be too far."

## 3. Where it blinds: the density argument and the alignment wall

The same lens shows *why* the clustered case is hard, concretely. Peel the multiple of 14, `w = 14k`,
and ask whether it can pull a covering family down to `M = 1/14`:

- **If `14k` is large**, its danger arcs are fine and evenly spaced — density `2/13` inside *any*
  good-component of the rest — so they leave an `11/13`-fraction uncovered: `M ≥ 1/13 > 1/14`. A large
  multiple of 14 **cannot** produce tightness.
- **If `14k` is small (`= 14`)**, its arcs are coarse (`width 1/91`, spacing `1/14`) and *can* align to
  cover the good-components of the other twelve — but only if those components sit exactly at the
  `j/14` lattice. That alignment is precisely the `2·AP` dilation (`{2,…,26}` is lonely at `1/28`
  because it *is* the AP scaled). **Primitivity forbids the dilation**, and whether a genuinely
  primitive family can align its twelve-speed good-set to the coarse `14`-lattice **is the apex-7
  wall** (7 danger sets of measure `1/7` tiling, union bound collapsing at `k=7`; §boundary reflection).

So the geometry lens converts the residual into a crisp dichotomy: *tightness with a multiple of 14
requires the coarse `14`-arcs to resonantly tile the rest's good-set — a dilation phenomenon that
primitivity is exactly designed to exclude.* The lens sees the mechanism perfectly; it cannot decide
the one Diophantine alignment question, because that question is LRC(14).

## 4. The margin landscape (empirical, this session)

Minimizing `M` over covering families (strong multi-start descent, `Vmax ≤ 40`) bottoms at `M = 1/11`
on the primitive family `{2,5,6,8,11,13,14,16,19,20,24,25,27}` at witness `t = 1/22` (denominator
`2·11`). The minimizers are **low-denominator mod-`p` resonances** (`t = a/(2p)`), family-specific —
which is why the `1/7` "double-threshold" was withdrawn (it was the `p=7` instance). The honest reading:
the covering margin is *substantial* (`≳ 1/11`, ~1.3–1.6× threshold) and *resonance-structured*, not a
razor-thin `1/14 + ε`. Whatever proof eventually closes the residual should explain this margin
landscape, not just bare positivity.

## 5. What remains, and the two live faces

The clustered/aligned residual is the apex-7 nucleus = LRC(14). Two faces, per the finish map:

- **Route B (covering):** the alignment question above — sharp, tight, needs the multi-linear/E₃
  cancellation. The geometry lens is exhausted here (it caps at THM-1002, `≈ 1/(2(n−1)) < 1/n` for
  clustered families).
- **Route A (density):** the *softer* face — `Q_s = o(r²)`, any power-saving suffices. The margin
  landscape of §4 lives here (it is the singular-series / good-measure `μ₀`), and the resonance
  structure of the minimizers suggests the right coordinate is the family's Diophantine position against
  the small-denominator lattice — exactly the `w`-grid discrepancy `Q_s` (THM-729).

**Unbuilt tournament bridge worth the margin (D5):** "not the regular tournament ⟹ margin bounded
below." The AP is the marginless regular/tight case; every primitive covering family is "not regular"
and carries the `≳ 1/11` margin. The geometry lens now says the margin is a good-component-escape
measure; the tournament lens says it is a Dedekind-sum/cusp-form (`disc_v`, THM-732). Making those two
the same object — the margin as a discrepancy that vanishes *iff* the AP-lattice alignment holds — is
the concrete next synthesis, and it is where the density face and the tournament face meet.

**Empirical test of D5 (this session) — the naive form breaks on GW.** Computing the winding
tournament `T(t*)` (14 vertices incl. the observer; `i→j` iff `frac((s_i−s_j)t*) ∈ (0,½)`) at each
family's optimal time:

| family | `M` | score sequence | balanced? |
|---|---|---|---|
| **AP** | `1/14` | all `6` (7 antipodal pairs tied) | **yes — perfectly regular** |
| GW | `1/14` | `5,5,6…6,7` | no (range 2) |
| deep well | `14/183` | `6…6,7…7` | near (range 1) |
| covering-min | `1/11` | `4…8` | no (range 4) |
| `{2..14}` | `1/8` | `5…7` | no (range 2) |

The AP is the **unique** family whose optimal-time tournament is perfectly balanced — the regular
"heptagon" tournament with the `n/2 = 7` antipodal (distance-`7`) pairs tied (the `z⁷=−1` structure at
`t=1/14`). But **GW is tight yet irregular**, and the deep well is near-regular yet *not* tight: so
"regular ⟺ tight/marginless" is **false**. The correct refinement: tournament-regularity singles out
the **AP specifically** (the difference-closed tight family), not the tight locus. This matches THM-996
§III (the census is flat across `{AP,GW}` at the threshold) from the tournament side — the balanced
tournament is another threshold-level invariant blind to the AP↔GW distinction. So the D5 margin bound,
if it exists, must key on a *sub-threshold* tournament statistic (the `disc_v` Dedekind sum, not the
score balance), consistent with the whole session's finding that rigidity lives below the threshold.

→ THM-1002, THM-1000, THM-999, THM-997, THM-995 (IX/X), THM-729/732, HYP-7305; boundary reflection
`bounded-spread-classification-is-the-apex-7-wall`.
