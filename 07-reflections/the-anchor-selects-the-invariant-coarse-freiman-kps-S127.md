# The anchor selects the invariant — the seven-sector residue is coarse Freiman

*kind-pasteur-2026-07-11-S127. Owner (to several machines at once): "attack the joint Φ-consec-extremality
lemma, see the bigger picture of what we're approaching." opus-S221 took the variance route the same hour;
this is the complementary half — which additive invariant governs, and why the answer flips from the one
opus proved for the density floor.*

---

## Where the whole problem now sits

Every thread of LRC(14) — detuned dispatch, the measure floor, the seven-sector cap, the wide-spread
recursion — has funnelled into one sentence: **the arithmetic progression `consec = {0,1,…,m−1}` maximizes
`Φ = p0 + (1/3)p1` over bounded integer cores.** mac-mini (THM-703) and opus (S220) reduced `Φ` to a
factorial-moment ladder of the empty-sector count `N`; the binding piece is that consec **minimizes**
`F2 = 4m1 − m2 = E[N(5−N)]`. opus-S221 rewrote this cleanly as consec **maximizing** the coverage-variance
`B = (E[N] − 5/2)² + Var(N)` — the AP wins by *bimodal resonance*: its circle-orbit covers perfectly on
resonance and clusters off it, so low mean miss-count and high variance at once.

That is the *why*. The question I took is the *what-measures-it*: which additive invariant of the core is
`F2`, so that "consec is extremal" becomes a known theorem?

## Two candidates, and a warning from our own files

There are two natural additive invariants of a set, and this project has already had a fight about them:

- **E2 = additive energy** `= #{a+b=c+d}` — the AP uniquely maximizes it (Freiman's `|S+S| ≥ 2n−1`,
  equality iff AP; opus **HYP-5681, proven**). mac-mini reaches for "additive-energy bricks" to prove the
  degree-2 extremality.
- **E3 = Schur triples** `= #{a+b=c}`. opus **HYP-5683** proved that for the LRC *density floor*, the
  governing invariant is **E3, not E2** — because loneliness is scale- but **not** translation-invariant, and
  only E3 shares that symmetry group. The clean separator: the translated AP `{1,3,…,25}` has *maximal* E2
  yet is loose (E3 = 0). E2's translation-blindness makes it fatal there.

So which is it here? `Φ` is scale- but not translation-invariant (THM-536) — the *same* symmetry as
loneliness. A naive import of HYP-5683 says: E3 should govern, and mac-mini is aiming at the wrong invariant.

## The measurement — and it flips

It does not. Over a battery of 0-anchored 8-cores (`lrc14_moment_invariant_symmetry_kps_S127.py`):

> `corr(F2, E2) = −0.978` vs `corr(F2, E3) = −0.632`; at **fixed spread** (max = 15) still `−0.941`.

`F2 ≈ 6.46 − 0.0106·E2`. **E2 — additive energy — governs the seven-sector residue, decisively, reversing
the density-floor verdict.** Why does the symmetry objection that killed E2 downstairs not apply upstairs?

**Because the seven-sector cores are 0-anchored.** The offset set always contains `0` (it auto-covers
sector 0; `0` is the top-speed offset). E2 is translation-invariant, so downstairs — where the speed set is a
free translation class — it cannot tell the tight AP from its loose translate, and fails. But the
seven-sector family contains **no translates of each other**: every core is pinned at `0`, and the only
identification within the family is scaling, under which E2 and `Φ` agree (verified: `consec` and its dilate
`{0,2,…,14}` give identical E2 = 344 and Φ = 0.4086). The anchor breaks translation, so E2's one blind spot
is never tested — and the simpler classical invariant is free to govern.

**The anchor selects the invariant.** Free set → E3 (Schur), the translation-sensitive one. 0-anchored set →
E2 (additive energy), the classical one. Same symmetry-match principle (opus's), opposite answer, because
pinning `0` changed the symmetry group of the family.

## The two halves are one object

opus's coverage-variance and my additive energy are the same quantity seen from two sides. `F2 = 25/4 − B`
(opus) and `F2 ≈ affine(−E2)` (here), so **`B ≈ affine(E2)`**: the coverage-variance *is* the additive energy
at leading order. opus's "bimodal resonance," his "AP rank-1 relation lattice," and Freiman's minimal-sumset
set are three names for one fact — the AP is the maximally self-resonant set, and that resonance is exactly
high E2. The variance route explains the mechanism; the energy route names the invariant and hands it to a
theorem.

Because that theorem — **AP uniquely maximizes E2 (HYP-5681, proven)** — the *leading term* of the
extremality is already settled. What remains open is only the **fidelity** `F2 ↔ E2`: the fit is `−0.98`,
not `−1`, and two cores with equal E2 (and equal E3) can still split `F2` slightly. That residual is the
sub-leading, degree-3 correction — the same triple-correlation `m3` the ladder needs for `k = 8`. So the open
problem sharpens to: **`F2 ≥ (affine in E2) − (small triple-correlation defect)`**, then invoke Freiman. That
is mac-mini's additive-energy lane, now confirmed to be aimed at the right invariant, with the target made
precise.

## The bigger picture — LRC at two scales, and why the coarse one is winnable

Here is what we are actually approaching. The Lonely Runner is an extremal statement about the arithmetic
progression, and it lives at two scales:

- **Fine scale (`1/14`, the real conjecture).** The tight locus is governed by **Schur triples E3**, because
  the speed set is a free translation class (opus HYP-5683). Translation-sensitive, and *open*.
- **Coarse scale (`1/7`, the seven-sector reduction).** The residue is governed by **additive energy E2**,
  because the offset set is **0-anchored**. Translation-neutralized — and E2's AP-extremality is **classical
  and proven** (Freiman).

The seven-sector reduction's quiet genius is exactly this: by pinning `0` (the maximum offset), it converts
the hard, translation-sensitive Schur-triple invariant into the easy, classical additive-energy one. We are
proving LRC(14) by dropping to a coarse scale on which its additive-combinatorial heart — "the AP is the
extremal set" — becomes a *solved* problem (Freiman), and then paying only a controlled triple-correlation
tax to climb back. Both scales are the same phenomenon (AP-extremality of a low-degree correlation); the
anchor is the change of variables that makes the coarse one a theorem instead of a conjecture.

And the "irreducibly aggregate" character everyone keeps hitting (THM-536's refuted local moves; opus's and
mac-mini's global-not-local) is no longer a mystery: E2 is a global four-fold correlation. You cannot see
additive energy one term at a time, so no local move or per-block argument can see the extremality. It was
never going to yield to compression — it was always going to yield to Freiman.

## The meta-lesson

The symmetry-match principle — *the governing invariant is the one that shares the problem's symmetry group* —
is a machine, and opus built it. Using it correctly means reading the symmetry of **this** problem, not
importing a sibling's verdict. Downstairs the group excluded translation, so E3. Upstairs the anchor already
broke translation, so the group is only scaling, and E2 — blind exactly where it no longer needs to see — is
free to govern. Same machine, opposite output, because the anchor moved the group. When two problems look
alike, check what their symmetries actually are before you carry an answer across; the anchor is small and it
decides everything.

*Files: `04-computation/lrc14_moment_invariant_symmetry_kps_S127.py` (+`.out`). HYP-5990. Complements
opus-S221 (coverage-variance), mac-mini THM-703 (moment ladder), opus HYP-5681/5683 (Freiman / the E3
density-floor verdict this one flips under the anchor). Extends
[[the-certificate-that-doesnt-factor-kps-S127]].*
