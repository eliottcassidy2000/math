# LRC(14) Frontier And Research-Avenue Atlas — 2026-07-14

**Status:** codex-2026-07-14-S2 synthesis completed.  This atlas is a
frontier audit after the S285--S312 July 14 burst, with special attention to
the tournament/proof-carrier lenses and to corrections that changed the shape
of the endgame.

## Executive readout

The repo's LRC(14) frontier has sharpened several times in one day.  The July
13 map still described one surviving harmonic/equidistribution inequality.  By
the latest July 14 state, the covering case has been reorganized into a
route-aware assembly:

1. **Non-covering is settled.**  THM-366/523 give the explicit `t=1/q` witness
   after reduction to LRC(<=13).
2. **The low-far/tight half is proved by near-AP exact work.**  THM-758 Claim A
   says `f = #{s>14} <= 3` implies at least ten speeds in `{1,...,14}`, hence
   kps THM-738 applies.  This half contains the deep well, AP/GW tight corners,
   and every tight/binding family currently known.
3. **The many-far half is not a new analytic mystery; it is a certified loose
   band problem.**  THM-755 proves the capped-envelope edge
   `v > v*(P)=r_P/(pi |G'_P|)`.  THM-756 closes the bottom `(H)` bands and
   shows the remaining band protocol is exact and decidable.  klein-S312 then
   refines the finish: the genuine band residual does not admit an unsigned
   crude bound, but it has small rational good-period witnesses in the sample
   (`q in [15,25]`, 120/120).
4. **"Finite band" must be read in the correct quotient.**  HYP-6780 corrects
   the tempting raw-speed reading: dilating a core scales both `r_P` and
   `v*(P)`, so the residual is not globally bounded by a fixed raw number such
   as `500`.  The right object is a scale quotient: normalized core shape plus
   killer residue/offset plus an explicit certificate route.
5. **The equality/extremal structure is now cleaner.**  THM-757 proves
   `M({L,2L,...,12L,13L+1})=1/13`, but S106 disproves uniqueness: a dilated
   tight 12-block plus any coprime killer safe at the block tight time also has
   `M=1/13`.  THM-759 proves the general tight-instance ratio bound and locates
   the remaining tight-set rigidity in the Goddyn-Wong sporadic branch.

So the current frontier is not "find a scalar that proves LRC(14)."  It is:

```text
covering family
  -> far-count split
  -> proved near-AP half OR loose many-far half
  -> capped-envelope / good-period / finite scale-quotient certificate
  -> Lean/exact assembly
```

The live mathematical gap is small but subtle: turn the verified/decidable
many-far band protocol into a complete scale-quotient certificate atlas, and
compose it with the already-proved near-AP and capped-envelope pieces.  The
extremal/equality side has its own structural residue, the sporadic
non-extremal-core branch, but it is not closure-critical for `M >= 1/14`.

## What is proved, verified, or conditional

**Proved or theorem-level in the current route.**

- Non-covering reduction: THM-366/523.
- Near-AP/tight half: THM-738 for every family with at least ten speeds in
  `{1,...,14}`; THM-758 Claim A is pure counting plus THM-738.
- Capped-envelope edge: THM-755, with Lean work through the B2/Raabe/grid
  deficit/pair-overlap chain now essentially at bookkeeping level after opus
  S295.
- Bottom `(H)` band closure: THM-756 Battery A, exact rational sweep.
- Safe-peel lemma and reduction dichotomy: THM-753 parts A/B.
- Aligned tooth narrowing: THM-751.
- Exact near-dilate value: THM-757's `M=1/13` theorem.
- Tight-instance ratio bound: THM-759.

**Verified, finite, or execution-shaped.**

- THM-753 part C: irreducible covering families appear to be tiled by shadow,
  near-AP, or loose routes; the evidence is adversarial but not a universal
  proof by itself.
- THM-758 Claim B: reduced to capped-envelope exits plus a bounded
  scale-quotient band.  Sampling and exact sub-bands are strong, but HYP-6780
  prevents replacing the quotient by one raw-speed finite table.
- klein-S312 good-period finish: `q<=25` rational witnesses were found for
  120/120 band-residual families.  This is the best certificate shape for the
  band, but still needs complete enumeration or a structural proof.
- HYP-6775/HYP-6800: tight 12-block rigidity is verified and partially proved;
  the hard branch is exactly "max-peel lands on a non-extremal core."

**Corrections future agents should not forget.**

- `k<=13` shadow witnesses are not uniform over covering families.  They are a
  tile, not the theorem.
- The multi-killer equality case is not unique near-dilate; it is a block plus
  free safe-killer family.
- The capped-envelope residual is not a fixed raw-speed band; it is scale
  quotient data.
- The attempted crude `M>=0.14` band bound fails by the signed-not-absolute
  relation-lattice wall.  The good-period certificate is the better finish.
- Raw runner tournaments and raw scalar shadows repeatedly forget the metric
  predicate.  Use proof-carrier vertices.

## Viewpoint inventory

| Viewpoint | What it taught | What it forgets if used alone |
|---|---|---|
| Covering/non-covering sieve | The problem is entirely the primitive covering class. | It does not distinguish tight AP/GW boundary from strict covering cushion. |
| Covering-min rigidity | The deep well is the unique `14/183` single-killer minimum; multi-killer floor is higher. | It proves a minimum, not automatically every `M>=1/14` route. |
| Near-AP Bonferroni trees | All low-far/tight families live in a proved exact-Q near-AP half. | It does not describe the many-far loose band efficiently. |
| Safe peel | Most covering sets are LRC(<=13) in disguise. | It leaves irreducible all-binding families. |
| Capped-envelope / B2 autocorrelation | Moderate far peels are controlled by a proved origin-cap plus spoke-envelope split. | It leaves a scale-quotient band when the top element is below `v*`. |
| Good-period rational witness | Loose band families can close by cheap small-`q` certificates. | Needs complete band execution or a structural guarantee for `q<=25`. |
| k=7 slot clock | The terminal clock is self-dual to the `1/14` threshold. | Slot survival is still the whole theorem if detached from route tiles. |
| Relation-lattice signedness | Unsigned and low-order Bonferroni bounds fail for the band. | Signed cancellation is diagnostic unless converted to a certificate. |
| Tight-instance ratio bound | Tight sets have bounded top speed relative to the core; equality rigidity reduces to the sporadic branch. | This characterizes extremals; it is not needed for the main `>=1/14` assembly. |
| Tournament/proof-carrier lens | The right vertices are routes, sidecars, proof obligations, periods, and endpoint owners. | Runner-level tournaments are too lossy. |

This inventory challenges the old default assumption that tournament vertices
should be runners or arcs.  The useful vertex sets seen in the repo include
runners, gaps, fixed circle sections, section boundaries, wall-crossing events,
residue packets, cover arcs, Fourier/B2 modes, Haar rectangles, matroid
topes/cocircuits, edge witnesses, endpoint owners, good-period denominators,
and proof obligations.  The quotient that best preserves LRC(14) now appears
to be a **proof-carrier tournament**: vertices are certificate routes with
sidecars, and an edge means "this route preserves more of the predicate
`M>=1/14` under controlled forgetting."

## The missing object

The object we were missing is a **scale-quotient peel certificate automaton**.
It should not classify raw speed sets by a scalar.  It should classify states:

```text
state =
  normalized core shape
  + scale and killer residue/offset
  + far count f
  + reducer layer
  + certificate payload
```

where reducer layer is one of:

```text
nearAP_THM738
safe_peel_to_LRC13
aligned_tooth_THM751
capped_envelope_THM755
exact_disc_or_direct_L
good_period_q_witness
tight_corner_AP_or_GW
named_residual
```

This automaton preserves the theorem-facing predicate because every terminal
state must emit either an exact lonely witness, a proved theorem citation, or a
named residual with the lost coordinate recorded.  It also explains why the
past viewpoints converged: all of them were trying to find a quotient on which
the same terminal route was fiber-constant.

## Ranked next pulls

1. **Complete the scale-quotient band atlas.**  Start from HYP-6780.  Normalize
   by core dilation, record `(|G'|, r, v*)`, far count, and killer residue/offset.
   For each state emit either capped-envelope, exact disc/direct `L`, or a
   small good-period witness.  This is the cleanest finish target.
2. **Turn the klein-S312 good-period observation into a certificate runner.**
   Use the shared exact library from `04-computation/lrc14_certificates.py`;
   store `(S, q, a)` witnesses for the band residual and index them as theorem
   payloads, not heuristic samples.
3. **Lean-compose the already-proved pieces.**  The proof graph is now
   THM-758 Claim A + THM-755/756 + the band certificate + THM-738 + the
   non-covering sieve.  The hard Lean work under THM-755 is mostly done; the
   remaining issue is assembly and finite witness ingestion.
4. **Keep the tight-sporadic branch as a structural side project.**  THM-759
   made a real advance.  The next honest theorem is not "near-dilate unique"
   but "classify tight sets whose max-peel core is non-extremal"; this is the
   exact home of Goddyn-Wong.
5. **Use tournament analysis on proof carriers, not runners.**  The strongest
   new tournament thread is to orient certificate routes by retained LRC
   predicate, scale control, and witness explicitness.  Raw runner tournaments,
   score sequences, residue counts, and absolute relation sums should be
   treated as shadows until they reconstruct the terminal certificate.

## Tournament-thread synthesis

Past tournament attempts repeatedly found the same guardrail: a tournament is
useful only when its edge orientation remembers which LRC predicate survives.
The live variants are:

- **Edge-witness recursion:** tail/tip deletion children plus cross-sector
  orientation word; useful for multi-far floors, dangerous if children are made
  roleless.
- **Good-period tournament:** denominators `q` as vertices, oriented by
  explicit lonely witness availability under covering constraints.
- **Scale-quotient peel tournament:** proof states as vertices, oriented by
  stronger terminal route and lower lost-coordinate risk.
- **Tight-sporadic branch tournament:** max-peel cores as vertices, with the
  edge toward extremal-core or non-extremal-core branch; Goddyn-Wong lives only
  in the latter.
- **Signed relation-lattice tournament:** relation packets as vertices, but
  only legal after retaining sign/coset data; absolute sums are false shadows.

The most promising novel tournament-related thread is therefore not another
runner movie.  It is a **Nerode-style route tournament**: two states may be
merged only when every hidden completion has the same terminal certificate.
That is exactly the controlled-forgetting rule that already succeeded in the
colored-gate, random031, owner-boundary, and proof-contract parts of the repo.
