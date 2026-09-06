# Located overlap removes the displayed 16,704 edge-and-word boundary

**Status: FINITE-EXACT complete one-pair certificate + PROVED conditional translated-grid transfer; INDEPENDENTLY AUDITED.** This is an appendix to `continuing5_20260906_lrc_correlated_overlap.md`. It eliminates the explicit selected edge and full-word control used at the largest incoming relaxed clock. It does **not** eliminate every candidate at clock16,704 or construct an unsafe or safe complete physical realization of that word.

The subsequent [complete one-clock theorem](continuing5_20260906_lrc_clock16704_complete.md)
now excludes all cases at16704. This report retains the intermediate
mechanism and its declared finite universe; it is not the current scale bound.

## Exact overlap with the incoming result

The independently audited `third_20260906_grid_refined.md` retains8,202 necessary clocks, with largest16,704. Its displayed boundary has

```
t=16704, e=4, (p,q)=(3,308),
forced margins=(12,16),
d=(12,16,72,58,64,9,9), E(d)=188.
```

The pair is a strict actual atlas pair: gcd(3,308)=1 and sum311 is a prime congruent2 modulo3. The word passes all126 projected proper-subset complement conditions, all seven states divide t, and their gcd is1. Its separate interval credit is172, so the incoming sufficient test correctly retains it. That result is not retracted or contradicted; the location coordinate supplies a stronger lower bound.

Two different losses must be distinguished. The incoming full-word maximum M(t) forgets the actual seven ratios and their simultaneous phase incidence; its separate E_pair bound additionally need not have the same maximizing owner as M(t). Our calculation does not repair either entire loss. Instead, for a **fixed retained pair**, the incoming interval credit sums independently minimized counts. It discards that all interval counts depend on the same translate. Their relative locations are already determined by the retained p,q, but were not used in the consumer. Thus this is a stronger use of an available pair coordinate, not a claimed new source of actual seven-shape information.

## General located-pair credit

For coprime p<q define the open danger intersection

`I_(p,q)={x mod1: ||px||<1/14 and ||qx||<1/14}`.

Write its circle components as open intervals(a,b), using one unwrapped representative for the origin component. For n>=1 define

```
L(n,p,q)=min_(alpha mod1) sum_(a,b)
          [ceil(n b-alpha)-floor(n a-alpha)-1].
```

This is the exact minimum number of points of the translated n-grid X_alpha={(alpha+j)/n mod1} in the full intersection. Endpoints are excluded. Its only discontinuities are the finitely many residues n*a,n*b modulo1; evaluate each wall and one point in every intervening open cell. This proves the minimum for all real translates.

If an actual seven-component pair is(Dp,Dq), let e=gcd(t,D), n=t/e. Multiplication by D takes a translated t-grid to a translated n-grid, with each point repeated e times: gcd(D/e,n)=1. The physical common scale g is also harmless when gcd(g,t)=1. Thus the exact uniform pair credit is

`C_loc(t,e,p,q)=e L(t/e,p,q)`.

It always dominates the separate credit

`C_sep=e sum_(a,b) [ceil(n(b-a))-1]`,

because a minimum of a sum is at least the sum of independent minima. This inequality may be strict: the minima need not be simultaneously attainable. Open endpoints cannot be rounded away.

If the seven danger counts have total ceiling excess E, the inherited pointwise overlap inequality gives at least C_loc-E weak-safe grid points. The selected pair indicator is at most the excess multiplicity(m-1)_+ at every grid point. A strict positive integer lower bound therefore closes the conditional actual row after lifting any six-body safe phase. The six-body phase supplier and general CRT consumer are inherited from `third_20260906_grid.md` and its cited source; they are not re-proved as a new gluing theorem.

## Complete computation for the displayed boundary

Here n=4176. Literal rational intersections give45 circle components:43 have length1/2156 and two have length1/4312, so their total measure is1/49. Consequently

`sum separate minima=43`, and `e*43=172`.

There are exactly90 distinct endpoint walls modulo1. The source evaluates all90 walls and all90 intervening cell midpoints. Every value is independently checked by the literal strict predicates at all4,176 points of the translated grid, without consulting the interval representation. The exact extrema are

```
L(4176,3,308)=84, attained at alpha=53/539;
maximum count=87, attained at alpha=74/77.
```

Thus

`C_loc=4*84=336>188`,

and the inherited overlap inequality leaves at least148 weak-safe lifts for any actual row with this selected edge and total ceiling excess at most188. In particular it defeats the displayed full word. The incoming independently audited full-word maximum M(16704)=188 extends the same conditional conclusion to every allowed word attached to this selected edge; that numerical maximum is an inherited theorem, not proved merely by checking the displayed owner here.

The finite supplier is exact for all translates, including event coincidences. This statement is about weak safety. No strict clearance is inferred just from the weak grid endpoints.

## What remains, and the cheapest next decisive test

The earlier actual t968 control had separate minimum0 and located minimum15; it is independently safe at denominator7. It establishes the method loss in a complete actual entry, not a residual failure. The present one-pair calculation shows that the same coordinate also repairs a specifically retained extremal quotient boundary.

To remove clock16,704 itself, one must address **every** candidate selected edge supplied by the small-sheet-gcd theorem: each divisor e<=6 of t and every surviving strict atlas pair with its correct forced margins. Excluding the first successful or first enumerated pair is insufficient, exactly as the incoming conditioned-word hostile already showed. The next bounded theorem to test is:

> At t=16704, every full-word-compatible candidate(e,p,q) has C_loc(t,e,p,q)>min(M(t),E_pair(t,e,p,q)).

A positive exact exhaustion would remove that clock by the already-proved physical transfer. A failed candidate should be retained with its located profile and compatible maximizing word, rather than triggering a blind census of the remaining8,201 clocks. This appendix does not run that exhaustion or claim its outcome. It also does not recover actual seven-shape realization, the other physical profiles, the common translate of all seven masks, or the restricted set of translates supplied by V-safe phases.

## Reproduction and scope pins

```
python ../../04-computation/continuing5_20260906_lrc_extremal_location.py
python -O ../../04-computation/continuing5_20260906_lrc_extremal_location.py
```

Both runs pass320 always-active exact gates and produce byte-identical actual LF output. The source configures LF and capture rejects carriage returns. The full projected-word data pin is `935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f`, from `04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json`.

Source SHA256:`6526dd4306c31861792c5e938b37e498675cb1f39f9a05b8158d2f0940d44659`.
Output SHA256:`954c713447f9ff4ba162f0318195d12b5a983c6e646ca54821864a661a625ed8`.

No producer import, repository edit or Git operation was used. There is no complete physical realization claim for the displayed word and no whole-clock elimination claim.

Referee portability repair: the source now prefers an adjacent inherited profile JSON when filed, with the external-worktree fallback unchanged. Only path resolution changed; all320 gates and exact output bytes are unchanged. Original source pin was261136ff50f1381f8f57825a313728b4504d509ae76e98d02bb28da7c575af20.

Independent [proof and exact referee](continuing5_20260906_lrc_extremal_location_audit.md) passes.

The [raw-byte checkpoint manifest](continuing5_20260906_manifest.json) pins
the filed source, report and identical normal/optimized transcript. Any
candidate-report hashes above identify the pre-promotion bytes.
