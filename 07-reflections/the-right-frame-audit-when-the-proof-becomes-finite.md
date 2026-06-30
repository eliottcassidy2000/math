# The right-frame audit: when the proof becomes finite, and the beautiful mechanism that wasn't

*klein-2026-06-29-S9. The owner asked to find the right frames for the other open LRC problems and to understand what we are missing, given all the connections. The audit, and the honest answer.*

## The principle, sharpened

A "right frame" is the one in which the hard inequality becomes a **finite check**. The wrong frames are
the ones that keep it analytic (per-set bounds, asymptotic estimates, Cauchy-Schwarz, averaging). The
floor's reframe (HYP-3566/3571) was the first instance; pushing it one step (HYP-3581) revealed the rule:
**the LRC proof's hard quantities live on finite mod-`7` (mod-`14`) data, so in the right frame they are
finite cyclotomic computations, not estimates.**

## The audit

| open problem | wrong frame (analytic, stuck) | RIGHT frame (finite, closed) |
|---|---|---|
| **`rho_j >= c`** (descent decorrelation) | per-level Cauchy-Schwarz; `Z_7\*`-averaging (INVALID, overshoots 30/127, HYP-3581) | **finite min** over the `2^7` cores `O subset Z_7` of the cyclotomic Gram gap = `4cos^2(3pi/7) = 0.198`, binding at the **doublet** (THM-578); covering boundary at `O=Z_7` |
| **`inf R' >= 1/(2 zeta(2))`** (global floor) | the per-set variance `CV(N_R)^2` (unbounded, set-dependent, HYP-3554) | the actual correlation `|R'-1|` in the doublet/`Gamma_0(14)` local frame; `inf R' = 0.344` (HYP-3571), closed form pending from the `p=2,7` local densities |
| **`m_Q`** lower bound | the per-set pairwise cap `C(14-r,2)/91` (THM-576) | the `Q`-core mod-`14` density -- a finite/multiplicative quantity, the same kind of object as the `R`-core; should match the `phi/psi/J_2` arithmetic |
| **gap line `M >= 7/89`** | searching for a disproof "far in gap value" | **multiplicative**: `M = 1/(smallest surviving denominator)` (HYP-3550) -- a DIVISIBILITY frame, *different* from the cyclotomic-gap frame; do not conflate the two (a frame-polysemy caution) |

The first row is the one that closes: in the right frame `rho_j >= c` is not an estimate at all, it is
`min` over 128 subsets of `Z_7`, equal to a single algebraic number `4cos^2(3pi/7)`, attained at the
doublet. The "BOUNDED" clause of the proof sentence is a finite cyclotomic fact.

## What we were missing (three things)

**1. Finiteness.** The relevant objects are mod-`7` residue cores -- a *finite* set. We kept reaching for
analytic machinery (CV, Reynolds averaging, asymptotic SOS) for a quantity that is a minimum over `2^7`
cases. The right frame doesn't bound the hard inequality; it *dissolves* it into a finite computation. The
hardness was an artifact of the wrong (set-indexed, infinite) frame.

**2. A beautiful mechanism that wasn't valid.** `Z_7\*`-averaging (manufacture transitivity by Reynolds
projection) is elegant, fits every connection (transitivity, the metagraph reference-collapse, Kaczmarz),
and is the natural-sounding move. It is also *invalid* as a lower bound -- it overshoots the gap for 30 of
127 cores, because Jensen gives `<= mean`, not `<= min`. Only an exhaustive check caught it. The lesson is
the persistence/validity-test discipline (the same one that killed the `b_1^-` octonion, HYP-3563), now
applied to a *mechanism*, not a number: **a plausible move that respects all the connections can still be
the wrong move; verify it on the finite case set before believing it.** The valid mechanism was sitting
next to it -- the Fejer-Bochner minorant (= the finite min), which mac-mini's S27 already had.

**3. The connections are the map, not the territory.** The project's web -- metagraph = finite Siegel
transform, complement = antipodal, octonion apex, perfect numbers, the descent, `X_0(14)` cusps -- is a
*map* that told us where the proof lives: on finite mod-`7`/`14` cyclotomic data, on the `sigma`-even floor,
binding at the doublet. But most of those connections are `sigma`-odd / structural (the witness side), and
are *orthogonal* to the finite check that is the actual proof (HYP-3548 said this; HYP-3581 confirms it --
the floor closes by a `2^7`-case minimum that uses none of the octonion/perfect-number structure). We were
missing the willingness to stop collecting connections and *do the finite computation* in the frame the
connections pointed to.

## The proof, as it now reads

The covering floor is one sentence about the danger relation (HYP-3571): it does not factor (essential;
the extremal lonely set is the units `(Z/14)\*` in `phi/2` antipodal Borsuk-Ulam pairs), and composed with
itself it stays small (`rho_j >= 4cos^2(3pi/7) > 0`, a finite min over the `2^7` mod-`7` cores, binding at
the doublet; `inf R' >= 1/(2 zeta(2))`). Both clauses, in the right frame, are finite cyclotomic facts. The
remaining work is bookkeeping in that frame: confirm the descent lands in `Z_7`-cores (S28), fix the `rho_j`
normalization, and write `inf R'` in closed form from the `p=2,7` local densities. The analysis is gone;
what is left is finite.

See HYP-3581 (the finite-min floor + averaging-invalid), HYP-3571 (the proof sentence), HYP-3566 (the
transitivity reframe), HYP-3575/3576 (mac-mini's S27 minorant / S28 averaging), HYP-3554 (the CV trap),
[[polysemous-constants-bridges-traps-and-homonyms]] (the validity-test discipline), THM-578 (the doublet),
HYP-3550 (the gap-line multiplicative frame).
