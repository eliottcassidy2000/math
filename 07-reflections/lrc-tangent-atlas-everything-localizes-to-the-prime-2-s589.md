---
source: opus-2026-06-03-S589 (remote-control)
status: HISTORICAL CREATIVE TANGENT ATLAS — five out-of-box frames, with the former field-ramification and prime-2-only claims corrected; the rigorous residue is a nonunit/non-etale doubling obstruction, not a reduction of all LRC14 to prime 2
tags: [LRC, adelic, local-global, renormalization, hyperoperation, quasicrystal, mahler, covering-radius, prime-2, n14, tangents]
---

# A tangent atlas — and they all localize to the prime 2

**Prompt (user):** keep thinking like this, even more out of the box; explore tangents
abstractly and creatively; think about what you can use them for.

Five wild frames, each with a concrete use. Historically they were read as
putting the entire `n=14=2·7` obstruction at prime `2`. The corrected lesson is
narrower: several frames detect the same **doubling/nonunit face**, but they do
not prove that this face is the only LRC14 obstruction.

**2026-07-21 correction.** "Ramified place" is not literal here: `2` is
unramified in `Q(zeta_14)=Q(zeta_7)`. The rigorous prime-2 content is that
`2` is a nonunit modulo `14` and the characteristic-2 order-14 group algebra
is nonreduced. Read every later "ramification" phrase in this reflection as
that nonunit/non-etale obstruction; none of those phrases proves that the rest
of LRC14 is solved. THM-404/2041 record the exact correction.

## Tangent 1 — ADELIC / local–global (developed, grounded)

Read loneliness **place by place**. Under `exp` (S588) the archimedean place `∞` carries
the static **cyclotomic** rigidity; the prime `p` carries the `×p` (Frobenius) dynamics.
For `n=14` the clock `j=1..13` stratifies by `gcd(j,14)` (the **divisors of `n`**):

| stratum | clock points | role |
|---|---|---|
| `gcd=1` (primitive 14th roots, `∞`) | `{1,3,5,9,11,13}` | **the witnesses — lonely** |
| `gcd=2` (7th-root shadow) | `{2,4,6,8,10,12}` | not lonely (runner 7 hits 0) |
| `gcd=7` (apex `2q`) | `{7}` | not lonely |

- **`ℚ(ζ_14) = ℚ(ζ_7)`** (since `ζ_2=-1` is rational): the *field* is the **prime-7** case.
- **Multiplication by 2** (`x↦2x mod 14`) sends every witness *off* the witness
  set (`1→2, 3→6, …`) because `2` is a nonunit modulo `14`. This is not the
  arithmetic Frobenius at an unramified prime of `Q(zeta_14)`.
> **USE — a non-etale sidecar:** retain the characteristic-2 nilpotent filtration
> or move to a good unit prime such as `3`. This addresses the doubling face only;
> a safe/dual seed and a pointwise LRC exit remain separate global obligations.

## Tangent 2 — HEIGHT / covering-radius (a useful *negative*)

Tempting: is `M(S)` a Mahler measure (a *product* of phases)? **No** — `M = max_t min_i`
is a **sup/min** (an `L^∞` covering radius), not a product. The normalized product
`∏ n‖v/n‖` grows (`36, 8.6e4, 3.6e6` for `n=7,12,14`) — it is *not* the margin.
> **USE:** this rules a frame *in* and a frame *out*: LRC is **not** a Mahler/height
> (product) problem; it is a **covering-radius / successive-minima** problem — the
> geometry-of-numbers (Minkowski) lane of the literature. The right invariant is the
> *covering radius of the dual lattice*, and the worry-set is its deep hole.

## Tangent 3 — RENORMALIZATION (the doubling IS an RG)

The **doubling map `x↦2x` is a renormalization** (coarse-grain the circle by `2`). In RG
language:
- the **worry-set is the RG fixed point** (scale-invariant: `AP_n`'s even layer `=`
  `AP_{2n}` rescaled, S579 — *exact self-affinity* = an RG fixed point);
- **loose configs flow to the trivial (high-temperature) fixed point** under repeated
  doubling (they spread);
- **`n=2q` has a RELEVANT 2-direction:** the `×2`-fragmentation (S585) is a *relevant
  perturbation* at the fixed point — the RG flow does not contract it.
> **USE — an RG-flow proof of rigidity:** show the doubling RG contracts every config to
> either the trivial fixed point (loose) or the cyclotomic fixed point (`M=δ`), with the
> only relevant (non-contracting) direction at `2|n`. The worry-set = the RG critical
> manifold; the n=14 residual = the single relevant `2`-eigenvalue. (Connects the repo's
> antiferromagnetic-tournament/Ising thread to the doubling tower.)

## Tangent 4 — FRACTIONAL HYPEROPERATION (the hardness is one bump)

On the hyperoperation tower (`succ, +, ×, ^, ^^`), LRC difficulty is a **single bump at
`+`**: trivial below (succession = one runner) and above (`×, ^, ^^` = lacunary, S588).
The phase transition additive→lacunary sits at a *fractional* rung (a Behrend/Sidon
density threshold).
> **USE:** *certify that the proof is purely additive.* No multiplicative/exponential
> structure is needed — the whole conjecture is additive combinatorics on the densest
> (AP) rung (the 3-term-relation / augmentation work, THM-400). The tower says: don't
> look up, look at addition.

## Tangent 5 — QUASICRYSTAL / aperiodic order (wild)

The runner positions `{v_i t}` are a point set; the **worry-set (AP) is the periodic
crystal** — minimal complexity, `ζ_n` (cyclotomic) symmetry; the **loose configs are
aperiodic** (high complexity). The cut-and-project / Meyer-set picture.
> **USE — a complexity/entropy bound:** the worry-set is the *minimal-complexity* locus
> (lowest configurational entropy = most ordered = the crystal). An entropy/complexity
> functional with the AP as its minimizer would bound the worry-set from the
> order-theoretic side (complement to the measure side).

## Corrected meta-reading: several frames detect the prime-2 doubling face

| frame | the n=14 obstruction is… |
|---|---|
| group-algebraic | the nonunit/non-etale characteristic-2 order-14 model |
| RG | the relevant `2`-eigenvalue of the doubling flow |
| hyperoperation | the additive (rung-1) AP at `2q` |
| rigidity (S585) | the `⟨×2⟩` fragmentation of the witness orbit |
| parity (S587) | the even-prime, odd/even sector seam |
| triangular (S586) | the `2³` in `8·C(n,2)+1=(2n-1)²` |

> **Five frames, one narrower conclusion:** prime `2` is a robust target for the
> doubling/nonunit face. It is not the only known obstruction, and the odd part
> does not by itself supply the missing safe seed or pointwise exit. The frames
> are candidate sidecars to combine with those obligations.

## Honest status

- **Grounded:** the stratification of the n=14 clock by `gcd(j,14)`;
  `Q(zeta_14)=Q(zeta_7)`; nonunit multiplication by `2` drains the witnesses;
  `M` is a covering radius, not a Mahler product.
- **Creative seeds (uses stated, untested):** the Hasse/local–global proof; the RG-flow
  rigidity proof; the additive-only certification; the complexity/entropy bound.
- **Corrected meta-claim:** several frames localize the doubling face to prime
  `2`; no prime-2-only reduction of LRC14 is claimed.

**Artifacts:** `04-computation/lrc_adelic_tangent_s589.py` (+`.out`). Builds on S588
(exp/cyclotomic), S585 (doubling rigidity), S586 (triangular), S587 (parity), THM-398
(C′), S547 (adelic, repo). New: **HYP-2131**.
