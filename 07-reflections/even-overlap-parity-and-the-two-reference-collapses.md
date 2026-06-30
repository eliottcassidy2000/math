# Even-overlap parity and the two reference-collapses: how H's variance mirrors the LRC floor

*klein-2026-06-29-S5. Extending the owner's CV(H)^2 identity (mac-mini proved it -> 0; THM-589 gives the exact even-run form and the 2/n rate). Two relations between the Hamiltonian-path second moment and the LRC floor fall out, both worth stating.*

## Relation 1: even-overlap survival IS the 2-adic descent, on the witness side

THM-589: two Hamiltonian paths sharing a run of `m` consecutive-integer arcs contribute to `Var(H)`
with a factor `1 + (-1)^m` — `2` if `m` is even, `0` if `m` is odd. **Odd overlap cancels; even
survives.** That is exactly the signature of the LRC covering floor's 2-adic descent (THM-580), where
the doubling map peels the even part cleanly and the odd part is what couples. The project has the
"even survives, odd is the obstruction/peels" pattern on three faces now:

- **floor side** (THM-580): the 2-adic descent peels even speeds; the residue is the odd (apex-7) core.
- **metagraph side** (THM-584/587): complement = the antipodal map; `R`-even (even hypercube levels) is
  the Brouwer/SOS bulk, `R`-odd the Borsuk-Ulam obstruction.
- **witness side** (THM-589): `H`'s variance survives only on EVEN overlap runs; odd runs cancel by
  orientation parity.

All three are one 2-adic / antipodal `Z_2`. The witness-side version is the cleanest to compute (it is
just a permutation succession statistic) and it concentrates: `CV(H)^2 ~ 2/n`.

## Relation 2: two reference-collapses — S_n (clean) vs Z_14 (set-dependent)

The owner's identity collapses the pair-sum `E[H^2] = sum_{pi,sigma} P(both paths)` to a SINGLE sum
against a reference path, using `S_n` relabeling symmetry. The LRC floor gatekeeper (THM-579) is the
SAME shape — a pair-sum `Var(N_R) = sum_{a,b} autocorr(lonely(R), (a-b)/14)` collapsed against the 14
sheets, using `Z_14` shift symmetry. Both are **reference-collapsed pair-overlap second moments.** The
difference is everything:

| | collapse group | result | behavior |
|---|---|---|---|
| `CV(H)^2` (witness) | `S_n` (transitive, no vanishing fiber) | one reference-path succession sum | clean, SET-INDEPENDENT, `~2/n -> 0` |
| `CV(N_R)^2` (floor) | `Z_14` (per-set sheets) | per-set autocorrelation sum | SET-DEPENDENT, unbounded as `m_R->0` (HYP-3554) |

So the Hamiltonian-path side is the rehearsal where the collapse works perfectly: `S_n` is transitive,
so the second moment becomes a single set-independent object and concentrates. The covering side fails
to be set-independent precisely because `Z_14` is not transitive over the speed structure and `m_R` can
vanish. **What the LRC floor needs is a collapse that restores `S_n`-like set-independence** — which is
exactly mac-mini's `Gamma_0(N)` congruence-subgroup route (HYP-3553): replace the per-set variance by a
subgroup-density quantity that depends only on `N=14`. The owner's CV(H)^2 formula is the proof-of-
concept that such a clean collapse exists when the symmetry is transitive; the floor's job is to
manufacture the transitive symmetry it lacks.

## What is genuinely new here (vs mac-mini HYP-3560)

mac-mini proved `CV(H)^2 -> 0` via Poisson(1) adjacencies. The Poisson limit hides the mechanism;
THM-589's even-run form exposes it: the variance is a sum over EVEN-overlap configurations, the odd ones
cancel, and the leading even config (a single length-2 run) gives the exact `2/n`. The parity
cancellation is what ties `H`'s fluctuation to the 2-adic descent and to the antipodal metagraph — it is
not visible as "Poisson(1)". And `A_n(2) = 1,2,8,32,158,928,6350,49752,439670,...` is a new integer
sequence (not in OEIS) — the suitably-normalized labeled-tournament second moment of the Hamiltonian-path
count.

See THM-589, [[the-finite-rehearsal-h-concentrates-and-poisson-gives-existence]] (HYP-3560),
[[the-variance-blows-up-where-the-fiber-vanishes]] (HYP-3554), [[the-2adic-descent...]] (THM-580),
[[complement-is-the-antipodal-map-of-the-arc-hypercube]] (THM-584), [[the-covering-is-a-congruence-subgroup]]
(HYP-3553).
