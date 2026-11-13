# Session293: neutral cycle-opening returns and the price of reactivation

**Proved in notes and independently reviewed. Chan–Pak remains open.**

The main advance is a distinction between returning the full count now and
returning the response to a future operation. A directed Hamilton cycle can be
opened indefinitely without changing the full Hamiltonian-path count, but the
renewal moves a named entrance port and accumulates a seam. A feedback arc can
see that seam even when every old-to-old endpoint count is unchanged.

| Result | Exact scope and proof |
|---|---|
| Neutral conveyor | [Classification and construction](opening/PROOF.md). The newborn sits immediately after the first strong component, the current entrance port moves to the newborn, and all old-to-old endpoint counts return. Preserving the conditioned opening count requires the endpoint-column support condition; a fixed entrance name prevents a second positive neutral opening. |
| Feedback budget | [Arbitrary middle cores and feedback sets](opening/MULTI_FEEDBACK.md). With a nonempty feedback set and a seam of order n, the full count is at least `1 + |E| 2^n`, sharp universally. Changing only the middle tournament cannot beat the transitive middle. Exact nonnegative coefficients count exterior slots, including slots needed to repair missing edges. |
| Compensated physical switches | [Exact two-switch carry law](switches/PROOF.md). Count-preserving switches in this carrier are converse-core transport. Unit repairs occur along one dyadic family apart from the short exceptions. The same exteriors around a general nonempty module pay the exact tariff `ΔH=2H(module)`, so simply inflating the module does not retain unit repair. |
| Exact arithmetic context ideal | [Uniform module probes](control/MODULE_IDEAL.md). The gcd of all uniform-context counts is the gcd of `(k−1)! k! pc_k`. Permuting exterior excursions supplies the extra factorial; alternating transitive exteriors give a unitriangular physical probe system. An asymmetric eight-vertex core has ordered-cover gcd 9 but context ideal 27. |
| Exact period and hidden state | [Sixfold transport](six_history/PROOF.md). Separate the rational quadratic three-cycle and its sixfold signed lift from an actual rational Möbius six-cycle, and from the clocks of 3x−1 Collatz. A positive conserved-mass realization carrying a primitive sixth mode needs at least six hidden states even when a three-coordinate observer suffices. Not a tournament realization. |

## Why this changes the constructive task

The existing five-core and fixed-defect constructions supply physical sources
and explicit endpoint data, not an arbitrary-target selector. Session24 already
shows that a fixed count admits only finitely many nontrivial strongly connected
skeletons. The conveyor uses the remaining freedom in singleton gaps. The
feedback bound now rules out even **cyclic** middle inflation as a way to
reactivate that conveyor indefinitely at a fixed target with the same uniform
exterior. This does not prohibit repeatable returns that replace an old seam,
change the exterior, or intentionally change targets between returns.

Session292 supplies an exact way to relax quadratic port guards. The next
proof obligation is a replacement that changes a nonuniform incidence or
replaces a core while preserving the prescribed full count and a usable next
interface. The coefficient compiler gives the repair equation. The module ideal
gives an all-context obstruction when an obstructing module stays fixed. A
primitive ideal permits signed arithmetic combinations, not positive coverage;
that distinction is where a new construction is needed.

## Independent review and reproduction

- [Opening, feedback, and switches](audit/REVIEW_OPENING_SWITCHES.md)
- [Uniform-module ideal](audit/REVIEW_IDEAL.md)
- [Rational and positive sixfold transport](audit/REVIEW_SIX.md)
- [Initial plan and inherited sources](PLAN.md)

Run the Python scripts listed in the reviews with both `python -B` and
`python -O -B`. The independent normal and optimized outputs agree. The
certificate universes test the formulas; their all-order scope comes from the
bijections, free actions, and integer inversions in the notes.
