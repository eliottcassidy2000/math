# The coding lens confirms the LRC support structure but not its analytic content

**Source:** kind-pasteur-2026-06-21. Dispatch: the owner's MDS-code / projective-arc lead
("an m-arc <-> a [m,N+1,m-N] MDS code"; "size 3 has 56 challenger shapes = 56 tournaments on
6 vertices"; arXiv:2501.19125; Spielman tc.pdf), applied to OPEN-Q-108 (the LRC(14) wide-cover
crux). Canon touched: HYP-2724 (this session), HYP-2719 (support-size seam), THM-538 (the
support-6 floor, corrected).

## What the lens gives

The carrier error `corr(E) = measS7(E) - iid_k = Sum_{n in Lambda(E)} K(n)` is a signed sum over
the **relation lattice** `Lambda(E) = {n in Z^k : sum n_i e_i = 0}`. Viewing `Lambda(E)` as a
linear **code** `[k, k-1, d]` with `d` = minimum support of a relation makes three things crisp,
all verified this session:
- **Extremal dichotomy:** the AP/consecutive set is the **anti-MDS** member (minimum distance
  collapses to 2; densest in low-weight relations) and is the verified hardest config; the
  Sidon/arc set is the **MDS** member (general position, no low-support relations) and is the
  easiest (`corr ~ 0`). The known "consec is the argmax" gets a coding-theoretic *mechanism*.
- **The driver is support-3** (3-APs `(1,-2,1)` + Schur triples `(1,1,-1)` = additive energy):
  `Pearson(A3, corr) = +0.93`. This matches mac-mini's HYP-2719a exactly.
- **A sign seam by support:** support-2 relations carry small/mixed sign, support-3 a large
  positive mass, support-4+ alternate — the combinatorial face of the "signed Erdos-Turan" the
  proof needs.

## The trap, and the lesson it re-teaches

A fan-out workflow confidently "REFUTED" the support-3 framing: it found `K(n)=0` for every
support `<= 5` (a clean algebraic vanishing `(1-1)^{6-|U|}=0`) and concluded the carrier is
support-6. **This is the already-conceded THM-538 error** (CASE-thm538-support6-floor-zero-padding,
MISTAKE-078). The floor holds for the **bare** vector (the active-coordinate sum `Q`), but `corr`
sums the **zero-padded** length-`k` kernel, whose zero coordinates carry `chat(0,T)=(1-|T|/7)`
factors that break the cancellation. Adjudicated on the workflow's own kernel:
`Kk([1,1,-1]) = 0` but `Kk([1,1,-1,0,0,0,0]) = +0.00066` — exactly the canon value. Support-3
relations *do* contribute; the AP's correction *is* support-3-dominated.

The lesson (twice now): **a clean vanishing can be the artifact of dropping a factor, and an
"exhaustive, max 5e-17" check can silently measure a different object than the one in the
statement.** The defense is not more compute — it is checking the verified quantity against a
case where the two candidate objects *must* differ (here, bare vs zero-padded), and against the
canon's record of past concessions. The workflow had the compute; it lacked the THM-538 memory.

## The transcendent point

Reframing reorganizes; it does not manufacture analytic content. Every lens this project has put
on OPEN-Q-108 — sector measure, singular series, Freiman dimension, the cut/cycle support seam,
and now the MDS/arc code — lands on the **same wall**: `corr(E)` is a *conditionally* convergent
sum, and its bulk lives in the **high-coefficient-height tail** (the slow reciprocal sum
`Sum 1/prod|n_j|`), not in any finite low-height enumeration. Concretely, with the correct kernel
the entire support`<=6`, `|coef|<=2` mass is only **~20%** of `corr` at the AP. The coding lens
sharpens *which* relations matter (support-3, anti-MDS) and *why* the AP is extremal, but the
remaining 80% — the conditional height tail — is the invariant nut: the R6-density / far-element
plateau (HYP-2645/2644), with the single-far case already closed (THM-546/547) and the
balanced-wide residual still open.

So the arc/MDS picture is the right *language* for the support structure (and it cleanly killed
the "56 challenger shapes = 56 tournaments" hope: the support-3 shape count is unbounded,
`56 = C(8,3)` is a coincidence with `A000568(6)`). It is not a new *estimate*. The proof still
needs the one analytic input every reframing has pointed to: a signed bound on the conditionally
convergent height tail of the relation-lattice sum. That the structure is now describable in five
independent languages, all agreeing, is itself evidence the wall is real and singular — the place
where the seven sectors, the additive energy, and the relation code meet.
