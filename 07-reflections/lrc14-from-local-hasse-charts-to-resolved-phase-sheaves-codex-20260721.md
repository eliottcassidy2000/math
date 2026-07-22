# LRC14: from local Hasse charts to resolved-phase sheaves

**codex-2026-07-21.** This pass began with the proved NC2 whole-face
Frobenius mechanism (THM-2022/2041), pulled the new support-transform and
JC/LRC work, and then searched the older monodromy, recurrence, confluent,
tree, code, and tournament-composition lanes. The concrete result is
THM-2043. LRC(14) remains open.

## What actually landed

Characteristic-seven Hasse jets are exactly the right coordinates for a
function already living on `Z/14Z`:

```text
F_7[C_14]
  ~= F_7[X]/(X-1)^7 x F_7[X]/(X+1)^7.
```

They are not a way to undo the quotient from integer speeds to residues.
The sharp counterexample is

```text
A   = {1,...,13},
T_n = {1,...,11,13,96+3444n}.
```

Every `T_n` has the AP's full owner-indexed mod-14 data, all raw phase
functions and Hasse jets, the same blockedness mask for `q=2,...,13`, and the
same `q_threshold=14`. Yet `A` is tight and every `T_n` has

```text
C_(41,17)(T_n)=1,
min_v ||17v/41||=3/41>1/14.
```

Even a fixed finite number of `7`-adic lift-height digits can be made equal
to the AP while this margin-one exit survives. Thus the bad-prime idea has
been sharpened from a vague possible route into two exact statements:

1. Hasse jets are complete **local coordinates**.
2. No fixed local truncation is a **global LRC carrier**.

The positive replacement is an exact owner height or an adaptive resolved
phase `(q,a,C_(q,a))`. On bounded heights `v<=181`, owner-labelled residues
modulo `13` and `14` recover `v` by CRT. On unbounded speeds, a fixed pair of
moduli cannot do this.

## Why the unrelated support work mattered

THM-2000's central correction was that an integer sequence is its atomic
support measure, not one scalar mass or one indexed presentation. The LRC
application is literal. If

```text
mu_S = sum_(owner i) e_i delta_(v_i),
```

then every period-14 phase object is computed after the pushforward
`pi_14* mu_S`. The direction `delta_(v+14)-delta_v` is already in the kernel.
Fourier, Ramanujan, and Hasse transforms merely give different coordinates
on the pushed-forward object; none can restore that kernel direction.

This suggests that the faithful object is not one phase packet but a sheaf of
owner-labelled phase views over denominators. The local sections are

```text
P_q(S) = owner-labelled residues and endpoint currents at denominator q,
```

and a transition is useful only when it preserves a positive integer slack

```text
C_(q,a)(S)=min_v(14 dist(av,q)-q).
```

The `q=14` Hasse section describes the bad-prime boundary; the `q=41` section
exits the strongest alias family. This is the cleanest current interpretation
of the repo's recurrent `{14,27,41}` denominator atlas.

## Unrelated mechanisms tested against the seed/exit gates

| Source | Exact mechanism worth importing | What blocks a direct LRC proof |
|---|---|---|
| THM-1605, orbit-product TNC proof | Transport a local branch identity around a connected monodromy orbit, multiply all conjugates, and contradict a global Vieta invariant without choosing one branch. | A positive LRC margin is not presently a multiplicative fiber invariant, and the endpoint cover can disconnect on tie strata. |
| THM-2033, confluent Vandermonde | Replace a vanished distinct-node determinant by a derivative/Hasse row at the collision. This is a rigorous boundary-to-next-jet operation. | A nonzero signed determinant does not imply a positive safe phase unless endpoint owner and height survive. |
| THM-1775, tournament rational-moment case | Cayley--Hamilton turns infinitely many trace vanishings into finitely many checks. A nonsingular recurrence for a resolved LRC certificate would be a genuine finite-depth theorem. | LRC's minimum over runners and supremum over phase is not a linear moment sequence; no safe nonsingular recurrence is known. |
| THM-1960, regular-seed substitution | Modular decomposition and regular-seed spectral splitting reduce a large object to prime seeds while preserving all even spectral moments. | Irregular seeds mix block means with internal modes, and no residue-module quotient has yet been proved to preserve loneliness. |
| THM-1745-leaf-graded-arborescence-filtration-721-shadow | A prime-order automorphism acts freely below a leaf-depth threshold, forcing divisibility of low-defect witness-tree strata. | Endpoint-labelled LRC witnesses can have stabilizers; a faithful witness-tree construction is missing. |
| THM-856, positive minimum-tree initial form | Keep an entire tied minimum face; positivity makes its aggregate seed automatically nonzero. | LRC owner currents are signed, so the automatic positive seed disappears. This isolates the missing theorem cleanly. |
| THM-1735, resultant good-place selection | Once a characteristic-zero polynomial seed is known, a resultant leaves only finitely many bad reduction primes. | It selects a good place but does not create a safe/dual seed. |
| THM-1435, cotangent/dual sidecar | Attach a dual variable so an obstruction enters a category with stronger algebraic certificates. | A naive LRC dualization can preserve a signed shadow while losing pointwise positivity; the analogue must retain `C_(q,a)>0`. |
| THM-346, prime-step cube transport | Whole aggregate transport survives a prime power even for non-equitable quotient buckets. | It is a preservation theorem, not seed production, and characteristic two is nonsemisimple for the natural cube algebra. |

The strongest two imports are therefore support-measure discipline and
orbit-product transport. The former already produced THM-2043. The latter is
a serious possible global mechanism if a multiplicative resolved-margin
object can be found.

## A general bad-prime theorem for other repo problems

THM-2043 extends verbatim. If `N=M p^a` and `p` does not divide `M`, then

```text
F_p[C_N] = F_p[X]/((X^M-1)^(p^a)).
```

Over a splitting field, every `M`th root of unity carries `p^a` Hasse
coordinates. This is the natural packet for repeated-root cyclic codes and
bad-characteristic periodic signals. It also clarifies two other lanes:

- For the binary cube, with `y_i=g_i-1`,
  `F_2[(C_2)^m] ~= F_2[y_1,...,y_m]/(y_1^2,...,y_m^2)`. The analogue is a
  multivariate nilpotent/Hasse filtration, compatible with THM-346's
  aggregate walk congruence.
- For higher-dimensional Gaussian moments, multivariate Hasse data may retain
  a tied reduction layer, but THM-2022's three-factor counterexample still
  forbids replacing the required vector valuation by one scalar face.
- For Poisson/Jacobian work, characteristic-`p` Frobenius makes `p`th powers
  central and erases the bracket. Jets may retain an associated-graded layer,
  but they cannot by themselves supply a canonical pair or quantization.

The reusable lesson is the same in every domain: a complete coordinate system
on a quotient is not a section of the quotient map.

## The new LRC proof target

The shortest live target remains THM-671's resolved-modulus `B5` supply:

```text
every reduced residual S
  -> some resolved (q,a) or dual Fejer packet
  -> positive integer/measure slack
  -> strict lonely phase,
```

with AP/Goddyn--Wong equality atoms and K33/state-lift packets split off by
labels. THM-2043 contributes a precise local handoff:

```text
q=14 parity-Hasse boundary chart
  -> retain owner and exact height, or move to an adaptive q
  -> q=41 unit-excess margin when available.
```

The HYP-2979 stress bank makes this plausible but not proved globally: all
21,906 tested AP-neighborhood rows had a primitive weak phase by `q<=42`, and
only AP/GW lacked a strict one; the eleven named rows show that raw q=14 jets
mix boundary, direct-q, petal, and K33 routes in one fiber. The theorem needed
is a structural supply theorem for the resolved phase, not a larger census.

## The incoming JC/LRC comparison, kept honest

The new S205 reflection correctly emphasizes the AP as a shared cold extremal
and correctly warns that Frobenius is a low-rank preservation mechanism. Its
useful LRC consequence is diagnostic: rigidity should be proved on the AP
fiber before attempting generic harmonic estimates. It is not, by itself, a
functional transfer from Jacobian or Gaussian nullcones to LRC. The
oscillating inequality and integer height data are exactly what those
algebraic nullcones forget.

## The incoming S206 disproof search, corrected and reused

S206 contributes a useful Wall-A search coordinate: THM-730 makes the AP the
additive-triple extremizer, so a higher-order good-set autocorrelation profile
is a plausible place to look for either AP rigidity or a competing shape. It
does not follow that every hypothetical counterexample must win through that
mechanism, and AP extraction is a sufficient route rather than an iff
reformulation of LRC(14).

MISTAKE-221 also corrects the finite search instrument. A bounded rational
scan gives `L_Q(S)<=M(S)`; it is exact under the simple complete-breakpoint
condition `Q>=2 max(S)`, because every maximizer lies on a pair-sum ruler
`t=p/(v_i+v_j)` by THM-1002-pair-sum-denominator-bound §1. Exact pair-sum replay shows the old large-blocker value was
strictly low:
`L_1200({1,...,12,5460})=92/1197 < M=420/5461`. The corrected program now
computes all fifteen rows exactly and finds them safe. This finite bank still
cannot establish a global minimizer.

THM-2043 makes the carrier choice concrete: the full raw `q=14` quotient does
not distinguish the AP from the strict `q=41` family, so an unlabelled
autocorrelation summary cannot be assumed faithful either. The resolved slack
`C_(q,a)` is the theorem-facing output; autocorrelation is potentially a
supplier of such a phase, not a substitute for the exit certificate.

## Ranked next experiments

1. Rewrite the HYP-3036 scheduler in terms of exact slacks `C_(q,a)`, not only
   primitive packet counts, and retain endpoint owners across `q=14,27,41`.
2. Search for a multiplicative orbit product of resolved margins or endpoint
   factors to which THM-1605's transitive-monodromy argument applies.
3. On the zero-open `q=14` fiber, test a confluent owner-current determinant
   inspired by THM-2033, with the `q=41` margin as the required exit oracle.
4. Attempt a finite recurrence only for the labelled denominator scheduler,
   not for the raw min--max function; audit singular leading coefficients.
5. Formalize the elementary THM-2043 CRT/Hasse/no-go core. It is small and
   would prevent future work from silently upgrading local completeness to
   global faithfulness.

## Assumption challenge and Tournament Analysis

Candidate vertices were runners, residue classes, denominators, phases,
Hasse orders, owner-height germs, monodromy branches, tree strata, and proof
obligations. The useful vertices are proof carriers. The comparison observable
records strict-witness status, endpoint ownership, magnitude, full local
phase, and jet depth. Its declared tournament is transitive; that is a
retention hierarchy, not evidence of LRC transitivity.

The challenged assumption was that more local depth eventually becomes
global information. THM-2043 refutes it even after adjoining any fixed number
of `7`-adic height digits. The exit must move sideways to another denominator
or retain the exact lift.
