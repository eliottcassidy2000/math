# Unknotting number one: corrected status and certificate engine

**Status (audited 2026-07-28):** general decidability of `u(K)=1` remains
OPEN.  The in-repo program is a sound three-valued certificate engine only
after the exact-one/at-most-one repair below; it is not the Boolean,
under-one-hour decider requested by the Epoch open problem.

Primary status source:
https://epoch.ai/frontiermath/open-problems/unknotting-number

## The distinction the algorithm must preserve

For a diagram `D` representing `K`, let `u(D)` be the fewest crossing
changes in that particular diagram and

```text
u(K) = min {u(D') : D' is any diagram representing K}.
```

Changing a visible crossing of `D` and certifying the result is an unknot
proves only `u(K) <= 1`.  Exact equality also requires a certificate that
the input knot is nontrivial.  Conversely, finding no visible unknotting
crossing proves only `u(D) > 1`, not `u(K) > 1`.

Taniyama proves that for every nontrivial knot `K` and every natural `n`
there is a diagram `D` of `K` with `u(D) >= n`; hence visible-crossing
enumeration cannot be complete even when `u(K)=1`.
https://arxiv.org/abs/0805.3174

McCoy supplies an important exact exception: if `K` is alternating and
`u(K)=1`, every alternating diagram has an unknotting crossing.
https://arxiv.org/abs/1312.1278

Coward--Lackenby classify genus-one `u=1` knots as doubled knots and
classify their unknotting crossing changes.  The cited paper does not by
itself establish the previously claimed general-purpose practical
arbitrary-PD decider for the genus-one input promise.
https://arxiv.org/abs/0809.4142

## Exact soundness correction

The original `TRUE_CERTIFIED` rule was:

```text
input greedy R1/R2 stalls
+ one visible crossing change greedily reduces to the empty diagram
=> u(K)=1.
```

The implication is false because a stalled greedy simplifier does not
prove that the input is knotted.  Exact six-crossing hostile witness:

```text
[[1,11,2,10],[6,10,7,9],[3,8,4,9],
 [11,5,12,4],[7,2,8,3],[5,1,6,12]]
```

This diagram is obtained from the one-crossing unknot by the legal reverse
Reidemeister path `R1+, R2+, R2+, R3, R3`, so it represents the unknot and
has `u(K)=0`.  The original engine stalls on the input, then changes
crossing 4 and reduces by `R2, R1, R2, R1`, returning the false verdict
`TRUE_CERTIFIED`.

Reproduction and full intermediate PD path:

```bash
python3 .scratch/unknotting_mu3_20260728/unknot1_true_soundness_witness.py
```

**Minimal safe repair:** emit `TRUE_CERTIFIED` only when both:

1. the changed diagram has a valid unknot certificate; and
2. the input has an independent nontriviality certificate.

The repaired stdlib core currently uses `det(K) != 1` or `sigma(K) != 0`
for (2).  If a visible change is found but both invariants equal their
unknot values, retain the change as an `u(K)<=1` upper-bound certificate
and return `UNKNOWN`.  An exact unknot oracle on the input would remove
this conservative loss.

Proposed patch and tests:

```text
.scratch/unknotting_mu3_20260728/proposed/unknot1_decider-soundness.patch
.scratch/unknotting_mu3_20260728/proposed/test-normal.out
.scratch/unknotting_mu3_20260728/proposed/test-O.out
```

All 13 original checks and 4 hostile checks pass (17/17); normal Python
and `python3 -O` outputs are byte-identical.  Executable `assert`
statements were replaced by optimization-stable `require` checks.

## What the current engine actually computes

1. PD validation, faces, and checkerboard coloring.
2. Exact Goeritz determinant and Gordon--Litherland signature.
3. Sound theorem-level lower obstructions, conditional on the exact
   implementation:
   - `|sigma(K)| >= 4` implies `u(K) >= 2`;
   - noncyclic `H_1(Sigma(K))` contradicts the Montesinos surgery form;
   - a cyclic linking form that does not represent `+/-2/d` contradicts
     Lickorish's necessary condition.
4. Visible crossing changes, screened by determinant/signature, followed
   by a positive greedy R1/R2 unknot certificate.

The code does **not** compute the Alexander polynomial or a Smith normal
form despite the earlier pipeline summary.  It uses an adjugate/determinant
divisor test for cyclicity.  Its bounded generator candidate search can be
incomplete, but fails safe: it returns no obstruction if it cannot find a
generator.  An independent audit checked 400 PD relabel/order variants and
1,000 symmetric integer presentations against SymPy Smith normal form.

A three-knot regression gate is useful but does not by itself confer
universal “decision rights”; the decision is justified by the cited
theorem plus a correct exact implementation.  The finite audit is:

```bash
python3 .scratch/unknotting_mu3_20260728/false_certificate_audit.py
```

## The owner's example is K11n3

An independent traversal gives signed DT code

```text
[4,8,10,-14,2,-16,-20,-6,-22,-12,-18],
```

which encodes as `bdeGaHJCKFI` and is entry `K11n3` in the official
KnotTheory table.  Knot Atlas gives the identical PD and DT code, with
`det=43`, `sigma=-2`, and current unknotting-number range `{1,2}`:
https://katlas.org/wiki/K11n3

Thus `UNKNOWN` is exactly the honest current answer.  All eleven visible
flips have determinant different from one, proving `u(D)>1` for this
diagram, but not deciding whether `u(K)` is one or two.

Reproduction:

```bash
python3 .scratch/unknotting_mu3_20260728/k11n3_dt_audit.py
python3 04-computation/unknot1_decider.py --example
```

## Strongest honest algorithmic statement

There is a complete positive semidecision for exact `u(K)=1`:

1. run an exact unknot recognizer on the input and return false if `u=0`;
2. breadth-first enumerate every finite Reidemeister sequence from the
   input diagram;
3. at every reached diagram, change each crossing and run the exact unknot
   recognizer;
4. return true with the move/change certificate when one succeeds.

Soundness is immediate.  If `u(K)=1`, the defining diagram is connected to
the input by a finite Reidemeister sequence, so the search eventually
halts.  For `u(K)>=2` it need not halt.

Lackenby's new hierarchy algorithm improves the exact unknot-recognition
oracle in steps (1) and (3), but supplies no finite bound on all relevant
crossing arcs or Reidemeister exploration:
https://arxiv.org/abs/2607.23350

No cited result yields a worst-case one-hour guarantee at 100 crossings.
The reinforcement-learning work often finds upper bounds up to 200
crossings but is not an exact negative decider:
https://arxiv.org/abs/2409.09032
