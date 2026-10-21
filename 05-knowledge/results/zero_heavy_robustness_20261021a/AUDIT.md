# Same-agent audit and exact controls

**Status:** same-agent proof/implementation review complete; independent
root review is pending. No claim of independent audit is made here.

## Proof checks

- Total-weight deficit, heavy-coordinate deficiency, and feedback-vertex
  number are kept distinct throughout. The arc weights are positive
  integers; r=k is forced only while total deficit is <2k+1.
- The corrected collision charge is k+1 per missing endpoint. At deficit k
  there is no P3 term; at k+1 there is exactly one P3 plus matching edges.
- All isomorphisms preserve weights. Large vertices are fixed in their
  intrinsic transitive order; the only nontrivial unit-cell symmetry in
  the proved ranges is the explicitly retained internal C3 rotation.
- An empty matched gap is rejected because it creates homogeneous pairs.
  An early checker overasserted that there was exactly one such pair.
  The rejection exposed the larger transitive clone chain involving the
  neighboring gaps. That unnecessary assertion was removed before the
  successful receipt; the theorem requires and proves non-pair-freeness,
  not a unique pair. No bad record is retained or silently contracted.
- An exterior empty gap creates an exterior unit singleton, not an absent
  unit cell. The proof's recovery sentence uses the former, correct form.
- Strongness is crossing every proper cut of the endpoint order, not
  ordinary connectedness of the backward graph. The P3 recurrence is
  checked by exact graph generation independently of the class generator.
- The two displayed boundary formulas' one-unit-gap proofs have explicit
  k>=3 and k>=4 scopes. Their unified formula has k>=3 and
  0<=delta<=2k-3 scope, with absent negative-index binomial terms.
  Additional k=2 and k=3 boundary cases are marked FINITE-EXACT rather than
  silently importing that proof into multiple-unit-gap cells.
- The initial hand calculation of P_5 was corrected to 3012 before receipt
  publication; its recurrence and independent generation agree.

## Computation universe and controls

zero_heavy_check.py generates all three presentation types at deficits
k and k+1 for k=2,3,4, with and without the proper-cut filter. It tests
actual homogeneous pairs, sums actual weights, tests actual strongness,
and canonicalizes by every permutation within complete unit-incidence
cells. Distinct accepted records have distinct keys. Reversing labels
checks covariance; subset DP independently checks weighted feedback on
small retained examples. Separately, all perfect matchings and all
P3-plus-matching graphs through k=5 verify the irreducible recurrences and
adjacent-pair incidence. The k=2 consequences agree with the pinned
independent complete compiler receipt. No fit or interpolation is used.

palette_check.py verifies dependency hashes, reconstructs every edge of
the source gate from exact rational biquadratic coordinates, enumerates
all six root fibres before and after one explicit port change, and emits
literal colourings. It builds an actual six-copy geometry at rotations
2,3,4,5,6, testing all 4,095 unordered distances and exactly 186 unit edges.
Four-colour and nearest-repair six-colour witnesses are verified on those
actual edges. There are 20,210,1540,8855 palette multisets at N=1,2,3,4;
all repaired-copy subsets are checked separately from the load formula.
This is not an exhaustive enumeration of every arbitrary boundary word,
which is not needed for the proved lower bound.

Both normal and python -O replays completed with byte-identical JSON
receipts. All gates use explicit exceptions; no removable assert appears.
The exact receipt hashes and source hashes are pinned in manifest.json.
The first palette script attempts exposed a Python generator syntax error
and a missing math.comb import; both were fixed before either successful
receipt. No partial output is offered as a completed result.

## Remaining frontiers

- At delta>=2(k-1), multiple unit gaps can interact. The small instances
  here do not prove a general exact owner grammar in that range.
- At delta>=2k+1, r<k enters, with endpoint weights no longer forced unit.
- The unrestricted zero-heavy interface and A000568 remain open.
- Root-only repetition produces unbounded conditional robustness in a
  four-colourable graph. An unconditional six-colour-forcing context or
  a source satisfying every five-colouring circle-saturation predicate
  is not supplied.
