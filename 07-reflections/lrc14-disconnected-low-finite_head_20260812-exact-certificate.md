# Disconnected-low LRC14: exact finite-head certificate

**Status.** FINITE-EXACT.  This note discharges the raw `p<264` head in the
small-ruler, moderate-ratio, primitive non-`3:5`, common-dilation `g<=3`
lane.  It does **not** certify any of the `22,890` affine rays with raw
`p>=264`, and therefore does not by itself close disconnected-low or LRC(14).

## Inherited boundary

The closest proved mechanism is the exact midpoint lower bound in
`lrc14_connected_low_uniform_high_forest_tail_thm3350.py`; the physical
overlap semantics are those of
`lrc_general_reflected_pair_mass_thm3352.py`.  MISTAKE-376 is the relevant
correction: neither the finite head nor the affine-ray bank was allowed to be
routed as closed before an actual certificate existed.

Two hostile controls are retained.  The excluded `3:5` lane has its separate
all-dilation proof and near-minimum `158/46397`; a low, reversed `(10,2)`
channel at `(L,j,e,f)=(168,90,1,3)` has physical mass zero and checks that the
compiler is not silently imposing the desired lower bound on every pair.

## Frozen universe

The compiler extracts, from the hash-pinned THM-3350 atlas, every body-safe
ordered context `(L,j,e,f)` with `L<4592`.  The resulting universe has

* `29` rulers, `2,530` contexts, and `1,304` groups after forgetting only the
  cell coordinate `j` used by the exact mass evaluation;
* primitive channels `P<Q<8P`, `P+Q>=8`, `(P,Q)!=(3,5)`;
* common dilation `g in {1,2,3}` and raw level `p=gP<264`.

There are `148,110`, `36,978`, and `16,286` channels at the three respective
dilations, hence `201,374` channels and exactly

`509,476,220`

physical context-channel rows.  The context and channel universes have
independent SHA-256 digests in the frozen output.

## Exact branch-and-bound compiler

For a grouped context `(L,e,f)` and channel `(g,P,Q)`, put

`C=|Qe-Pf|`

and use the already-proved midpoint estimate

`I >= Phi(P,Q)-E`,

where

`Phi(P,Q)=max(1/105, 1/49-12/(49PQ))`

and

`E = gamma(P)eP/(gLP-e) + gamma(Q)fQ/(gLQ-f)
     + C(floor(C/L)+1)/(2g^2LPQ)`,

with `gamma(k)=1/2` for even `k` and
`gamma(k)=(k^2+1)/(2k^2)` for odd `k`.

No floating-point decision is made.  With `S=2^56`, the C++ compiler forms

`floor(S Phi) - ceil(S E_1) - ceil(S E_2) - ceil(S E_3)`.

Dividing this integer by `S` gives a rigorous lower bound for `I`.  It clears
`498,526,403` rows directly at the target `1/294`.  Every one of the remaining
`10,949,817` rows is evaluated by a separate `__int128` Euclidean floor-moment
port of the physical mass engine.

Certifying only the target would not certify the reported global minimum.  An
actual row in the declared universe has value `92/7645`; it supplies an upper
cap for the minimum.  The compiler retains every target-safe grouped row whose
dyadic lower bound lies below that cap, then evaluates the corresponding
`36,276,868` physical rows.  Any omitted row has a proved lower bound at least
the current cap; if the candidate falls, that omission remains valid.  Thus
the minimum claim is exact branch-and-bound, not a minimum over only the
target failures.

The ordered `201,374`-row channel summary hashes each literal mass evaluation
and is frozen as

`85ff8655e647585c9113bdf204807d62e1c9e0243e50ec93971b592f9b46949e`.

## Result

There are no failures below `1/294`.  The exact global minimum is

`I = 92/7645`

at

`(g,P,Q;L,j,e,f)=(1,4,5;168,90,12,6)`.

Its strict margin is

`92/7645 - 1/294 = 19403/2247630 > 0`.

The 99 deterministic controls agree along three routes: the C++ compiler,
the canonical Python THM-3352 engine, and a definition-level two-pointer
intersection of the literal clipped tooth intervals.  They include both the
`3:5` near-minimum and the zero-mass low-channel control.

After this replay, concurrent commit `9936848e6` arrived on
`origin/codex/lrc-math-20260812`.  Its separate integer-only C scanner performs
a literal bulk scan of the strict superset obtained by restoring the three
excluded `3:5` dilations.  It has the same context digest, reports zero
failures, and finds the superset minimum `158/46397` at `g=1`, `3:5`; this is
`1/294+55/13640718`.  That computation is an independent corroboration, not a
dependency of the certificate here.  Its inclusion of `3:5` explains its
`201,377` channel count versus the present non-`3:5` count `201,374`.

## Proof scope and remaining horn

This certificate composes with the existing large-ruler, `q>=8p`, all-`3:5`,
and non-`3:5` `g>=4` results, but it deliberately does not duplicate their
proofs.  In the arbitrary-ratio Dirichlet partition it removes the complete
raw `p<264` head.  The disconnected-low frontier is now concentrated on the
separately defined `22,890` nonzero-resonance affine rays with raw `p>=264`.
Those rays remain OPEN until their own exact/analytic certificate is supplied.

Reproduce from the repository root:

```text
python 04-computation/lrc14_disconnected_low_finite_head_20260812.py --threads 12
python -O 04-computation/lrc14_disconnected_low_finite_head_20260812.py --threads 12
```

The two modes reproduce the frozen semantic digests and exact output.
