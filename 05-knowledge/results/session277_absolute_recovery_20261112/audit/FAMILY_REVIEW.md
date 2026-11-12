# Independent family-lane review

Status: PASS WITH SCOPED CORRECTIONS (12 November 2026).
Scope: `families/PROOF.md` and `families/recurrence_bridge.py` only; this is
not a review of inherited repository conjectures or all historical results.

## Mathematical checks

1. The quotient-walk formula counts each path once. A run assignment is an
   ordered set partition of each transitive module; its blocks need not be
   intervals. The factor is therefore `k! S(r,k)`, not a binomial interval
   count. Zero-size modules and the empty tournament convention are retained.
2. The aggregate is explicitly multinomial-weighted over named modules. It
   must not be called a tournament-isomorphism census or a target spectrum.
3. The recurrence dichotomy is valid: a nonnegative rational walk series for
   a quotient containing a directed cycle has a genuine positive finite pole.
   Substitution `u=exp(t)-1` gives infinitely many genuine poles, because the
   derivative is nonzero at every preimage. This excludes D-finiteness and
   thus eventual P-recursiveness. No unsupported cancellation assumption is
   needed. The transitive case has aggregate `q^n`.
4. The regular-quotient geometric moment and translated-pole formulas have
   the stated factor `q/(d(d+1))`. The pole-tail integral is an absolute bound
   for every n>=1, so interval uniqueness, not heuristic rounding, certifies
   a recovered integer. The signed leading-pole error can oscillate.
5. For a nontrivial strong quotient, every vertex is in a directed triangle.
   Deterministic insertion is injective on Hamiltonian paths, so induced
   `C3[TT_r,1,1]` supplies the stated `2^r+1` lower bound. This proves the
   finite target box and fixed-bank sparsity, not universal Chan–Pak coverage.

## Corrections identified

- The target-inverse theorem originally said only "strong quotient".
  The singleton is strong under the usual convention and would contradict
  the exponential bound: `TT_r` has H=1 for every r. Restrict the theorem to
  q>=3, handle singleton separately, and explicitly reject c<=1 there.
- The initial pole routine used the exact answer to choose its numerical
  precision. The author is replacing this with a priori `B_n<=q^n n!`;
  exact Stirling evaluation remains an independent control, not certificate
  input. This distinction matters for describing the algorithm.
- The initial routine also had an arbitrary K<10000 assertion and a fixed
  precision allocation. The all-order mathematical convergence theorem does
  not by itself prove that this particular operational cap or fixed precision
  always succeeds. Either label the executable a bounded demonstration or
  add adaptive precision and remove the arbitrary K cap before claiming an
  executable all-order termination guarantee.

The first item is a mathematical boundary correction; the latter two are
implementation/claim-scope corrections. No failure of the principal compiler,
recurrence dichotomy, pole identity or nontrivial-quotient inverse was found.

## Independent controls

`python3 audit/family_independent.py` uses direct permutation enumeration,
not the author's Hamiltonian subset-DP routine, on72 seeded random quotient
and module-size instances (quotient orders1..6; total inflated order<=7).
All72 agree with the run compiler. Additional42 transitive aggregate checks
and120 regular-quotient checks using a separate binomial-convolution
recurrence agree. `FAMILY_CONTROLS.json` records the results.

Boundary pole controls q=3,5,21 and n=1..5 all returned exact certified
integers with one symmetric pole pair. These controls test boundaries, not
the all-order theorem by extrapolation.

## Priority and limits

The quotient-walk substitution input is explicitly inherited. The geometric
moment identity is a classical rational-exponential transform, not a claim
of discovering a new general transform. The useful new session junction is
its scoped recurrence-class classification, certified correction realization,
and explicit target-inversion consequence for this family. No claim that
these results settle Chan–Pak or enumerate A000568 is licensed.

## Final repair confirmation

The final author source has now been inspected: q>=3 and the singleton
exception are explicit; precision is chosen from the a priori bound
`B_n<=q^n n!`; the exact Stirling sum runs only after interval certification;
the arbitrary pole cap is removed; and ambiguous intervals with analytic
tail<1/8 trigger a restart with doubled extra precision. These repairs resolve
the three listed findings. For each fixed n,q, the analytic tail eventually
falls below1/8, and finite-precision enclosures for that finite prefix converge
as precision increases, proving termination in the arbitrary-precision
mathematical model (ordinary finite-machine resources still apply).

The subsequently added Mersenne face has the correct four run multiplicities
3,1,1,1 and formula `2+(2^a-1)(2^b-1)`. Twenty-five separate subset-DP controls
for a,b=1..5 pass. Complete independently rerun finite boxes exclude25,27,29,31
from positive transitive C3 inflations. The box theorem makes these all-size
family exclusions; they are expressly not tournament-spectrum exclusions.
