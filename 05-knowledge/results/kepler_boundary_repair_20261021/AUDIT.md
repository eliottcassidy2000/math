# Scoped independent audit: THM-4737

Date: 2026-10-21. Target: `01-canon/theorems/THM-4737-kepler-polygon-reformulation-of-the-lonely-runner-conjecture-observer-angular-momentum-vertex-times-and-the-relation-torus.md` as inspected at 01:42 UTC. No repository/source edits or Git mutations were made by this audit.

## K2: demonstrated closed-endpoint loss; exact repair

**REFUTED as an everywhere iff:** with n=3, D={1,2}, t=1/3, the phases are 1/3 and 2/3. Both have circular distance 1/3, so the observer is lonely. The document assigns phase 2/3 to side 2 using its half-open side convention. That side is incident to the observer, hence the assigned right-hand areal rate is zero. Indeed, for D={1,2}, the true lonely set modulo 1 consists of {1/3,2/3}, while the set where both assigned right-hand rates are positive is empty. Equal measures do not preserve nonemptiness.

**PROVED repair for n>=3 and positive speeds:** let A_i^-(t) and A_i^+(t) be the left and right derivatives of the observer-centred area clock. Then

    lonely at t iff max(A_i^-(t), A_i^+(t)) > 0 for every runner i.

At each far incident vertex exactly one derivative is positive. At the observer's vertex both are zero. On an incident-side interior both are zero, and on all other side interiors both are positive. At other vertices both are positive. This exhausts all phases and proves the equivalence for every real time, without a measure-zero exception.

The order of quantifiers matters: each runner may use a different side of t. The criterion is not that all runners share a common positive left or right derivative. In the minimal tight example the two runners require opposite sides. Replace the title and concluding informal conjecture rewrite too; qualifying only the body's old iff almost everywhere does not repair their exact existence claims.

The original right-rate set remains a valid almost-everywhere representative for computing measure. For integer speeds, statements that the tight lonely set is finite must be restricted to one period, not the entire real time axis.

## K3: rational maxima are not automatically lonely maxima

**PROVED repair:** define N_min(D) in {2,3,...} union {infinity}, with infinity when there is no lonely rational time. For positive integer speeds,

    N_min(D) < infinity iff kappa(D) >= 1/(k+1).

Rational attainment of kappa does not itself prove the right-hand side, so the document's unconditional finiteness argument assumes the conjecture it is supposed to reformulate.

**Stronger conditional survivor:** if V=max(D), kappa is attained at a time whose reduced denominator is at most 2V. Indeed, f(t)=min_i||v_i t|| is a continuous 1-periodic piecewise-linear function, with all linear-piece slopes among the nonzero values +/-v_i. Its maximum is positive. At a maximizing kink, an active increasing branch v_i t-a meets an active decreasing branch b-v_j t, for integers a,b. Thus t=(a+b)/(v_i+v_j); i=j includes a single tent peak. Consequently, conditional on kappa>=1/(k+1), N_min(D)<=2V. This is not an unconditional proof that N_min is finite, nor a speed-independent bound.

Alternatively, the closed lonely-time set in [0,1] is a finite union of intervals with rational endpoints, obtained by intersecting the unions [(m+1/n)/v_i,(m+1-1/n)/v_i]. A real witness therefore implies a rational witness for integer speeds, even when every component is a singleton.

## K4: keep the valid closure statement; separate actual attainment and density

**PROVED survivor:** the displayed equivalence between kappa_sup(D)>=1/n, intersection of the orbit closure K_D with the closed box B, and an integer point in M B is valid. Continuity on the compact closure identifies the supremum along the dense orbit with the maximum on that closure.

**GAP, not a claimed counterexample:** density alone does not establish that the actual orbit hits a closed target met only on its boundary. Therefore the generic displayed closure equivalence should be labelled a supremum formulation unless an attainment argument is supplied. This audit supplies no real-speed counterexample to the actual-time LRC and does not claim such a counterexample exists.

**Rank<=1 actual-time conclusion survives:** rank zero has an interior hit. For rank one and n>=4, the interval image of the open box has length ||m||_1(1-2/n)>=3(1-2/n)>1; it contains an integer. The orbit closure therefore meets the box interior, and density gives an actual hit. For n=3, k=2 and rank one, the speeds are commensurable and the orbit is itself closed, so a closed-box hit is attained.

For noncommensurable real speeds, the rank-zero formula (1-2/n)^k denotes Haar measure/asymptotic time density, not the unit-window Lebesgue measure defined for integer speeds in K2. Example: D={1/100,sqrt(2)/100}, n=3 has no lonely time in [0,1], because both phases stay below 1/3, but its asymptotic lonely density is 1/9. The formula is valid with the correct meaning.

## Reproduction and finite scope

Run `python3 session19-kepler-boundary-audit/check_kepler_boundaries.py` from `/workspace/scratch/07d5532f2fce`.

- 1,050 exact phase checks: n=3..32, every vertex and an interior representative of every side. The repaired equivalence always passes; the original right-rate criterion has exactly one endpoint mismatch per n.
- The minimal tight witness and both singleton lonely components are verified exactly with fractions.
- 154 rows: all distinct integer k-sets in {1,...,8}, k=2..4. Independent full piecewise-linear-envelope maximization agrees with pair-sum candidate maximization, and the conditional rational-witness statements pass.
- Hostile target control: D={1,2}, target 2/5 has an empty closed-time set, guarding against unconditional witness existence.

These bounded checks corroborate the proofs and expose the precise failure mechanism. They are not a new LRC proof or a broad audit of the theorem's other claims.

## Optimization-safe checker and frozen packet

All 13 original assertion statements were converted by Python AST transformation to explicit conditional raises. Each condition is evaluated once; a failure message is evaluated only on failure. No Assert AST nodes remain, so optimization cannot remove these checks. Normal and `-O` runs both pass and produce byte-identical standard output and JSON. A temporary forced-failure copy also raises AssertionError under `-O`, confirming the failure path remains active.

Portable reproduction, from this packet directory:

```sh
python3 check_kepler_boundaries.py
python3 -O check_kepler_boundaries.py
```

Portable byte-comparison command, also from this directory:

```sh
tmpdir=$(mktemp -d)
python3 check_kepler_boundaries.py > "$tmpdir/normal.json"
python3 -O check_kepler_boundaries.py > "$tmpdir/optimized.json"
cmp "$tmpdir/normal.json" "$tmpdir/optimized.json"
rm -r "$tmpdir"
```

`manifest.json` pins the four deliverable files and preserves the original checker hash and inspected-source lineage. It excludes its own hash; the publisher should pin the manifest independently. The mathematics and finite verification universe did not change during this hardening.
