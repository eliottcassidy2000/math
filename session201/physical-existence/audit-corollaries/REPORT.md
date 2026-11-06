# Independent audit: asymmetric physical budgets and J crossings

**PASS.** The real/integer asymmetric inventory criterion and conditional J-crossing bound are valid in their stated scopes. The fixed-word J sharpness family is physical for k>=1, as the frozen note now explicitly requires. These are independent corollaries of the already reviewed physical adjoint and adjacent-cut framework, not existence results for absent targets.

## Frozen inputs

`INPUTS.json` records all input hashes, including the imported main proof and verifier. The companion files are:

- `COROLLARIES.md`: `e84a9f070844a830825e485de786a9eb1ae56270d6b197aefbcad0bb9539bf2d`.
- `verify_corollaries.py`: `55dd20f9bdbe5ee613f18c0ead47bad1f95b9ec8fd6eae1b2822362ff1786670`.
- Both author companion receipts: `154fb484b2aafafe03dc228e551880e8d2ea6f12e842ecc63fcf0c9bd1e18880`.

The main proof and checker have exactly the frozen hashes independently reviewed in the separate adjacent-cut audit. The author companion checker imports that main checker, so its bytes are included as a reproducibility dependency.

## Asymmetric budgets

The source-word cut sequence has initial left height two and strictly increasing left heights. Thus L>=2 ensures that the last cut with left height at most L exists. If it is terminal, its right boundary has seed height h0<=R. Otherwise the next left height exceeds L. If the preceding right height exceeded R, the product of these two adjacent heights would exceed LR>=3(n-1), contradicting the all-order adjacent-letter inequality and positive integral deletion product.

This proof uses the next left and previous right heights at one common letter; it does not apply a product bound to independent cuts. Consequently that last admissible-left cut simultaneously meets both budgets. The same source-word split preserves total J cost. In the joint case it transports b and z together, so both target readouts and the full auxiliary state are retained.

Real budgets and integer inventory cutoffs are distinguished correctly. Since source heights are integral, U<=L is equivalent to U<=floor(L), but the hypothesis is the product LR before rounding. The floors need not have product at least 3(n-1), and the proof never assumes they do. For example, the full empty-word target n=3 permits L=5/2 and R=12/5, although floor(L)floor(R)=4<6; the physical pair t,t still fits. The separate balanced singleton target three is outside the triangle-precursor clause.

The criterion is exact only as a representation-preserving search reduction for a target known to have a physical word, together with the ordinary converse that a found physical pair is a witness. It does not make a positive fiber exist.

## J-only coefficient and sharpness

Independent raw cone-basis multiplication gives the displayed matrix C=P^T J^T QP and row A=(2,1,4,3). The residual matrix 16C-7A^T A is exactly the displayed entrywise nonnegative matrix. Hence

    U(Jv)U(Jw) <= (16/7)K(Jv,w)

holds on the stated positive cone and therefore on all actual physical sources. No positive-cone vector is treated as physically realizable solely by this calculation.

At an actual height crossing whose common letter is J, the two adjacent maximum heights are U(Jv) and U(Jw), since J is self-adjoint. The smaller square is bounded by their product, proving the conditional cutoff floor(sqrt(16n/7)). The integer expression isqrt(floor(16n/7)) is identical. The proof does not assert that every source word or target has a J crossing.

For k>=1, the source t.Z^k is `(2^k+1,1,2^k+1,1)`, and its J-image is exactly `(2d-1,d-1,2d,d)` for d=2^k+1. Then `K(Jv,v)=7d^2-2d` and U(Jv)=4d-1. The ratio tends to 16/7, establishing sharpness on actual sources. At k=0, t.J=(3,1,3,1), so the displayed d-parametrization would be false; the frozen note explicitly excludes that case and warns about it.

For the specified word Z^k JJJ O^k, the two central cuts are (v,Jv) and (Jv,v). Both have maximum 4d-1. Earlier cuts have right heights at least that large, and later cuts have left heights at least that large, by strict monotonicity. Therefore its best possible cut maximum is exactly 4d-1. This is fixed-word sharpness, not a lower bound over all different representations of the same target.

## Independent controls

Both author runs were replayed and match both frozen normal/optimized receipts exactly.

The independent checker imports no author code. It reconstructs generators, physical cuts and cone matrices directly, uses exact rational budgets, and enumerates literal directed paths by DFS. Normal and optimized independent receipts are identical.

- Complete full/balanced word universe through length seven: 6,560 source-boundary cases.
- 136,692 rational rectangle tests, choosing thresholds at every left cover and immediately below the next integer cover.
- 119,088 noninteger cases, including valid real budgets whose floored products are below the real threshold.
- 68,346 synchronized rectangle controls retaining the auxiliary readout and coordinatewise bound.
- 1,960 actual J crossings satisfying the letter-specific cutoff.
- Sharpness family checked for k=1..32, plus an explicit k=0 regression.
- Literal J-sharp graphs for k=1,2 have respectively 57 and 165 Hamilton paths and exactly one Hamilton cycle, on 8 and 10 vertices.

The all-order claims rely on the reviewed inequalities and cut proofs, not these finite checks.

## Reproduction

From this audit directory:

    python3 independent.py > check-normal.json
    python3 -O independent.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json INDEPENDENT_NORMAL.json
    sha256sum -c MANIFEST.sha256

A release status promotion should archive the reviewed companion note and record its exact changes. This approval neither alters the main packet's independent audit nor proves universal target coverage or globally optimal target-representation constants.
