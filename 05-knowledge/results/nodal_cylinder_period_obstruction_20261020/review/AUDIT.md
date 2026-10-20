# Independent audit: nodal-cylinder period obstruction

Date: 20 October 2026.

## Verdict

**PASS: the proof of Theorem 1 and the stated generic conductor-degree-two exclusion are sound.** The only requested wording clarification concerns ambient equivalence of curves/defining ideals, rather than exact equality of arbitrarily rescaled equations. The producer incorporated it in Section 2 before the reviewed hash was frozen. No mathematical repair is required.

The conclusion is about generating pairs of the entire algebra C[x]+(y²−1)C[x,y]. Combined with the actual-source conductor packet, it excludes generic source-conductor degree two under its monic, primitive, nonzero-constant-Jacobian hypotheses. It neither excludes an exceptional degree-two slice in a higher-generic-degree family nor proves the planar Jacobian conjecture.

Reviewed source paths and SHA-256 hashes are in `source_hashes.json`. All outputs written by this audit are confined to this review directory. The producer files and frozen conductor packet were not edited by the reviewer.

## Inheritance and audit board

Closest earlier proved mechanism: THM-3586, `math/01-canon/theorems/THM-3586-nodal-cylinder-cap38-width-period-and-second-conductor-keller-gates.md`, Section 4.3, gives the period derivative identity for arbitrary pairs in the nodal ring. Its arbitrary-pair PDE remains a different, weaker-hypothesis problem.

Canonical hostile: P=y²−1, Q=x+y(y²−1), which generates the ring and is not Keller. Corrected near miss: endpoint equality or a unit conormal determinant alone does not supply the ambient embedding needed to make the period constant. Least-used sidecar restored by the new note: global algebra generation, and consequently a global invertible conormal matrix.

The audit board was: normalization/conductor; relative generation versus fibers; abstract versus parametrized embedding; conormal orientation; polynomial closed forms; generic versus exceptional source degree. The independent endpoint-jet formula below provides a second representation of the same invariant, rather than an unrelated exploratory claim.

## 1. Exact algebra and the determinant

Write R=C[x], w=y²−1, z=yw. The identity

E=R+wR[y]=R[w,z] ≅ R[W,Z]/(Z²−W³−W²)

is correct. Every multiple of w is expressed by separating even and odd powers of y, and the displayed kernel is the irreducible relation of z over C(x,w). Evaluation at either endpoint has kernel I=(w,z)=wR[y].

If E=R[P,Q], subtracting the polynomial's evaluation at the endpoint proves I=(P−p,Q−q). This is an equality of ideals over R, not merely a pointwise claim. The relation has no linear term at (W,Z)=(0,0), so I/I² is free of rank two with ordered basis (w,z). Its two displayed generators therefore give a surjective square matrix L over R. A right inverse implies det L is a unit; the units of R are exactly C*. This justifies one global constant λ and does not require a relative ambient lifting theorem.

A useful independent description is obtained from endpoint derivatives. For f∈E,

- its w coefficient modulo I² is (f_y(1)−f_y(−1))/4;
- its z coefficient modulo I² is (f_y(1)+f_y(−1))/4.

Indeed, I² pulls back into (y²−1)²R[y], whose first derivatives vanish at both endpoints. Consequently, with v_±=(P_y,Q_y)(x,±1) as columns,

    det[v_+,v_-]=8 det L,
    A=(1/15)det[v_+,v_-].

The second identity here is equivalent to the audited period theorem; the first is elementary and supplied the independent computational conormal calculation.

## 2. The imported theorem has the needed quantifier

The fresh primary-source check in `PRIMARY_SOURCE_CHECK.md` confirms that the exact classical input applies to the nodal cubic and quantifies over all alternate embeddings.

The source is correctly labeled CITED, not reproved by the finite script. Direct PDF access did not succeed; indexed text of primary-paper pages 1061–1063 was independently retrieved. This access limitation is explicit in the source record.

The indexed text's equation-level equivalence notation should not be used to claim an arbitrary nonzero scalar multiple of an equation is the exact image under an automorphism. For example the finite nonzero critical value of h scales under multiplication of h. The relevant curve/ideal conclusion is sufficient and is now stated in that form by the producer.

## 3. Abstract equivalence really gives the prescribed parametrized lift

For each a∈C, the algebra equality specializes to C[P(a,y),Q(a,y)]=E₀. The resulting surjection C[U,V]→E₀ has a height-one prime kernel, hence one irreducible defining equation. Thus this is genuinely a closed embedding of the abstract node, and the cited theorem applies.

Ambient equivalence of embedded images alone does not prescribe which abstract isomorphism is used. The note correctly repairs this. The normalization of E₀ is C[y], with fraction field C(y), and every automorphism extends uniquely to normalization. The conductor is (y²−1)C[y]. To check it directly, both f and yf must have equal values at ±1, forcing both endpoint values of f to be zero; conversely that vanishing suffices for fC[y]⊂E₀.

An automorphism y↦by+d preserving the conductor must preserve the unordered pair {−1,1}, so d=0 and b=±1. Both automorphisms extend to (W,Z)↦(W,±Z). Composing the initial equivalence with the appropriate extension therefore yields precisely P(a,y),Q(a,y), not a parametrization differing by an uncontrolled sign.

Because the nodal relation has order two, the derivative at the node of this corrected ambient automorphism is exactly L(a). Representatives differing by a multiple of the relation have the same first-order coefficients. Its ambient Jacobian is constant, equal to its value at that point, namely λ. No rational, continuous, or polynomial choice of these automorphisms as a varies is assumed.

## 4. Closed-form calculation and the boundary sign

For each corrected ambient automorphism (F_a,G_a),

    d(F_a dG_a−λ W dZ)=(det DΦ_a−λ)dW∧dZ=0.

Every closed polynomial one-form on the plane is exact, by integrating one coefficient and then the residual coefficient in the other variable. Pulling back its potential gives equal endpoint values because both normalization endpoints map to the same node. The remaining model moment is

    ∫[-1,1](y²−1)(3y²−1)dy = 8/15.

Hence A(a)=8λ/15 for every a. Since A is a polynomial in x, this is the stated polynomial identity. This is a pointwise-to-polynomial argument, not a relative lifting assertion.

Differentiation and integration by parts give exactly

    A'=∫[-1,1](P_xQ_y−P_yQ_x)dy+[P Q_x]_-1^1.

Both P and Q_x have equal endpoint values, so the boundary term is zero. If the Jacobian belongs to C[x], its integral is twice that polynomial; it must therefore vanish. In particular it cannot be a nonzero constant. The signs, orientation, and factor 8/15 are correct.

## 5. Actual-source consequence

The frozen conductor packet's Theorem D supplies exactly what is needed, rather than only equal endpoint values. Under noninvertibility and generic conductor degree two:

1. The even, nonzero slice-degree rule forces every slice to have degree two; the leading coefficient is a nonzero constant.
2. Immersion makes every local branch smooth. Total delta one therefore gives precisely one transverse node, not a cusp; its two conductor roots are simple.
3. The quadratic discriminant has no complex zero, hence is a nonzero constant. Completing the square and applying a constant scaling gives conductor ideal (y²−1), by a source automorphism with constant nonzero Jacobian.
4. Endpoint equality gives S⊂E. The actual conductor identity supplies (y²−1)B⊂S, while R⊂S. This gives equality S=E, not merely membership of P,Q in E.

The transformed pair is still Keller, so Theorem 1 contradicts item 4. Generic source-conductor degree is even; degree zero gives collision-free slices and invokes the conductor packet's cited one-line theorem. Thus the displayed generic lower bound four and the equivalent generic defect lower bound two follow. No exceptional-slice lower bound two for delta is inferred.

This audit checked the use and local reasoning of the frozen conductor packet, not a replacement independent re-audit of every earlier dependency.

## 6. Independent exact controls

The producer's script was rerun with normal and optimized Python. It reproduced its stored output exactly: **285 checks over 17 declared generating families**.

The separate `independent_checks.py` does not import producer code. It uses coefficient moments

    ∫ y^k dy = 0 (k odd), 2/(k+1) (k even)

rather than symbolic integration, and computes conormals from endpoint derivatives rather than ambient linear coefficients. It checks **491 exact identities over 36 independently selected generating families**, including nonlinear x-dependent triangular maps, reversals, and constant scalings beyond ±1. Normal and optimized outputs are byte-identical.

Hostile controls separately verify:

- localization-only generation: (w,xz) has determinant x and nonconstant area;
- unit conormal without generation: (w,z+xw²z) has determinant one but area 8/15+64x/315, with the extra collision at x=−1, y=±√2;
- a nonmember boundary trap: P=y²,Q=xy has A'=2/3, Jacobian integral −4/3, and boundary term 2. Dropping the boundary term outside E is detected;
- adding multiples of the nodal relation changes neither pullback nor conormal first-order data.

These are FINITE-EXACT normalization and implementation checks. They do not prove generation for arbitrary pairs, uniqueness of all embeddings, or the conjecture.

Reproduce from the workspace root:

    python3 session13-nodal-review/independent_checks.py --output /tmp/nodal-independent.json
    python3 -O session13-nodal-review/independent_checks.py --output /tmp/nodal-independent-optimized.json
    cmp /tmp/nodal-independent.json /tmp/nodal-independent-optimized.json
    python3 session13-nodal-exclusion/verify_nodal_period.py --output /tmp/nodal-producer.json

## 7. Transfer boundary and review close-out

The mechanism is:

    global algebra generation
      → invertible conormal frame over C[x]
      → constant determinant
      + all-embedding uniqueness on each nodal fiber
      + correction of the finite node automorphism group
      → constant polynomial loop period
      → zero mean Jacobian.

Passing only to fiberwise embedded images loses the prescribed parametrization; the ±y correction restores it. Passing only to conormal jets loses global embedding; the nongenerator hostile shows why the entire algebra equality must be retained. Passing only to generic generation permits the determinant to vary and fail at an exceptional fiber. These are genuine failure boundaries, not proof omissions.

Recommended status: **PROVED from the specifically CITED classical embedding theorem**, with the generic conductor corollary conditional only on the explicitly stated actual-source hypotheses and its already audited conductor input. No priority claim for the classical embedding theorem or conductor identities is warranted. The broader arbitrary-pair nodal PDE and higher generic conductor degrees remain unclosed here.

Repository hygiene: the reviewer left `math` clean and made no publication changes. A final fetch found HEAD e34b69dcad1eff43cfef2a809e06a4c4bda84a21 zero commits ahead of origin/main and one commit behind; the incoming commit only reserves THM-4724 for unrelated affine/circulant work and adds an index entry. The parent owns publication of this review, as explicitly instructed.

Final status-only acknowledgment: changing the reviewed proof's single status line from “PROOF COMPLETE; independent review pending” to “PROVED; independently reviewed” gives final SHA-256 `dd0d343128a77f51fced21eee047558e0992d6fa94566e50b598771fd3d3a46a`. Reverting exactly that line recovers reviewed hash `3b5f9da355f961f5b2f705eed5bc7b37e161189e94edd0fa26685cbd6102c84a`; this byte-level test passed. The mathematical content is unchanged.
