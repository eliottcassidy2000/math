# New operations, actual geometry, and all-row debt

**Session254, 11 November 2026. COMPLETE RESEARCH SYNTHESIS.**
This session merges prior planar-Jacobian, Hadwiger–Nelson and Lonely Runner work with Chan–Pak through actual operations and their retained information. The output includes new proofs in the notes, exact symbolic certificates, and bounded controls. Inherited claims were used at their stated scope, not subjected to a new general audit. No literature-priority claim or canon promotion is made. **Chan–Pak, JC(2), general LRC and the targeted plane-colouring improvement remain open.**

## 1. Chan–Pak: a genuine count-preserving repair, not just a new obstruction

The path-colour polynomial Q carries the complete prime-power substitution response. Adding r universal sinks multiplies Q by t^r, while preserving the exact Hamilton count H=Q(1) and every old vertex and arc.

For 3-unit sources this legal operation now has exact normal forms:

- Modulo9, at most7 vertices erase the extra derivative state. Only H and chosen order parity remain.
- Modulo27, at most20 vertices normalize the first derivative, but one curvature class K=2HE−D²+HD modulo3 survives, with D=Q′(1), E=Q″(1)/2.
- Modulo81, at most57 vertices normalize the correct observable derivative, leaving one nine-valued coordinate. The observable is D+9F modulo27, with F the third Hasse derivative: integer-function aliasing u³≡u modulo3 prevents treating raw derivatives as independent measurements.

These are normal forms for EVERY iterated substitution context modulo the specified power, not only a small set of test quotients. Padding can increase a later output enormously; the precise two-ear cost H+2^rQ(2) remains attached to the state.

The constructive step is stronger than the normalization. Two strong six-vertex H37 sources have different curvature. Reversing one specified arc in an explicit H37 graph changes K from2 to0 while preserving all vertices and H. Its polynomial change is exactly

    2t²(1−t²).

It vanishes at1 but changes the derivative state. This supplies an actual positive-fibre surgery rather than a signed formal coefficient combination. The operation extends through ordinal joins, producing an infinite count-preserving family at counts37H(R) for3-unit R. It does not claim every realization of those counts contains the patch, or preserve every marked endpoint needed by a ternary constructor.

[Proof, explicit graphs and operation cost](chan_pak/PROOF.md).
[Modulo81 and the precise derivative alias](chan_pak/NEXT81.md).

**Next constructive target:** transport this neutral-arc repair into a five-core source while retaining its endpoint/readiness data. The H141 source remains outside these3-unit normalizations; it needs its own nonunit analysis.

## 2. Planar Jacobian: the homogeneous interface is now known at every row

The initial task was a row16 quotient. Its solution exposed a structural semigroup theorem. In the inherited source-normal chart

    u=x²t, p=t(1+x²t), y=xtp,
    P_d=sum_(a+b≤d) x^a u^b K[p,y],

the exact leading coefficient space at t-row n, after all earlier rows vanish, is

    K[x] of degree at most floor((n+3d)/2), for n≥d.

The proof retains two different vanishing points: u=0 and u=−1. This gives a sharp degree budget that cancellation cannot evade. Individual monomials attain every permitted leading exponent. It is valid over every field; the Jacobian application remains characteristic zero.

Consequently, the homogeneous terminal tangent at every row n≥3 has degree at most floor(n/2)+2. For the next bracket row m≥5, the complete quotient has dimension floor(m/2)−1. Only ONE of those directions is reached by the entire prefix-preserving same-row Hamiltonian source packet. The other floor(m/2)−2 coefficients are rigid under that packet.

At row16 the exact seven-coordinate ledger consists of six high-degree coefficients x^11,...,x^16 plus one Student moment. The source response to the latter is

    4495s0+870s2+540s4+600s6+1008s8.

The least-weight visible scalar payer has weight2m−2floor(m/4) at every m≥5. Thus rows4k−1 and4k share weight6k; rows15/16 both admit weight24 scalar payments. This does NOT solve their other debts.

The all-row formula also defeats a tempting extrapolation: row15/16 happen to fit degree n−6, but row17 has maximum degree10, not11. The mechanism, rather than the short pattern, supplies the continuation.

[All-row theorem and proof](jacobian/ALL_ROW_DEPTH.md).
[Complete row16 response and remaining actual-source obligations](jacobian/PROOF.md).

**Next native target:** evaluate the six explicit rigid row16 debts on the actual inherited G_m family, then the inhomogeneous depth condition. Existence of a compatible prefix, actual row16 lifting, source entry, and polynomial termination remain unproved here. The new theorem makes those missing conditions explicit at every scale; it does not make them vanish.

## 3. Lonely Runner: a strict, phase-uniform seven-pack that attaches physically

A nonnegative integer occupancy C satisfies

    E|C−1|=2Pr(C=0)+EC−1.

A zero-mean bounded observable therefore converts one uncancelled Fourier coefficient into a positive empty-occupancy set. For a shifted runner comb, a frequency receives contributions only from slopes dividing it.

If one of seven positive integer slopes occurs once and no other slope divides it, that frequency has a single contributor. Arbitrary phases cannot cancel it. The resulting seven-pack has safe measure greater than311/18750 even at radius2/25, strictly larger than1/14. Other slopes may repeat, so the theorem includes (1,2,2,2,2,2,2).

The strict margin pays for affine perturbations. Let six core speeds have maximum M, and write seven further speeds as h a_i+b_i with B=max|b_i| and the isolated-slope condition. If

    h≥7M,    3h≥175B,

and all thirteen physical speeds are positive and distinct, LRC(14) holds for this family. The proof uses the cited lower-dimensional theorem only to obtain a safe interval for the six-speed core. It then gives an actual simultaneous safe set inside EVERY named block of length1/h contained in that interval, with measure greater than311/(18750h).

For example, it covers the tightly clustered set

    {1,2,3,4,5,6,175,347,348,349,351,352,353}.

Seven equal slopes with equally spaced phases show why phase and multiplicity cannot be discarded: at radius1/14 only isolated safe points remain, and every larger radius fails. There is no all-entry or arbitrary thirteen-speed theorem here.

[Fourier proof, physical attachment and exact interval controls](lonely_runner/PROOF.md).

**Next native target:** allow several contributors to the lowest divisor frequency while using the actual offset-phase path to prevent cancellation somewhere inside a protected core interval. This is a specific new obligation, not another scalar density comparison.

## 4. Hadwiger–Nelson: the geometric multiplication exists, but support still matters

The old two-common-neighbour cap obstructed unlimited length-two parallel paths. Longer paths are different. There is now an explicit algebraic construction of arbitrarily many internally disjoint length-four paths between two fixed ports at distance TWO, with no unwanted unit-distance edges. The entire point set's induced unit-distance graph is exactly that parallel-path graph. Both equal and unequal port colourings are genuine geometric states.

At six colours, its exact conditional channels are

    (105^r,104^r).

Thus geometrically lawful parallel multiplication can make unequal ports arbitrarily rare. Yet104^r is positive for every finite r. The family is bipartite, not a new high-chromatic graph. This separates a real geometric advance from the support obligation that amplification cannot discharge.

The proof uses rational-circle parameters and one square root per path, and excludes every unwanted cross edge analytically. The strengthened distance-two construction avoids the misleading variant in which a unit edge between ports must first be deleted.

[Explicit coordinates, all-pair incidence proof and exact counts](hadwiger_nelson/PROOF.md).

**Next native target:** controlled cross-path interactions that delete unequal states without deleting every equal state. Adding further independent paths cannot do that. Search should change the incidence, not merely amplify the same positive channels.

## 5. Laurent niche: a general factor theorem and four further full-matrix strips

The Session253 complementary-resultant frontier now has a universal forced-factor theorem and a sharper degree bound. At complementary width ell, the resultant is divisible by

    product_(j=1)^ell (r+j)^j(2r+2j+1)^j,

and its total degree is at most4ell². The residual degree is at most3ell²−ell. The factor divisibility follows from exact rank loss at two families of parameter values; it is not inferred from sampled factorizations.

Complete exact coefficient certificates show positive residual polynomials for ell4,5,6,7, of degrees44,70,102,140. Together with the inherited strips0–3 and small-block theorem, the FULL carried-response Smith partition is now known for all r≥2 at every widthell≤7:

    (0 repeated ell+1 times, 1 repeated r−1 times).

This is unbounded in r, not a bounded height census. Arbitrary width and physical positive-x noncancellation remain open.

[Proof and full rational certificates](niche/PROOF.md).

## 6. What genuinely transfers

Three distinct forms of transfer emerged, with their boundaries retained:

1. **Exact operator/gauge calculation.** Chan–Pak's curvature is (t d/dt)²log Q at1; sink padding changes only the first logarithmic derivative. The JC Student map has a fully computed annihilating moment and rigid quotient. The transferable procedure is to compute the annihilator of a lawful operation BEFORE inferring repairability. No direct identification of their state spaces or positive realizations is asserted.
2. **Exact incidence identity.** The LRC integer-occupancy formula transfers verbatim to finite incidence tables of actual Chan–Pak constructions. It supplies an empty-fibre obstruction under a biased observable, not positive target coverage. Legality, integer target and cost must stay in the table.
3. **Realizable positive multiplication.** The HN construction realizes the Cartesian product of conditional colouring sets at fixed geometric ports. Tournament ordinal join realizes product of path-colour polynomials. Both have genuine positive combinatorial product operations, but no bijection equating their native statistics is claimed. Their support losses differ and cannot be erased by sharing an algebraic formula.

The previous literal Laurent/Chan–Pak centralizer identification remains an inherited exact bridge. These new results do not turn it into a geometric embedding or a positive tournament constructor.

The practical outcome is a coordinated research program with new native theorems, not four conjectures collapsed into one. The best near-term cross-problem bet is **operation-sensitive repair**: retain the precise unpaid coordinate, then find a realizable local move that changes it at controlled integer cost. The H37 arc reversal demonstrates that this can produce an actual construction rather than merely another invariant.

## 7. Evidence and publication scope

Seven new-claim engines reproduce the two Chan–Pak prime-power packets, both JC packets, HN geometry/counts, LRC exact intervals, and Laurent symbolic certificates. The source trail is pinned in each lane and aggregated in SOURCE_PINS.json. The root replay checks only this session's new claims; it does not rerun inherited large censuses or audit outside conjecture claims. Exact proofs, symbolic identities, finite controls, and open obligations are labelled separately.

The starting main was7b7c5965cc1. The original shared dirty checkout was left untouched. All publication checkpoints stage only the session packet and its small navigation links; concurrent work is retained. See MANIFEST.json and SHA256SUMS for artifact integrity.
