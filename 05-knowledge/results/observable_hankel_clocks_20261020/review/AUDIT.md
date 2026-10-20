# Independent audit: observable Hankel clocks

2026-10-20 UTC. **PASS for the mathematical claims reviewed.**

Reviewed final target: `session9-hankel/OBSERVABLE_HANKEL_CLOCKS.md`, SHA-256 `268e5ec48f7174b39579f514f533af06d6c25cb2e8fa9f271acdd010d615b358`. This final revision includes the staircase proved below, the torsor clarification, companion intertwining, and the distinction between reduction quotients.

Scope: the observable-image theorem, exact integer sequence, full Smith lists, binary period, six-register precision boundary, compatible orbit completion, and full state-module completion. The source paper's own analytic/coloring arguments are outside this audit. No physics claim is inferred.

The independent program reconstructs the sequence and matrices directly. It imports no root checker. Its symbolic lane identities, integer Bareiss determinant, Smith computation, prime-field elimination, and modular recurrence checks agree. Normal and optimized outputs are byte-identical. The final root checker was separately inspected and replayed; its output matches its frozen receipt, and all full and rectangular-prefix local-pivot lists agree with the independent construction.

## 1. Observable image theorem

The proof is sound over every residue ring R=Z/m, including composite m. Let χ be a globally valid monic degree-r annihilator. The coefficient map R^r→R^N sends c to Σc_jE^j a. First-r-value evaluation is Hc. A sequence in this shift span that has r initial zeros vanishes globally, because it obeys the same monic recurrence. Hence the image is exactly the sequence module M, with no field or freeness assumption.

To make the stability statement explicit, let C be multiplication by x modulo χ in the coefficient basis, and let S advance a first-r-value state. The exact identity is

\[
HC=SH.
\]

Thus ker H is C-stable. This is the precise finite quotient; an arbitrary low-dimensional measurement list would not suffice.

The invertibility claim correctly uses χ(0)∈R*. If χ(0) is not a unit, the note says only that the inference fails, not that tails are necessary. That distinction is correct: a nonminimal annihilator may have nonunit constant term even for an invertible realized shift.

The phrase “smallest faithful linear state module” is justified in the observable/reachable sense: any linear realization of a restricts to the shift-generated reachable submodule; quotienting its states with identically zero future output gives M. This says nothing about minimum nonlinear encodings of one orbit.

## 2. Independent exact data

For the normalized three polynomial lanes, symbolic expansion verifies χ(E)a=0 identically in each residue lane. This is a global identity, not a finite-prefix fit. An independent integer Bareiss elimination gives

\[
\det H=-2^{206}3^{191}5^{22}7^8 61^7\ne0.
\]

The nonzero 22×22 Hankel determinant also independently certifies that no rational recurrence of order below 22 exists. All five stated Smith valuation lists agree with the independent calculation. Separate modular row elimination gives ranks 6,3,14,21,21 at p=2,3,5,7,61 respectively. Good-prime controls p=11,13,17,19 have rank 22.

The program checks HC=SH exactly and det S=1. Therefore the whole integer shift is unimodular; the observed finite quotients retain an invertible action.

The order-six candidates are verified exactly:

- Modulo 64: (7,48,16,2,16,48); all first 22 residuals vanish.
- Modulo 128: (7,112,16,66,80,48); the first failures are residual 64 at n=6 and n=7.

The leading 6×6 Hankel minor is invertible modulo 128, so the latter candidate is unique. The first 22 vanishing residuals in the mod-64 case propagate globally because the residual itself obeys χ. Consequently these certificates prove the stated sharp boundary, rather than just suggesting it numerically.

## 3. Period and completion

The proof P_k=3·2^k is valid for every k≥1. Reduction modulo two forces any period to be divisible by six. The residue-two polynomial has the form −X+4Q(X), so for any such candidate period m its difference has valuation exactly v₂(m). The upper bound follows from integer-polynomial congruences and preservation of the residue lane.

An independent finite check verifies the least period for k=1,…,12 and rejects every proper divisor of 3·2^k. The all-k conclusion rests on the proof above, not on this range.

The group of shifts acting on M_(2^k) has exactly this order: equality on the original sequence implies equality on all its shifts and their linear span; conversely equality on all states includes the original sequence. The reduction maps preserve phase labels. Thus the group completion is

\[
\varprojlim_k\mathbb Z/(3\cdot2^k)\cong\mathbb Z/3\times\mathbb Z_2,
\]

with supernatural order 3·2^∞. The factor at 3 is exactly one, not ∞.

The revised note now correctly distinguishes this group from the orbit itself. The orbit closure is its torsor, identified with the group only after a base state is chosen. Its map into the additive state lattice is not a group homomorphism. In particular, the additive lattice has no order-three torsion and the orbit never contains zero, already visible at precision one.

The inverse limit of the observed images really is L=HZ₂^22. For example, each finite compatible observed state defines a point of Z₂^22 that is congruent to an element of L at every precision; L is compact and hence closed, so that point lies in L. Conversely every element of L gives compatible observed states. The index is 2^206 and the abstract Z₂-rank is 22.

The equivariant isomorphism H:Z₂^22→L supplies the completed companion action. Its finite observation filtration is

\[
I_k=H^{-1}(2^k\mathbb Z_2^{22}),
\qquad M_{2^k}\cong\mathbb Z_2^{22}/I_k,
\]

which need not be the standard filtration 2^kZ₂^22. This revision correctly identifies the nontrivial finite-precision question.

“Larger” here concerns finite observable state counts and linear module structure. Both the orbit closure and L are uncountable zero-dimensional compact spaces; raw cardinality is not the distinction. As dynamical systems they are different: the whole linear state space contains the fixed point zero, whereas the orbit torsor's adding-one action is minimal and has no fixed point.

## 4. Verified extension: the complete binary register-order staircase

Let r(k) denote the smallest order of a monic scalar recurrence for a modulo 2^k. The following strengthens the order-six boundary:

\[
r(k)=\begin{cases}
6,&1\le k\le6,\\
12,&7\le k\le14,\\
14,&k=15,\\
18,&16\le k\le20,\\
20,&k=21,\\
22,&k\ge22.
\end{cases}
\]

**Status:** PROVED conditional on the finite-exact prefix Smith invariants below, independently verified by the accompanying computation. No all-k extrapolation is used.

Write v=(0⁶,6⁶,14²,15⁴,20²,21²). For each r in {6,12,14,18,20,22}, the 22×r matrix H_r consisting of the first r columns has 2-adic Smith valuations exactly the first r entries of v. The full integer diagonals are reproducible from the program; the needed valuation certificates are frozen in the JSON.

**Lower bound.** The additive decomposition in the note implies

\[
\dim_{\mathbb F_2}(M_{2^k}/2M_{2^k})=\#\{i:v_i<k\}=:d(k).
\]

A monic recurrence of order r makes the first r shifts generate M_(2^k), so their images generate this vector-space quotient. Hence r≥d(k).

**Upper bound.** Set r=d(k), which is one of the six checked prefix lengths. The image of H_r is contained in that of H. Their logarithmic sizes are respectively

\[
\sum_{i\le r}\max(k-v_i,0),\qquad
\sum_{i\le22}\max(k-v_i,0).
\]

The omitted entries satisfy v_i≥k, so these sizes coincide. Equality of finite contained images follows. Thus the first r shifts span M_(2^k), and E^r a is a linear combination of them, giving a global monic order-r recurrence. Therefore r=d(k). The finite list of v_i proves the displayed formula for every k, including all k≥22.

These recurrences can be chosen with odd constant feedback coefficient. The first r generators induce a basis of M/2M, and E induces an invertible companion action on that quotient. Its determinant is nonzero, forcing the coefficient of a_n to be a unit. Thus the staircase is achievable by invertible binary-modular register models, not just noninvertible scalar recurrences.

**Important quotient distinction:** at k=7, M_128/2M_128 has dimension 12. This is not the same as reducing its observed-coordinate vectors modulo two, whose image is M_2 and has dimension 6. Confusing these operations would erase the very precision layer responsible for the recurrence-order jump.

## 5. Deliverables and boundaries

- `independent_hankel_check.py`: independent construction and validation.
- `independent_check.json`: exact normal-run receipt.
- `independent_check.optimized.json`: optimized-run receipt, byte-identical.
- `REPLAY_RECEIPT.json`: commands, dependency versions, hashes, target revision.
- `SHA256SUMS`: integrity manifest.

No edits were made to the repository checkout or target note by this reviewer. The audit does not assert historical novelty, formal Lean verification, or any transfer to the source paper's analytic theorem or to a sphaleron dynamics claim.
