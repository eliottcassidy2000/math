# Independent audit: three-state three-run unit gate

**Verdict: PASS, within the explicitly stated scope.** No correction required. Audited on2026-11-04. No global completeness, new globally attained counts, unrestricted ternary return construction, or historical novelty is established.

## Independently reproduced finite data

The audit script imports no author code. It uses the displayed adjacency masks directly, verifies their common directed cycle and strongness with Floyd–Warshall reachability, and explicitly identifies the sole changed unordered pair on each edge: {2,4}, then {0,1}. Both pairs lie outside the replaced vertex5.

For each of the120 orders of the other five names, the checker constructs the admissible-gap set G and compulsory-gap set B, then uses the binomial formula in the proof. This independently checks every entry of the five-row grouped profile table, not just its column totals. The profiles are exactly(27,21,4), (29,21,4), (31,21,4). The seed's only out-neighbors are0 and1, so no Hamilton path can have four separate module runs: at least three nonterminal runs would need three distinct such followers.

A separate substitution builder and direct depth-first enumeration of directed paths reproduce all21 inflated counts for r=1,...,7, independently of the author's endpoint subset DP. These include the stated first six intervals and{9897,9899,9901} at r=7. Reachability verifies every resulting graph is strong.

For every labelled module tournament of orders1 through4, independent permutation-and-cut enumeration counts ordered path covers q2 and q3. Direct Hamilton-path enumeration of each substitution checks the arbitrary-module formula in all225 cases. This is independent of the author's color-assignment cover checker. All substituted graphs are strong.

Both author's normal and optimized outputs agree with both saved receipts. Both independent-audit runs agree as well. Audit evidence: `audit.py`, `receipt.json`, `normal.out`, `optimized.out`.

## All-order mathematical review

The occupied-gap bijection is correctly inherited from Session176 `spectrum/PROOF.md`, sections1–2, inspected during this audit. It applies to arbitrary module size: a nonempty gap receives one run, and each module path uniquely determines its other-vertex permutation, occupied gaps, and ordered cover. For transitive modules, the ordered cover counts are k!S(r,k). Thus the verified finite skeleton table proves for every r≥1

    H(Qj(r))=(27+2j)+21(2^r−2)+4(3^r−3·2^r+3)
             =4·3^r+9·2^r−3+2j.

The three-run coefficient remains nonzero, and the two higher channels balance across both literal edge reversals. General substitution yields(27+2j)H(A)+21q2(A)+4q3(A), so both count changes are exactly2H(A). The claimed unit-gate rigidity follows from H(A)=1 iff A is transitive; the proof explicitly supplies the triangle/insertion argument for the nontransitive direction.

Strongness of a nonempty substitution into a vertex of a strong quotient is valid even when the substituted module itself is not strong: a directed quotient cycle through the replaced vertex lets one leave and re-enter at any requested module vertex. Quotient paths handle distinct modules.

The recurrence obstruction is exact. A fixed shift identity for all sufficiently large r would force3^q=3 and2^q=3 by independence of the two exponential sequences and the constant sequence. No integer q satisfies both. The separate one-step difference is−9·2^r+6+2l−6j, whose maximum for r≥1 and j,l∈{0,1,2} is−8. These are restrictions on this particular three-state family, as stated, not on other constructions.

## The 7-adic claim

Independent exact modular exponentiation verifies2^42≡3^42≡1 mod49. The upper states at r4 and r6 have residues28 and14 mod49, respectively. Hence each corresponding positive progression modulo42 gives exact valuation1. Their starting counts469=7·67 and3493=7·499 have prime cofactors, checked by complete trial division through their square roots. No assertion of infinitely many prime cofactors is made. Since H grows exponentially with r and only three states occur at each r, O(log X) outputs below X is valid.

## Scope and inspected versions

The new claim is an explicit persistent three-point unit interval with a genuine three-run channel, not the general inflation machinery or arbitrary-feedback circuit already present in Session176. The finite exploratory searches are correctly excluded from unbounded nonexistence evidence. The scope caveats are sufficient.

SHA256 of inspected author files:

- PROOF.md:461ecf1b5498542fde54bee1faeaf23ce941a8a5dc942d43820dbdbc2c166063
- verify.py:5c48724c39bbed78a5c70f6fd83fa3fd1776cbf75033938c4dbb7287557327c5
- receipt.json:82262cb48ef8bbce560c624c7aa50253e4900322c2b4d24a7b0a8bd3e0d8213b
