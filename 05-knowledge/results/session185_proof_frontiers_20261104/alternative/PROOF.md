# A three-state unit edge gate with a genuine three-run channel

2026-11-04. **PROVED:** an explicit strong three-state interval at every transitive-module size, with two literal unit edge reversals and a nonzero three-run channel. **NOT PROVED:** spectrum completeness, reusable ternary digits, or new globally attained counts. Global exclusions 7 and 21 remain unchanged.

## Relation to inherited work

Session176 `spectrum/PROOF.md`, sections 1–2, supplies the general occupied-gap inflation identity and a two-state dyadic gate with profiles (11,6) and (13,6). Its section 4 already supplies the full arbitrary-feedback streaming circuit; this packet does not claim that circuit as new. The present explicit seed chain retains a nonzero three-run response while balancing both higher-run channels across two successive single-edge changes. It produces a three-point interval, not just one unit pair. Historical novelty is not established.

## Explicit construction

On names 0,1,2,3,4,5, use outgoing-neighbor masks

    Q0 = (16,5,33,39,46,3),
    Q1 = (16,5,49,39,42,3),
    Q2 = (18,4,49,39,42,3).

Q0→Q1 reverses only 4→2; Q1→Q2 reverses only 1→0. All three have the directed cycle

    0→4→3→1→2→5→0.

The common vertex 5 beats precisely 0 and 1. Replace vertex 5 by a transitive module X_r of size r≥1, with all external edges inherited from vertex 5. Denote the resulting tournament Qj(r).

**Theorem.** Each Qj(r) is strong, its order is r+5, and

    H(Qj(r)) = 4·3^r + 9·2^r − 3 + 2j,     j=0,1,2.

Both transitions remain literal single-edge reversals outside the inflated module. Thus the three values form a full three-point odd interval, connected by unit count changes at every r. The first intervals are {27,29,31}, {69,71,73}, {177,179,181}, {465,467,469}, {1257,1259,1261}, and {3489,3491,3493}.

## Exact finite skeleton certificate and all-order proof

Delete vertex 5 and fix an order p of the other five vertices. Let G(p) be the allowable gaps for a run of vertex-5 copies: its preceding vertex, if present, must beat 5, and 5 must beat its following vertex, if present. Let B(p) be the internal gaps at which the displayed direct edge of p points backwards. These gaps must be occupied. If B is not a subset of G this order contributes zero. Otherwise its coefficient for k occupied gaps is

    binomial(|G|−|B|, k−|B|),

interpreted as zero outside its ordinary range. This is a complete finite formula, with 120 orders for each seed. Grouping these coefficients by p's first vertex gives the following triples (a1,a2,a3):

| First vertex | Q0 | Q1 | Q2 |
|---|---|---|---|
| 0 | (4,4,1) | (4,4,1) | (6,5,1) |
| 1 | (5,4,1) | (7,5,1) | (5,4,1) |
| 2 | (3,2,0) | (4,2,0) | (3,1,0) |
| 3 | (12,10,2) | (11,9,2) | (12,9,2) |
| 4 | (3,1,0) | (3,1,0) | (5,2,0) |
| Total | (27,21,4) | (29,21,4) | (31,21,4) |

There can be no fourth run: each nonterminal run must be followed by one of the distinct vertices 0 and 1, the only out-neighbors of the replaced vertex. The inherited occupied-gap bijection therefore yields

    H(Qj(r)) = (27+2j) + 21(2^r−2) + 4(3^r−3·2^r+3).

Here k!S(r,k) counts assignments of the linearly ordered module vertices onto k ordered nonempty runs. Simplifying proves the formula for every r, not only the checked finite range. The common directed cycle proves strongness of the seeds. Inflating a vertex of a strong tournament by any nonempty tournament preserves strongness, since paths through other modules connect any pair of vertices in one module and paths between quotient vertices lift. This proves strongness for every r.

## General module law and its precise limitation

For any nonempty tournament A, let h=H(A), and let q_k(A) count ordered spanning covers by k nonempty directed paths. The identical skeleton proof gives

    H(Qj[A]) = (27+2j)h + 21q_2(A) + 4q_3(A).

Both single-edge transitions have difference 2h. Thus they are unit gates exactly when A is transitive. The converse uses the usual injection from the three directed-triangle Hamilton paths into full paths under successive deterministic vertex insertion: every nontransitive tournament has H≥3. This channel cannot be iterated on a nontransitive substituted gadget while retaining unit steps.

Nor does switching among these three gate states repair ternary digit closure. If, for fixed j,l and integer clock shift q,

    H(Ql(r+q)) = 3H(Qj(r)) + d

held for all sufficiently large r, independence of 3^r, 2^r and 1 would require both 3^q=3 and 2^q=3. There is no such integer q. In particular the single-clock increment has

    H(Ql(r+1))−3H(Qj(r)) = −9·2^r + 6 + 2l − 6j ≤ −8.

This exact obstruction applies only to the displayed three-state carrier, not to other interfaces or arbitrary tournament transformations.

## A bounded contribution to the 7p lane

For r≡4 (mod 42), the upper state has H≡469≡28 (mod 49). For r≡6 (mod 42), it has H≡3493≡14 (mod 49). Indeed 3^42≡1 and 2^42≡1 modulo 49. Consequently both classes yield infinitely many strong counts of exact 7-adic valuation one. Their first cofactors are 67 and 499, both prime. No infinite prime-cofactor claim is made. The family has only O(log X) outputs below X, so it cannot itself cover either unresolved prime lane.

## Reproduction and evidence boundaries

`verify.py` is standalone standard-library code. It independently enumerates literal permutation/gap patterns for the three profiles, checks inflated Hamilton counts by endpoint subset DP for r=1,...,7, verifies strongness, and checks the arbitrary-module identity for all labelled tournaments of orders 1,...,4 (225 substitutions). `receipt.json` and `receipt.optimized.json` agree; validation uses explicit exception gates, not removable asserts.

The exploratory files are not evidence of any unbounded nonexistence theorem. `explore.py` implements the already-proved Session176 two-feedback circuit, `search_returns.py` and `pure3.py` are bounded searches, and `gate_enum.py`/`larger_gates.py` enumerate small marked seeds. They motivated this explicit positive gate. Their empty return searches are not extrapolated to the full two-feedback family.
