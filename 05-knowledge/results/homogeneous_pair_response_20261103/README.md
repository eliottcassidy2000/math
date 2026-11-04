# Homogeneous pairs, two actual queries, and the three-state seed interval

**PROVED:** a homogeneous pair cannot reverse between equal full named
two-port responses, even if the other arcs change arbitrarily and even up
to a positive scalar multiple. Hence no module containing such a pair
supports a hidden converse in any ambient tournament.

This exclusion was already proved by **C2** in
[the incoming twin-sign proof](../two_port_eight_core_20261103/algebra/PROOF.md).
The present [P1-P3](PROOF.md) recover its mechanism and give a direct
two-query certificate, with no ordered-cover reconstruction or positive-cut
theorem required. This is a refinement and integration, not a priority claim.

For twins a->b, contracting them to z in Q gives

    F_T(e_b,0)-F_T(e_a,0)
      = 2(E_z(Q)+sum_w C_zw(Q)) >= 2H(Q) > 0.

The same named difference is negative for reversed twins. Homogeneity in
both tournaments is essential. This separates an actual full response;
it does not identify that response with H,S,E,C.

The first independent controls cover all 5,405 rooted expansions of all
1,099 labelled quotients of orders 1 through 5. The endpoints-only hostile
has identical H,S,E after reversal but query vectors (5,9) and (9,5).
See [check.json](check.json), [check.py](check.py) and its matching optimized
receipt. Finite controls corroborate the all-order combinatorial proof.

The incoming three strong seed states have
H(Q_j[TT_r])=4*3^r+9*2^r-3+2j. Their scalar path-cover law is converse
invariant, while the named pair test exposes reversal for every r>=2.
Two nonuniform ports can produce five module runs; the scalar three-run
truncation is not a full-response truncation. [The integration](INTEGRATION.md)
credits the incoming scalar and cover-sign formulas, proves the sharp
five-run bound, and retains [complete named response polynomials](response_tables.json)
for all three seed states at orders7 and12. Their whole-converse contrasts
are exactly cubic; the partial-module quartic hostile is preserved separately.

The [response receipt](response.json) records42,668 actual compiled queries,
24 independent forward-DP controls and three literal word checks.
Both engines and the complete tables match their optimized replays.
Run `python verify.py` to check frozen hashes and source pins; `--replay`
reruns both engines normally and under -O into a temporary directory.
The [review](REVIEW.md) records quantifiers, failure boundaries and provenance;
the [concept board](BOARD.md) connects the three research lanes.

The twelve-vertex hidden-converse residual is therefore strong, prime and
nonselfconverse, using the inherited finite exclusions and **C4** of the
incoming proof. The prime branch and universal reconstruction remain
**OPEN**. No exact full-response collision is claimed.

