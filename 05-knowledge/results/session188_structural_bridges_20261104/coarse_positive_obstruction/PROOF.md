# Positivity alone does not prevent complete coarse blindness

**Exact formal-profile obstruction, not a tournament collision.** There are three finite sequences q^A,q^B,q^C with q_0=0, q_1=1, and every entry on their positive support strictly positive, for which all 64 block-uniform two-port contrasts of the cyclic quotient versus its converse vanish in the ordered-cover profile functional.

The sequences are not realized by tournaments. This result rules out a proof of universal unequal-child coarse separation based only on finite support and positivity of the three cover sequences. It does not settle the realized-tournament problem, a full named-response collision, or any ambient question.

## 1. Exact data

The first two sequences, displayed from q_0 through q_10, are

    A=(0,1,8,2,5,2,9,5,4,7,2),
    B=(0,1,9,1,8,9,3,5,6,1,8).

The third has q_0=0, q_1=1 and positive rational support through q_13. Its complete exact rational entries are in `profile_data.json`; none is rounded or omitted. The same file records the positive common denominator D and the complete integer sequence Dq^C.

Multiplying one child's whole profile by D multiplies every query value by D, since the profile functional is linear separately in each child. Thus clearing denominators also gives three strictly positive integer sequences with the same all-64 equality, though the third first entry then equals D rather than one.

The lack of tournament realizability is explicit, not merely unverified: every tournament's ordered q_k is divisible by k!, because permuting the k nonempty cover slots acts freely. Here q_2(B)=9 is odd and q_3(A)=2 is not divisible by six. The top profile entries also do not have the required tournament values. Consequently these data cannot be promoted to a counterexample involving actual tournaments.

## 2. Complete query certificate

Use the exact physical-run-word generating function from `coarse_response_rank/PROOF.md`, Section 1. With R=1-XYZ and cyclic run matrix N, each actual block mask pair (a,c) has a full numerator W_(a,c)(N). Its formal query value is

    L(W_(a,c)(N)/R^3),
    L(X^rY^sZ^t)=q_r(A)q_s(B)q_t(C).

Transpose N for the converse quotient. The verifier evaluates both functions at every one of the 64 masks by exact rational/integer arithmetic. All pairs are equal and positive. Both values, using both the normalized rational third profile and its denominator-cleared integer version, are retained for every mask. The certificate includes all port placements and repeated child runs; no quadratic-only or aggregate-start approximation is used.

The ten-coordinate theorem in the adjacent rank packet would suffice to extend a ten-query equality to all 64, but this verifier computes all 64 directly rather than relying only on that implication.

## 3. Independent run-word verification

A separate automaton in this packet counts every admissible quotient word with both distinct ports used exactly once, each child visited between one and its declared support cap times, and no equal consecutive child colors. Port transitions use the selected actual block masks and x->y. Each terminal word with child visit counts r,s,t receives the exact product q_r(A)q_s(B)q_t(C). Terminal words contribute and can continue, so all possible run counts are included.

This automaton does not use the rational numerator formula. It independently checks both orientations at all 64 masks, giving 128 exact formal word counts. They agree with the numerator functional and each other. These are formal weighted word counts, not literal Hamilton counts of a realizing tournament.

## 4. Scope and reproduction

The obstruction shows why the equal-profile positive-square theorem cannot be generalized to three unrelated positive sequences merely by dropping the equality hypothesis. Additional restrictions satisfied by actual tournament profiles are necessary for such a proof strategy. The certificate does not show that those restrictions are sufficient, and does not exhibit a genuine tournament collision.

Run:

    python3 verify.py > normal.json
    python3 -O verify.py > optimized.json
    cmp normal.json optimized.json

Only the Python standard library and the adjacent `coarse_response_rank/verify.py` are used. The complete frozen source packet consists of `PROOF.md`, `verify.py`, `profile_data.json`, `normal.json`, and `optimized.json`. The source and dependency hashes are retained in both receipts. Explicit checks remain active under optimization.
