# Reusable cyclic-response framework and a minimal actual hostile

This is separate from the two already frozen Session138 zero-rigidity packets.

`GENERAL_CRITERION.md` gives exact finite **sufficient**, not necessary, tests:

1. For an actual module difference d, entrywise nonnegativity of its exact
   tridiagonal cyclic-response K proves nonnegative response for all nonempty
   side tournaments; K11>0 gives the sharp universal lower bound
   `K11 H(A)H(B)`.
2. A one-crossing row certificate, a log-concave actual base profile v, and
   `v^T K v<0` prove a whole negative chamber `P_side=P_base R`, R nonzero
   monomial-nonnegative. This includes eventual ordinal-power families from
   a single supplied finite base certificate.

`certify.py` implements these tests without treating `NOT_CERTIFIED` as a
negative theorem. Literal masks are validated before the actual-input wrapper
computes profiles. Formal profile interfaces retain their realizability limit.

A strong seven-vertex same-card pair proves the row condition cannot be dropped:
its cyclic response goes from-36 at singleton sides to+828 at transitive-pair
sides (and+67932 at triangle sides), despite a log-concave base and nonnegative polynomial multipliers.
An independent literal computation proves module order7 is minimal: every
same-order profile-pair kernel through order6 is entrywise one-signed.
No actual non-log-concave tournament profile is claimed; the short
`(2,1,2)` log-concavity countercontrol is explicitly formal.

## Replay

    python3 -B framework_audit.py
    python3 -B -O framework_audit.py
    g++ -std=c++17 -O2 -Wall -Wextra -Werror minimality_literal.cpp -o /tmp/minimality_literal
    /tmp/minimality_literal

Only stdlib Python and a C++17 compiler are needed. The local `profile_tools.py`
uses backward endpoint and anchored unordered-partition DPs; the C++ program
independently enumerates23,717,429 literal permutations. Sources, complete
witnesses, and matching outputs are included. Executables and the retired
`exploratory/` discovery files are excluded from the curated release.

The independent Session138 reconciliation audit accepted the framework,
minimum-order implication, exact actual witness and boundary conditions.
The separate audit packet is `independent_zero_framework_later`. Curated
release bytes are fixed by `public_manifest.json` and `SHA256SUMS`; no
failed certificate is interpreted as a realizable sign reversal.
