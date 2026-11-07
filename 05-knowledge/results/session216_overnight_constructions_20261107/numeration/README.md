# Three-core numeration lane

Read PROOF.md and audit/AUDIT.md. This portable packet has been independently
audited. The all-odd-at-least-65 inverse conjecture remains open.
This packet proves a unique balanced-prefix candidate lemma, implements exact
all-length membership via deterministic alpha stripping, and records a
restricted curvature obstruction and explicit hostile tests.

Run from this directory:

    g++ -O3 -std=c++17 balanced.cpp -o /tmp/session216-balanced
    /tmp/session216-balanced CORE_STATES.txt 0 range 10000
    /tmp/session216-balanced CORE_STATES.txt 1 range 10000
    /tmp/session216-balanced CORE_STATES.txt 2 range 10000
    /tmp/session216-balanced CORE_STATES.txt 1 1000000007 1000000033 1000000000039 1099511627775 1099511627777
    /tmp/session216-balanced CORE_STATES.txt 2 1000000000039 4295971 9948407
    python3 inverse_alpha.py
    python3 verify.py

A balanced result with core=-1 is only a failure of the tested band, including
its explicit b<27 resource ceiling. It is never an all-length absence claim.
Only inverse_alpha.py can return ABSENT_FROM_GRAMMAR after its complete finite
inverse exhaustion. No background run or uncontrolled resource search is
needed for these receipts.

CORE_STATES.txt is copied verbatim from the audited arithmetic three-core
packet. The producer uses the f recurrence; verify.py uses its independent
selected-position pair sum and also performs 320 literal subset-DP graph
counts. LARGE_GRAPHS.json retains all seven displayed large-target graphs
and explicit Hamilton cycles proving their strongness. Their large counts
are certified by the independently evaluated inherited construction formula,
not by infeasible whole-graph subset enumeration.

## Arithmetic domain and audit reproduction

The C++ evaluator requires all signed-64-bit intermediate values to fit. Its
__int128 threshold check does not widen its formula arithmetic. The b<27
resource ceiling alone therefore does not establish correctness for arbitrary
numeric inputs. All retained receipts are within the safe arithmetic domain.
The Python all-length membership decision uses arbitrary-precision integers;
its completeness is not an all-target occupancy theorem.

From this directory, additionally run:

    PYTHONDONTWRITEBYTECODE=1 python3 audit/check_independent.py
    sha256sum -c SHA256SUMS

The independent check compares the alpha-stripping membership procedure with
a separate bounded direct enumeration and returning coordinate compiler.
Large graph counts use the audited formula/compiler, not whole-graph subset
DP; the 320 subset-DP controls are small graphs.

original/ preserves the exact pre-release inputs and original independent
checker. RELEASE.diff gives all changes to proof, README and checker. The
checker change only relocates its input directory for this portable layout.
The historical audit and input manifest are preserved unchanged in audit/.
PUBLICATION_FILES.txt lists the portable files; SHA256SUMS hashes every listed
file other than itself. No cache files are part of this release.

The audited three-core proof dependency and its final independent verdict
are included under dependencies/, with original proof and manifest hashes.
Its historical pending-review wording is preserved unchanged. The supplied
sibling session216-three-core-portable/ contains the complete certificate.
The Session203 curvature citation records inheritance; this packet claims
only its explicitly restricted fixed-block curvature result.
