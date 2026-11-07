# Audited fixed-H5 register and two-sided pair

Read REGISTER.md and REVIEW.md. Exact register occupancy is not Hamilton-count occupancy. The register theorem and persistent pair do not prove a tail coverage theorem. The four-choice exploratory menu remains outside this audit.

From this directory run:

    g++ -O2 -std=c++17 audit.cpp -o /tmp/register-audit
    /tmp/register-audit

Compare the JSON with RESULT.json. The copied helper three_core_audit.cpp makes all subset-DP and graph routines self-contained. All published files are at most131072 bytes. FROZEN_REGISTER.md and STATUS.diff preserve input and status-only promotion. Fixed finite integer checks fit signed64 arithmetic; the all-length claims are proved algebraically, not by running this checker at arbitrary sizes.
