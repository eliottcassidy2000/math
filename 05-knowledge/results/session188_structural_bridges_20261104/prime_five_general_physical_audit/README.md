# Reproduce the complementary prime-five physical audit

Run `python3 audit.py > normal.json` and `python3 -O audit.py > optimized.json`, then compare. The script reads the two frozen sibling packets `prime_five_all_profiles` and `prime_five_all_ambient` as pinned data and imports no producer code. The installed C++17 compiler is used temporarily; nothing is installed.

The review covers normalization, singleton and parity branches, complete ambient tails, physical naming/lift, and source-relative family quantifiers. It independently recovers all22 requested source physical controls, the four additional exact core cover profiles, all93 singleton-position controls, and the264-prime quotient orbit coverage. The separate bank auditor owns full kernel interpolation and exhaustive type classification.

Frozen allowlist: `REVIEW.md`, `README.md`, `audit.py`, `count.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`.
