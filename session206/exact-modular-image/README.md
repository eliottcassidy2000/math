# Session 206: exact modular images of fixed O/Z tail sources

This checkpoint proves the full projection of a fixed physical source's O/Z scalar image onto every finite modulus. For even M and `g=gcd(R,d,M/2)`, the image consists exactly of the two classes `h` and `h+R−e` modulo `2g`; each permitted residue has arbitrarily large representatives. The root source consequently has every parity-compatible residue modulo every modulus.

The construction uses separated blocks of consecutive binary ones in the exact session205 formula `H(m)=h+Rm−dQ(m)−e(m mod2)`. It is an unconditional all-order proof, not an extrapolation from finite scans.

**Scope:** this rules out obstructions based solely on a missing fixed congruence class for a bank containing the root. It does not rule out growing-modulus arguments combined with the inverse bound, and it proves neither eventual coverage nor a general finite-bank impossibility.

## Contents

- `DYADIC_SURJECTIVITY.md`: self-contained proofs, culminating in the exact modular-image characterization.
- `verify.py`, `RESULT.json`: 54,102 dyadic witness checks on 142 physical sources.
- `verify_all_moduli.py`, `ALL_MODULI_RESULT.json`: 20,556 checks of the coprime-d version, all even moduli through 64.
- `verify_exact_moduli.py`, `EXACT_MODULI_RESULT.json`: 62,832 explicit low-prefix/block-pumping checks of the exact characterization, 21 physical sources and even moduli through 32.
- `SHA256SUMS`: frozen-file integrity manifest.

Run each Python file with Python 3. They use only the standard library. The scripts are evidence checks supplementary to the proofs, not proofs of the universal quantifiers.

## Source

[Session205 exact arithmetic reduction](https://github.com/eliottcassidy2000/math/blob/da58c1df8c424655f17b28a7b9ddfa6fde11737c/session205/finite-oz-tail-bank/ARITHMETIC_REDUCTION.md).
