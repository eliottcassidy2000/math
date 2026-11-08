# Session225 joint modular law: reconstructed recovery

**PROVED:** exact fixed-length two-J distribution at every fixed odd modulus, with O_q(1/L) error. **FINITE-EXACT:** fresh controls explicitly listed in VERIFICATION.md. This packet reconstructs retained Session225 reasoning after the original temporary artifacts were unavailable. It does not claim those artifacts survived, does not overwrite them, and makes no priority claim.

- [PROOF.md](PROOF.md): full state maps, quantifiers, finite-sum law, SCC line-transversal proof, and four-branch CRT mixture.
- [check.py](check.py): portable exact checker, Python standard library only.
- [VERIFICATION.md](VERIFICATION.md): fresh run scope and reproduction commands.
- [AUDIT.md](AUDIT.md): fresh independent review, added at the audited checkpoint.
- [KERNELS.md](KERNELS.md): independently checked prime-power branch formulas, added at the audited checkpoint.

The key distinction is that broad gcd classes can split into several SCCs. Uniform-line inputs nevertheless meet each SCC proportionally, and scalar output is uniform on the invariant coset. CRT preserves a common branch choice, producing a mixture of products rather than a product of averaged local factors.

Example: L_15(0)=5179/3375 while L_3(0)L_5(0)=5134/3375. Their difference is 1/75.

## Scope firewall

This is a law for uniformly sampled words of fixed length with exactly two Js. It is not a height-conditioned law, an independence assertion for real history, a finite-bank completeness result, a scalar coverage theorem, or a resolution of the Chan-Pak cutoff. No guardrail, maintained warning, or existing result is changed. Session225 navigation adds a separate, scoped link to this recovery.

## Provenance and dependency boundary

The reconstruction uses retained model definitions and the Session223 all-3-adic / THM-4916 odd-prime cascade as motivation. The current proof is self-contained. Original Session225 finite outputs are historical context only; fresh runs are distinguished in VERIFICATION.md. No theorem identifier is reserved or reused.

The packet is published additively on a fresh main base, preserving the separate graft recovery at 6872236362fbb7077314e25f3a66a7d0070c5402 and all concurrent work. Each checkpoint receives exact immutable file verification before publication.

## Verification limits

The mathematical checks listed in VERIFICATION.md were run. Repository-wide checks and full CI were not run; the focused checks are not a substitute for them.
