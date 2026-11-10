---
id: THM-5097
title: "Factorial-gauge integral cycle-response recurrence and safe composite-modulus A000568 evaluation"
status: "PROVED elementary change of basis of the inherited response quotient, independently audited. Timings are FINITE-EXACT implementation comparisons, not a universal speed theorem."
session: Session239, 2026-11-10
source: 05-knowledge/results/session239_sequence_engines_20261110/burnside/PROOF.md
audit: 05-knowledge/results/session239_sequence_engines_20261110/audit/README.md
---

# THM-5097: integer states for the exact cycle-response recurrence

[Full proof and implementation](../../05-knowledge/results/session239_sequence_engines_20261110/burnside/PROOF.md)
· [Independent audit](../../05-knowledge/results/session239_sequence_engines_20261110/audit/README.md).

At fixed horizon N, retain processed cycle mass d and all future gcd responses
R_j. Scale the inherited rational state by d!. Adding k cycles of a previously
unprocessed odd length l has the integral multiplier

    C(d,l,k) 2^Delta,
    C(d,l,k)=(d+lk)!/(d! l^k k!),
    Delta=k(l-1)/2+l k(k-1)/2+kR_l.

The responses update by R'_j=R_j+k gcd(l,j), and coordinates j>N-d-lk may be
discarded. Terminal degree n is n! A000568(n). This is the existing exact
future-response quotient in an integral divided-power basis, not a new
Burnside formula or a polynomial state-count bound.

Running in Z/(M N!)Z permits exact terminal division by n! and recovers every
prefix term modulo any positive M. For M=2^b, b+v2(N!) bits suffice: remove
2^v2(n!), then invert the odd factorial part. Zero states and shifts beyond
the working precision may be pruned in residue mode only.

The saved same-machine prefix comparison through200 is7.48s versus27.14s
for the inherited Fraction quotient, with exact equality. The16-bit prefix
uses1.12s. Historical compiled code wins at some smaller orders. These are
reproducible observations, not an external fastest-known claim.

Dependencies: Davis/Burnside as in [A000568](https://oeis.org/A000568),
THM-514 (`the-a000568-unit-tail-core-kernel-and-prime-tail-ladder`) and the
response quotient in `05-knowledge/results/tournament_future_boundary_20261020/README.md`.
The same integral gauge extends to the inherited reversal-palette family:
[palette proof](../../05-knowledge/results/session239_sequence_engines_20261110/burnside/PALETTE_PROOF.md) and
[independent audit](../../05-knowledge/results/session239_sequence_engines_20261110/audit/palette_and_asymmetry_audit.md).
It supplies A000088, A000273 and A001174 as well, without claiming a measured
speed ratio for those sequences. Scalar class counts do not determine
Hamiltonian-path-value support.
