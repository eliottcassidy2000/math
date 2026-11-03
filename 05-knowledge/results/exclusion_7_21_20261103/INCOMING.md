# Concurrent exclusion proofs and the limit of residue-only certificates

Status: **CITED all-modulus family / VERIFIED small sources and finite clocks /
PROVED conditional-on-cited-family consequence.** Incoming commit `63488fcae6`
contains Session157 work produced concurrently with this packet.

Its [gap certificate](../session157_gaps_20261103/gaps/gaps-note.md) imports
Busch's global lower bound, proves the needed Hamiltonian-cycle theorem,
and checks only512 fixed-cycle orientations at order6. Our alternative uses
an elementary strong-deletion proof and a larger fixed-path cover through7.
Both recertify the already-known exclusion; neither proves full coverage.
Their counts are for different presentation universes and must not be compared
as ordinary class multiplicities. The incoming
[odd-cycle packing reduction](../session157_gaps_20261103/cycles/README.md)
further explains the residual small cases, while retaining finite-check scope.

The incoming [prime-power/CRT family](../session157_gaps_20261103/prime_power/NOTE.md)
gives actual strong one-apex tournaments F_n, n>=1, with

    H(F_n)=h(n)=(10*4^n+6n-1)/9,
    |F_n|=2n+2,
    h(1)=5,  h(n+1)=4h(n)+1-2n.

Its increasing-prime CRT proof, together with ordinal sums of triangles,
realizes every residue modulo every fixed M, subject only to odd parity
when M is even. This is a **CITED theorem from that packet**; our checker
does not claim an independent all-modulus verification. It independently
constructs and counts F_1,...,F_4, giving5,19,73,287, and verifies that
n=1,...,14700 already covers all2450 odd residues modulo4900. In particular,

    h(1075)=21 (mod4900),  |F_1075|=2152,
    h(2634)= 7 (mod4900),  |F_2634|=5270.

These are residues of large positive counts, never exact21 or7. They match
the target's scalar H residues modulo49 and100 simultaneously. This does
not certify the preceding bulk-return lattice gates or a marked attachment.
The graph family and vertex clock are explicit; no small representative is
silently substituted.

**Consequence.** No sound exclusion criterion depending solely on
H modulo a fixed finite list of moduli can rule out exact7 or exact21.
Take their least common multiple M. The cited coverage theorem supplies
an actual tournament whose count has exactly the same residue data as the
forbidden target. Any rejection based on those data alone would reject
that actual tournament too.

This statement does not apply to order bounds, positive inequalities,
graph/source structure, marked responses, or an unbounded precision
procedure. Our exclusion succeeds by combining an all-order structural
reduction with a finite actual-graph classification. It does not succeed
by finding a single forbidden congruence. This distinction is useful for
the preceding mixed-return mod49/mod100 work: a failed gate may exclude
one operation language, but no fixed scalar congruence can be the global
reason for the7/21 gaps.
