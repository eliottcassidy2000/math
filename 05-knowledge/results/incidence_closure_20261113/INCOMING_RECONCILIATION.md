# Incoming exact terminal geometry and the remaining growth interface

**Status: inherited PROVED / FINITE-EXACT in the cited scopes; no new global
spectrum or count-recursion theorem.** Baseline a8b6d8dbf0c6b601813785efb0b5904b46c23d43
includes eight Session317 commits after the previous packet. All four replay
programs below were rerun successfully during this session. Their source bytes
are bound in SOURCE_PINS.json.

The [terminal atlas](../session317_exact_anchor_atlas_20261113/TERMINAL_ATLAS.md)
keeps a physical terminal profile (h,E,s,c), not merely its scalar Hamilton
count. Its binary-prefix law is

    H(B)=c+(h+E)B-h K(B)-s(B mod 2),
    K(0)=0, K(2B)=K(B), K(2B+1)=K(B)+B+(B mod 2).

The three strong terminal profiles (9,245,0,91), (9,245,0,93),
(9,245,0,95) give actual consecutive odd triples for every positive B.
The mixed second difference is -2h: nonzero source mass prevents affine
cancellation of the binary curvature. This is a concrete discrete derivative,
separate from both the reversible observer Gram matrix and stress derivatives.
Replay_terminal.py reproduced these profiles and its two (3,12,3,13/15)
controls.

The [one-return-star factorization](../session317_exact_anchor_atlas_20261113/ONE_RETURN_STAR.md)
retains the coupled boundary response

    H = alpha beta K(U) + (h_I d - h_W lambda) H(U),

where this K(U) counts ordered two-slot covers with empty slots permitted;
it is NOT the binary K(B) above. The negative exposed-boundary correction
lambda is essential. The inverse-design family

    H=(2^(k+1)+1)(2^k+1)2^m+1, k,m>=1

uses exact cancellation while retaining actual graph realization. Its replay
passed 7,623 literal factorization controls, 247 transitive-mask controls and
36,864 radix-seam controls. This is the relevant response-level counterpart
of our signed difference D: each records information erased by taking a sum,
although the carriers and operations are different.

The [exact inverse](../session317_exact_anchor_atlas_20261113/EXACT_INVERSE.md)
splits B into high and low bits, retains reversed low bits and one parity seam,
and recovers all prefixes for a fixed physical terminal in
O_u(sqrt(N+1) log^2(N+2)) integer operations. Its replay agreed with exhaustive
prefix values. It is target-aware within that terminal language. The missing
step for an unrestricted spectrum is selection or renewal of actual terminal
profiles, including source weights and shared ports.

The [radix/capacity distinction](../session317_exact_anchor_atlas_20261113/RADIX_AND_CAPACITY.md)
gives a complete arithmetic mixed-power decoder for coprime q>h>=2 and
large enough digit interval, while a fixed physical bank satisfying
sum_t M_t/q_t<1 has only O(X^alpha) outputs, alpha<1. Its replay checked
372,000 literal integer decodes. Arithmetic coverage is not a physical
spectrum theorem; the capacity obstruction explains why simply adding more
fixed affine charts can remain insufficient.

The new observer census complements these results. Within one order, retain
an invariant observation algebra or its full state partition. Across orders,
retain a constructor with a valid joint terminal response. The signed matrix
D restores the ordinary incidence erased by converse merging; the terminal
profile restores growth response erased by H. Neither restoration by itself
supplies the other's realization theorem.

The existing scalar A000568 recurrence is

    a_0=1, n a_n=sum_(k=1..n) g_k a_(n-k),

with independent connected gcd-cycle inventory coefficients g_k; see the
[prior recurrence summary](../joint_metagraph_dyadic_20261113/README.md).
Our fixed-order rational walk resolvent does not reduce the cost of computing
g_k. The concrete next interface is an inter-order constructor carrying
(A,D), source masses and enough shared-port response to deduplicate outputs.
A finite matrix identity without that constructor does not count new classes.
