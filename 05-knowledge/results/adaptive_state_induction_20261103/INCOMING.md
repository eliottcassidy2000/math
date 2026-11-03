# Incoming work and a new exact two-step valuation law

The twelve commits from base ff94a6d741 through
`44d94b4e50` (use the full revision in `sources.json`) were integrated
before publication. The incoming claims below are **CITED repository
results**, reviewed for their hypotheses and relation to this packet; their
large computations were not rerun. Section2 is a new elementary **PROVED**
consequence, independently audited by history, with a small exact receipt.

## 1. What the concurrent work changes

[Ternary repair](../ternary_repair_clock_20261103/README.md) supplies actual
95->287,289 and97->293,295 sources using two additional seven-vertex charts,
and all30 branches through four offset levels from H31 using the complete
seven-vertex bank. These are source replacements with integer carries,
not exact affine closure in the old three-coordinate state. Its two H293
sources have the same residue data but next q2 counts3557 and3865. The
closed carry module is7{(2u,x,y)}. This strengthens the actual-source and
joint-state requirement in our induction; it does not supply all-depth
branching. The283/287 native excursions in our packet are different literal
operations, and we do not claim priority for realizing the numbers themselves.

[Session154 closure](../session154_gadgets_20261103/closure/NOTE.md) bounds
the scalar support of a fixed branching family when uniform growth outpaces
branching. That obstruction has different hypotheses from our finite-bank
valuation obstruction. Neither applies to arbitrary adaptive operations.
[Mixed prefixes](../session154_gadgets_20261103/full_prefix/NOTE.md) fill the
primitive-pair orbit modulo7*13^2 in both parity layers, removing the earlier
pure-q joint phase obstruction at that modulus. The concrete recovered
7-divisible count has composite quotient, illustrating that modular coverage
does not give a prime/7p realization theorem.

The [13-adic image](../session154_gadgets_20261103/padic13/REFINEMENT.md),
[general lift test](../session154_gadgets_20261103/padic_general/GENERAL.md),
[7/13 joint image](../session154_gadgets_20261103/crt713/JOINT.md) and
[7/19 synchronization](../session154_gadgets_20261103/crt719/CRT.md) retain
singular phases and compatible exponent residues. They reinforce the scope
of our cap-two clock: an invertible register is not automatically a
surjective scalar observation, and separate residue coverage is not joint
CRT coverage. No conflict with the exact observer or return theorem was found.

The [uniform-triangle identity](../session154_gadgets_20261103/interfaces/README.md)
is an explicit specialization of homogeneous module substitution. Its
H'=3H+6P2+6P3 cannot supply nonzero ternary offsets. The displayed unequal
ports break that symmetry, agreeing with the order-five need for named
incidence. This is a connection of mechanisms, not an equivalence of the
different source languages.

## 2. A two-step s78 clock: exact3-adic progress and a remaining49-index gate

The incoming [CARRIES C4](../ternary_repair_clock_20261103/CARRIES.md)
recovers a symbol s78 preserving the **same** error line
tau=(2,-1,0) as b and q2. In our saturated bulk coordinates
beta=(Y,X2)=(y,h+2x-y), its row-action matrix is

    Q = [[77,152],[28,55]],   det Q=-21, trace Q=132.

Thus it is outside RETURN_BANKS'3-unit hypothesis. This is a same-channel
escape, unlike our q5 example, which changes the error channel. The exact
effect at every depth can be determined without a search.

**Proposition I1 (all nonzero integer bulk rows).** Let nu(beta) be the
minimum3-adic valuation of its coordinates. For every k>=0,

    nu(beta Q^(2k))   = nu(beta)+k,
    nu(beta Q^(2k+1)) = nu(beta Q)+k.

The one-step increment nu(beta Q)-nu(beta) is0 or1. Hence its increments
alternate0,1 or1,0. For the actual H31 seed beta=(25,34), the law is
nu(beta Q^j)=ceil(j/2).

**Proof.** Cayley-Hamilton gives

    Q^2 = 132Q+21I = 3U,
    U = 44Q+7I = [[3395,6688],[1232,2427]],  det U=49.

U is invertible over Z_3 and therefore preserves minimum coordinate
valuation. Since Q commutes with U, both displayed identities follow.
An integral Q cannot decrease the valuation; an increase of at least2
would contradict the exact two-step increase1. For the seed,
beta Q=(2877,5670) has valuation1, proving its formula. Nonzero rows stay
nonzero since det Q is nonzero. This is an all-depth proof, not an inference
from the finite powers in the receipt.

The removed3-adic obstruction exposes an integer lattice condition rather
than providing a return. The row lattice Z^2 U has index49 and exactly the
description

    (a,b) in Z^2 U  iff  a-21b = 0 mod49.             (I2)

Both rows satisfy the congruence. Its kernel has index49, equal to |det U|,
so inclusion is equality. Equivalently U has Smith invariants1,49.

Consequently a proposed full-state return

    w P_s78^2 = 3v + m tau

for an integer source state w, with any integer m, requires
beta(w)U=beta(v), hence (I2) for beta(v). The H31 seed bulk(25,34) has
residue46; the H95 and H97 seed bulk(75,102) has residue40. All three
therefore fail this proposed exact two-symbol return **before** actual
source realizability is considered. Longer mixed returns, carries and
other seeds remain outside this exclusion.

The connection map is correction-preserving source words -> saturated
bulk ->3-adic level plus integral row lattice. It preserves exact transverse
transport, forgets the correction coordinate and graph identity, and needs
the source word, carry and named-interface data for reconstruction. The
cheapest next test is the mod49 row-lattice gate, followed by the actual
source inverse. Merely observing the alternating3-adic increments is not
enough to close an induction.

`incoming_clock.py/json` checks the exact matrices, determinants, row-lattice
congruence, seed residues and twelve powers as controls, normally and under
optimization. It introduces no large support scan. The live concept board
now has a more precise return question: find a mixed word and actual source
family that jointly solve the3-level,49-index and full-state carry conditions.

## 3. Final incoming legality guard

Commit `7516d6f461c0c65adfc84cb11aff1854e3ed4867` adds
[return channels](../return_channels_20261103/README.md), reviewed after the
first checkpoint. Its [port-guard theorem](../return_channels_20261103/PORT_GUARD.md)
proves `S(T star B,q)=kappa(B;b,q)S(T,a)` for the full port-start count.
The native bank has positive kappa, so one extra bit exactly preserves
future start feasibility, alongside its three H-observer coordinates.
This is a different observer and language from our seventeen-query theorem.

The [context-entry hostile](../return_channels_20261103/CONTEXT_ENTRY.md)
makes the source-side distinction literal: a five-vertex source and a
four-vertex bank representative have the same triple(3,3,3) but port-start
counts0 and1. A context entering the observer image cannot justify replacing
the supplied graph or cancelling the context backward. For our proof plan,
native selected-path resets therefore need their own legality witnesses;
equal seventeen-query states certify future H under our declared operations,
not arbitrary native path availability or the vertex clock.

Its [three-birth exclusion](../return_channels_20261103/BIRTHS.md) strengthens
the old H55/111 hostile for the six specified bank presentations; it remains
a bounded induced-extension result. The common negative-offset family and
error-parity barrier in that packet do not establish an all-size return.
These are cited incoming results, not additional computations of this session.
