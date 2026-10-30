# Connections with exact maps and explicit losses

**PROVED** here: the elementary residue-bank capacities and the stated
interfaces to the new results. **INHERITED PROVED:** the specified lattice
and aliquot results. **CITED CLAIM, not independently validated here:** the
factoring paper's general complexity theorem. Numerical resemblance alone
is never a dependency.

## 1. What the supplied paper contributes

Rudich–Rousseau,
[*Deterministic Factoring in Time N^{o(1)}*](https://hexagonmath.org/pdf/2610.00600v1),
was read as a six-page PDF, including the displayed formal statement.
Its useful interfaces are prefix-count search, uniquely extended local
constraints, and modular zero tests with an integer-height bound. The claimed
counting engine and its Lean complexity proof were not independently built
or audited here. The stronger explicit runtime in its conjecture is not a
proved dependency. Incoming Session112 separately inspected the formal
statement and solution and records that scope in
[PAPER_SCOPE.md](../session112_counting_interfaces_20261030/paper/PAPER_SCOPE.md).

For this session the transfer is elementary: a residue-zero observation is
useful only with a known integer range, just as a modular Collatz identity
needs a clock/height bill and an overlapping construction count needs an
exact union. None of our results uses the paper's claimed exponent saving.

## 2. The requested number triples as a deficit certificate

For a labelled tournament union U on n vertices let
D=2^binom(n,2)-|U|. It is an integer between0 and2^binom(n,2). If all its
residues modulo a bank m_i vanish and

    lcm(m_i)>2^binom(n,2),

then D=0. This is an unscaled special case of the stronger incoming
[modular-observer theorem](../session112_counting_interfaces_20261030/arithmetic/MODULAR_OBSERVERS.md),
whose correct modulus for a scaled observation sD is L/gcd(L,s).

| chosen bank | exact lcm L | largest entire labelled tournament cube below L |
|---|---:|---:|
| 3,5,7 | 105 | n=4: 64 points |
| 11,19,27 | 5643 | n=5: 1024 points |
| 1001,343,57 | 2795793 | n=7: 2097152 points |

Each next tournament order exceeds that bank's capacity. This is a designed
height comparison, not a unique arithmetic characterization of the triples.
The last bank is **not coprime**: its product is7L. The nonzero integer D=L
has every residue zero and is smaller than the product, so product-based
capacity would falsely certify zero. This is a concrete reason to retain
valuations rather than just prime support. The actual first-layer five-vertex
union has deficit1024-280=744, detected by the banks.

Source: an exact union deficit. Target: a zero/nonzero certificate. Map:
residue evaluation. Preserved predicate: zero when the range and lcm guard
hold. Lost information: integers differing by L, or by L/gcd(L,s) after
scaling. Needed sidecar: scale and height. Cheap decisive test: D=L.
Computing the residues is not asserted to be cheap. Burnside division cannot
be performed by inverting a nonunit group order without retaining the scale.

## 3. The stronger existing map 11,19,27 -> 3,5,7

The recovered
[prime-power-plus-one arithmetic packet, L1–L2](../prime_power_plus_one_20261029/ARITHMETIC.md)
already supplies a structural explanation. For any doubly regular tournament
of order q, on its integral sum-zero root lattice L let theta=(I+S)/2.
With C=(q+1)/4 and r=(q-1)/2, it proves

    theta(1-theta)=C I,
    L/theta L=(Z/C)^r,

and a perfect pairing with the complementary cokernel. At q=11,19,27 the
pairs (C,r) are(3,5),(5,9),(7,13). This is a map between actual operators and
finite groups, stronger than noticing two arithmetic progressions.

The new intersection theorem tests a **different invariant**: whether two
transported regular cyclic subgroups generate an edge-transitive group.
Prime q=11 and19 have the point-intersection property after a pin; q=27
permits a common16-point chart in the explicit C3[circulant9] family.
The doubly regular Paley27 object has regular translation group C3^3 and is
not that C_27 cyclic family. The lattice map is uniform across these q;
the cyclic-intersection proof is not. Keeping the regular group prevents an
incorrect transfer across the shared integer27.

## 4. 1001 -> 343 -> 57: recover the older mechanism first

The elementary aliquot calculation is

    sigma(1001)=8*12*14=1344, s(1001)=343,
    sigma(343)=1+7+49+343=400, s(343)=57,
    343+57=20^2, 57=3*19.

The stronger characterization was already proved in
[Session82](../session82_decimal_aliquot_square_sum_20261028/ANCHOR_CHARACTERIZATION.md)
and inherited by
[the divisor-overlap packet](../aliquot_divisor_overlap_20261028/THEORY.md).
We do not replace it by rediscovering the edge. For a squarefree triple pqr
with aliquot sum343, the exact search compression is

    (p+q+1)(p+r+1)=343+p(p+1).

The p=7 factorization399=19*21 yields q=11,r=13. This is an arithmetic
factor-pair certificate, not a tournament intersection computation.
The newer
[Session111 preimage classification](../session111_scale_connections_20261030/ALIQUOT_MERGE.md)
has the complete19-source preimage list of343, including the even578=2*17^2;
the decimal/square-sum subfamily is not the whole fibre.

Incoming Session112 goes further in a different direction:
[paired aliquot residues](../session112_counting_interfaces_20261030/arithmetic/ALIQUOT_MODULAR.md)
recover n modulo p^k and produce coordinates(r,C) in which every prime-ray
step is(r,C)->(p^t r,C). This is a proved compatible inverse-limit
coordinate system with an invariant arithmetic core. The source residue can
become zero while that core persists. It is useful inspiration for preserving
Collatz's two affine coordinates through refinement, but the transition is
not conjugate to the Collatz map: multiplication by p contracts here, while
the Collatz coefficient updates are units at p=5.

The older
[57/27/16 atlas](../number_bridges_57_27_16_20261028/ATLAS.md)
also identifies57 as an Eisenstein norm, with
1+8omega=(1+2omega)(5+2omega)=omega conjugate(8+omega),
and explains which norm-character information is lost. This makes the
valuation-plus-phase requirement precise in another setting. The norm is
not a source-dependent Collatz payment bound.

## 5. A common research move, with distinct mechanisms

| lane | comparison object | restrictive invariant | surviving information |
|---|---|---|---|
| prime tournament charts | difference of two orientations | invariant under both cyclic parents | constants only; one pin kills them |
| fixed Collatz controller | two orders of affine composition | commutator numerator modulo content | full carry, entry word and phase |
| aliquot observer | two readings along a prime ray | compatible residue isometry | source residue plus arithmetic core |
| modular count | observations of one integer | lcm and scale kernel | a unique integer only below known height |

The proposed broad principle is **an approach, not a universal theorem**:
compare two constructions before forgetting their shared frame, compute the
smallest invariant that their agreement forces, and then test whether a
marked pin or a height bound eliminates the ambiguity. Common symmetries
can leave a large family (the3m hostile); large valuations can be inactive
(the dark-defect word); a zero residue can hide a nonzero multiple of L.
These are different ways the missing coordinate becomes visible.

The least-used sidecar that paid off here was the older faithful affine
monoid/carry cocycle, THM-4604, rather than another trace invariant. It yields
the marked-entry divisor in COLLATZ.md(10). The corresponding tournament
improvement came from applying an old group theorem to orientation
**differences**, not by assuming complements preserve automorphism groups.

## 6. Specific next targets

1. Prove or refute the phase-aware discounted epoch budget, preserving its
   initial word and unused low-height clock. A refutation must exhibit the
   exact word and distinguish budget failure from content-bound failure.
2. Construct corrections in the observed bounded K interval at every layer,
   using old edge slack. A local digit rule must transport the pair(m,b)
   together, not independently select two favourable lifts.
3. Bound the entry-dependent commutator valuation in changing-controller
   compositions, treating exact word resonance separately. A fixed-period
   bound is already supplied; the target is uniform entry control.
4. Replace composite cyclic parents by a recorded block quotient and internal
   difference-orbit modules. The3m family predicts precisely what the shared
   ports fail to pin. Test whether intersection terms can be indexed by those
   modules under repeated native moves, before claiming any unrestricted
   A000568 complexity improvement.
