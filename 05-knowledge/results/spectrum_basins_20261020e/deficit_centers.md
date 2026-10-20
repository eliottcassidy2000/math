# The center series beyond the single quartic

**PROVED, using the imported order bound and center classifications:**
every fixed-defect strong self-converse center series is a finite polynomial.
The deficit-two polynomial is exactly `z^4`, but the higher-defect series
need not be monomials: deficit six already occurs at orders five and six.
This is a refinement of the user's proposed power pattern, not a new
classification of every coefficient.

## Definitions and inherited mechanisms

Let `h=H/|Aut|`, let b be the number of blue Hamiltonian-path orbits, and
let `Delta=h-b`. Define

    D_d(z)=sum_[T] z^|T|,

where the sum is over ordinary isomorphism classes of **strong
self-converse** tournaments with `Delta=d`. The singleton is included;
the empty tournament is not. Thus the earlier D(z) is precisely D_2(z).
This is neither a raw-H stratum nor the feedback-arc deficiency series.

The closest mechanism is the cycle-cut rotation subgroup, now made
effective by a removable power-of-two cycle and two prefix attachments.
The [incoming cubic order bound](../effective_black_defect_order_bound_20261020/PROOF.txt)
proves, if `Delta<=d`,

    |T| <= 2(d+2) r_d < 16(d+2)^2(d+1),
    r_d = least power of two >= 4(d+2)(d+1).             (1)

The proof was read in full. Its rotation index q satisfies `q<=d+2`;
the q rotation orbits induce regular tournaments and force a large vertex
connectivity bound. A removable power-of-two cycle gives too many distinct
path orbits, at most one blue, if (1) fails. Moon's strong-tournament cycle
theorem is the CITED external input; the group and attachment steps are
proved in the imported packet. No large classification is needed here.

There are finitely many tournaments at each bounded order, so

    deg D_d <= 2(d+2) r_d                              (2)

whenever D_d is nonzero. This makes every individual coefficient problem
finite and effective; it does not make a census at this coarse bound small.
The effective bound strengthens the independent
[factorial-versus-exponential finiteness argument](../bounded_black_defect_finiteness_20261020/BOUNDED_BLACK_DEFECT_FINITENESS.md).
That earlier packet's statements that D_2 is still unknown are historical
scope boundaries, superseded by the separate classification below.

## Exact small centers and the failure of a monomial ladder

Reflection acts as an involution on the ordinary fixed-path fibre of a
self-converse tournament. Its fixed points are exactly the blue tilings;
its nonfixed points occur in pairs. Therefore Delta is even and

    D_(2j+1)=0.

The independently closed pure-blue alphabet gives

    D_0(z)=z+z^3+z^5.

See the [bounded proof composition](../rotation_owners_spectrum_20261020e/purity_closure.md)
and the separate uniform-prefix classification linked there. The audited
[deficit-two center theorem](../deficit_two_centers_20261020/session15-deficit-centers/DEFICIT_TWO_COMPLETE_CENSUS.md)
gives

    D_2(z)=z^4,       deg D_2=4.

The following two positive witnesses show that higher coefficients do not
continue as a sequence of single powers. A mask uses lexicographically
ordered pairs i<j, with bit one meaning i->j.

| order | mask | H | Aut order | h | b | Delta |
|---:|---:|---:|---:|---:|---:|---:|
|5|10|9|1|9|3|6|
|6|408|45|3|15|9|6|

Both tournaments are strong and self-converse. Hence, coefficientwise,

    D_6(z) >= z^5+z^6.                                 (3)

In particular D_6 is not a monomial. This does not assert equality in (3),
nor determine D_4 or any other full positive-defect center polynomial.
The coarse degree bounds from (2) are 1536 for D_4 and 4096 for D_6;
no enumeration at those orders is proposed.

The witnesses were recovered from the incoming independent bounded-defect
census. The separate [checker](deficit_centers_check.py) imports none of
its routines: it reconstructs the two matrices, checks reachability from
each vertex, and enumerates all permutations for automorphisms,
anti-automorphisms, Hamiltonian words, and blue words. The finite universe
is only these two explicitly specified positive controls. The parity and
degree statements are proofs, not extrapolations from these controls.

    python3 05-knowledge/results/spectrum_basins_20261020e/deficit_centers_check.py
    python3 -O 05-knowledge/results/spectrum_basins_20261020e/deficit_centers_check.py

The conceptual connection preserves the black-orbit budget and the full
center type. Passing to Delta alone loses vertex order and raw H, as the
two witnesses show. That lost coordinate prevents a single-power rule.
