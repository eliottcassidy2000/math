# Why a fixed finite bank of shared-module pairs cannot select all independent prime pairs

Suppose each of two Hamiltonian counts reads one common growing homogeneous module B through a fixed exterior, with possibly different exteriors for the two reads. A path in either quotient uses B in a bounded number of runs determined by the exterior. Write R_j(B) for ordered path covers by j **nonempty** paths. Contract each maximal B-run. The count in that exterior is a nonnegative linear combination of the R_j, with constant coefficients counting quotient words. Thus the first read is p=Σα_j R_j(B) and the second 7q=Σβ_j R_j(B), with no cancellations. Let r,s be their largest supported run ranks.

For n=|B| and j<n, splitting one of the j nonempty paths at a nonempty interior cut and putting the suffix in a new slot gives

    R_(j+1)(B) ≥ ((n−j)/j) R_j(B).

There are n−j choices for each ordered cover and at most j preimages of an output, by choosing the original slot and rejoining. For any fixed ranks, this shows that smaller ranks are uniformly negligible relative to larger ones as n grows. In particular

    p=α_r R_r(B)(1+O(1/n)),
    7q=β_s R_s(B)(1+O(1/n)).

If r≠s then q/p tends to zero or infinity. If r=s then q/p tends to the single rational number β_r/(7α_r), with error O(1/n). All these bounds are uniform over B of that order. Unbounded p implies unbounded order because each fixed exterior on bounded vertex sets has only finitely many realizations.

It follows that a fixed finite bank of pairs misses every independently prescribed prime pair in any sufficiently small positive ratio interval avoiding its finitely many equal-rank slopes, once p is large. Such intervals contain infinitely many independent prime pairs by the prime number theorem. This rules out a universal **simultaneous shared-module pairing** in that model. It does not rule out the requested unilateral conclusion p OR 7q: different sources for the two branches, variable exteriors, or a separate cross-axis implication can evade this argument.
