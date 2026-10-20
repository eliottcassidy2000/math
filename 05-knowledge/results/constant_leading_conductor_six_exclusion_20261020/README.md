# Constant-leading reduction excludes the conductor-six sector

20 October 2026. **PROVED; independently reviewed PASS.**

For P,Q∈C[x,y] with one coordinate monic in y and J(P,Q)∈C*, the generic source algebra C(x)[P,Q] cannot have degree semigroup <3,4>.

The new mechanism forces constant leading coefficients after every target degree reduction. The only terminal degrees are (3,4) and (4,6). In the exceptional class, the quartic coefficients force the generic canonical cubic to be polynomial, and immersion makes its recovery coefficient a unit. Thus the whole source algebra has monic cubic/quartic generators, contradicting the earlier moving obstruction. Independent direct Jacobian coefficient arguments exclude both terminal classes without the moving-period theorem.

Combined with the [prior low-delta theorem](../low_delta_source_conductor_exclusion_20261020/LOW_DELTA_QUADRATIC_AND_RELATIVE_GLOBALIZATION.md), a noninvertible source under its monic/primitive hypotheses has generic source-conductor degree at least eight. This is not a bound on the map's function-field degree and does not resolve the planar Jacobian conjecture.

- [Proof and formal-jet boundary](CONSTANT_LEADING_GLOBALIZATION.md).
- [Independent audit and direct exclusion](review/AUDIT.md).
- [Dependencies and primary-source checks](DEPENDENCIES.md).
- `verify_constant_leading.py`: 130 exact producer controls.
- `review/independent_checks.py`: 77 independent controls.

Both checkers pass normally and under Python -O with identical outputs. The formal-jet hostile shows every immersed polynomial curve has an all-orders formal Keller lift, but the nodal example's y-degrees grow without bound. It does not preserve the generic semigroup.
