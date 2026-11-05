# Independent review: universal central-factorial context positivity

**PASS**, no corrections requested, 2026-11-05.

Frozen source in `session197-odd-cycle-path-polynomials/`:

- UNIVERSAL_CONTEXTS.md: 8615845e7bf06a6cdb1925fdbd8198a77b1a8864aa9c22704b0459678ff2ed51
- verify_extension.py: 3b537a9cbfea66aabe8b1d4cf59c84de45cc41d2ea6f80d053df8c083b1cad76
- Normal/optimized output: 4bbb4df5462fbd95805fd4647a2514c291cd318172ddf0a0d46c5c1802b4f9fc

Both copied-source replays reproduce the output hash. Dependencies are the previously reviewed companion proof and verifier; no source was edited.

## Full operator and multiplicities

The maximal-module-run argument defines a linear map on the complete falling-factorial profile, not merely on H. Its recombination polynomial depends only on the number of input fragments because every outside incidence with the intact module is uniform. Fragments are distinguishable by their disjoint vertex labels, so the recombination count includes their ordering multiplicity. The monomials from transitive inputs, including x^0=1, span Q[x], establishing uniqueness and permitting the later extension to virtual basis vectors.

Every crossing odd cycle has exterior vertices and at least one marked clone gap. Canonical rotation at its least exterior vertex, followed by ordering disjoint objects by their anchors, avoids both rotational and collection-order overcounts. An assignment of named clones to marked nonempty gaps fixes every run's increasing order. The forward construction and contraction are inverse.

For one crossing skeleton, (S-I)^t counts assignments with all t run bins occupied, while (S^(-1)-I)^t weights them by the parity of the total number of assigned clones. The sign (-1)^(1-l) is therefore correct. For several cycles, products of these commuting shifts encode assignment to disjoint named run-bin groups. Each OCF factor2 cancels exactly one projector denominator2. The remaining exponent is unused exterior vertices plus exterior-only cycles plus crossing cycles. Thus X must act after the shifts, exactly as written.

At m=0, every crossing term vanishes because t>=1 and its difference operator kills constants. The noncrossing terms give precisely the deleted-root graph's profile. Negative intermediate shift arguments are signed polynomial evaluations, not negative physical bin sets. No forbidden exchange of X and a shift is used.

## Positivity and parity

The lowering identities include the essential separate constant basis element. In particular delta q_0=2 and Delta r_0=2; these cannot be discarded or assigned to a nonexistent q_(-1). All four general formulas and all boundaries are correct.

For y=S-I and z=S^(-1)-I, yz=-Delta gives the positive recurrences for y^t+z^t and (y^(j+1)-z^(j+1))/(y-z). The latter is explicitly a polynomial geometric sum, not an attempted inversion of a noninvertible operator. Every skeleton term is consequently a composition of positive X, delta, and Delta operators, with a nonnegative integer coefficient. This proves positivity on every basis vector, not just on physical inputs. The displayed parity sum equals the number of exterior vertices, including the exterior-only odd cycles.

The coefficientwise lower bound follows by retaining only the physical input's u_0 coordinate. It does not require arbitrary cone vectors to be realizable. No simultaneous multivariable positivity is asserted.

## F/G physical escapes

The pinned S193 multi-port theorem was checked in the earlier audit. For the triangle/two-port core, the two mixed responses are x+1 and x-1, and the start/end sums are twice x^2-x+1 and x^2+x+1. These give F and its manifestly positive form exactly.

For the reverse transitive four-core, the first-bin choices give the displayed L and R, while the two ordered mixed-port cases give x(x+1) and x(x-1). Thus G and its positive form also follow without fitting. The readouts 13h+36u_1 / 25h+72u_1 and 19h+60u_1 / 39h+120u_1 have the correct input parity. The explicit H13 and H19 rows have the stated orientation and are strong. These operators genuinely lie outside the earlier D/E grammar, so its omission statements remain consistent.

## Independent finite evidence

The new `path-polynomial-audit/universal_independent.py` does not import the source verifier or skeleton enumerator. It reconstructs full monomial columns from literal directed odd cycles of actual clone graphs. For contexts of order at most four, every cycle uses at least two exterior vertices, so no two cycles can be disjoint; the OCF profile is then reconstructed directly from that literal cycle list.

Across all285 rooted labelled contexts of orders1 through4 it checks clone sizes0 through8, yielding2,565 full physical columns. It applies those columns to signed q/r monomial expansions through index3, independently confirming2,280 virtual basis columns, all constant columns, and the required parity. This checks the positive operator claim on virtual inputs without assuming the skeleton positivity theorem in the expected calculation.

A separate generic odd-cycle packing computation verifies F/G full profiles on four different input graphs, including a directed triangle and a strong four-vertex graph, and reproduces the exact H13/H19 witnesses. The supplied normal/-O replay additionally verifies1,425 actual clone profiles,168 literal parity projectors,1,140 source basis controls,150 physical F/G profiles, and lowering identities through index30.

The all-order conclusion rests on the skeleton formula and positive operator recurrences. These finite controls independently test its conventions and boundaries; they are not an inference of universal positivity from a census.
