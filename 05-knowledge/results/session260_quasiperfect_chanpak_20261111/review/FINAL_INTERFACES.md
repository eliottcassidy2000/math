# Focused final-interface review

**Verdict: PASS, no material correction required.** Reviewed new arguments in `bridges/DEFECT_QUANTITATIVE.md` and `coverage/INCOMING_CONNECTIONS.md`. This is a mathematical/interface review, not a replay of the large census or an inherited-literature audit.

## 1. Positive original-arc cut law

The partition by whether deleting z leaves a Hamilton path is exhaustive and disjoint. For a fixed old path, the binary-word identity counts all legal insertion gaps, including both endpoints, as one plus the number of A-to-B adjacencies. For a failed deletion, the consecutive old vertices necessarily have the reverse old arc; reversing that single arc gives exactly the claimed repaired path. The marked arc is recoverable, so no multiplicity is silently lost. This proves the sum of L_uv+R_uv over forward old cut arcs directly.

The empty-old-tournament convention H(empty)=1 gives zero increment and is separately handled. Arc-flip parity proves even weights, and the inherited raw kernel identity supplies their strict positivity. None of this presupposes the conjectural stronger domination bound.

## 2. Coefficient comparison and circulation

The two cut laws agree on every Boolean subset. Both are multilinear polynomials of degree at most two, so uniqueness of Boolean multilinear coefficients justifies comparison. The squarefree coefficient at x_u x_v receives one original-arc term and both raw directed capacities, yielding K_uv=2(c_uv+c_vu). At x_u the original-arc outgoing sum then cancels against the raw outgoing sum exactly as written, leaving

    sum_(u→v) c_vu = sum_(v→u) c_uv.

With r(u,v)=c_vu this is precisely outgoing-equals-incoming flow, not the opposite orientation. It licenses moving total reverse capacity across a cut, not discarding a nonarc capacity individually. The transitive specialization has zero reverse flow and reproduces the stated distance weights. The distinction between raw endpoint statistics and the resummed arc weights is clear.

## 3. Selected fibres and spectator transport

The all-selected dominating set and the sets omitting one nonisolated vertex are distinct. Their cuts are respectively(B,A), (B−b,A+b), and(B+a,A−a), yielding the endpoint products in equation(7). Every selected set is dominating under precisely the indicated nonisolation guard. The inherited nonnegative fibre identity therefore permits discarding the other fibres. This bound is not implicitly substituted for H(U)−1.

For an ordinal spectator W before T, every spectator vertex lies on z's incoming shore, while each old outgoing-shore vertex loses to it. Thus it is isolated in the new bipartite graph. The dual spectator-after operation has the same property. Ordinal multiplication gives both H identities, and direct algebra yields the stated defect, increment, and conjectural-slack transports. It follows only that the conjectural bound is preserved from a base where it already holds.

**Required scope preserved:** D≥H(U)−1 and D≥2t_3(U) remain OPEN. The note explicitly identifies the missing global injection/Hall obligation and treats exhaustive/random survival only as evidence. There is no promotion of either inequality in the reviewed text.

## 4. Two-root coarsening and gap birth

The outgoing masks of B verify the two exterior-twin pairs2→1 and4→3. Contracting them in the stated exterior orders gives exactly Q_I=(6,8,10,1) and Q_II=(8,1,3,6). The u≥1 and u=0 descriptions map to the previously supplied marked families with their stated parameter shifts; the zero-power boundary term is retained.

For g=0, q+delta=5C·2^u3^v, and one intervening clone changes the factor5 to13. Hence q'=(13q+8delta)/5 is integral and is the actual marked gap insertion. Direct substitution verifies the L/G, R/G, and L/R response commutations. The next insertion gives D=3−12·2^(u+v)<0, so G is correctly a one-use guarded operation, not a repeatable digit.

The parity distinction is exact: u=0 retains positive odd delta and odd q, whereas u≥1 has delta0 and even q. Therefore same-q rearming is impossible inside this carrier. At u=1 the inverse L boundary subtracts e; the examples60→29 and100→47 correctly expose why ordinary halving would be wrong. The eight slope sequences follow by stripping powers of two and retaining the four odd boundary rows. This is a complete statement inside the named carrier, not an all-tournament exclusion.

## 5. Conditional head replacement

If every even multiplier q≥2B exists, its one-tail outputs cover1+4p for p≥B. Separate direct constructors for1+2p for every odd p≥B cover the other residue class beyond the claimed sufficient threshold4B+1. The argument requires actual direct Hamilton constructors but not seed deflation, exactly as stated. Incoming R10 is used only as an example on its specified residual family; the missing universal hypothesis is not inferred from it.

No new large computation was needed for this pass. The verdict covers these proof interfaces and their present quantifiers, not future changes or a Chan–Pak completion.
