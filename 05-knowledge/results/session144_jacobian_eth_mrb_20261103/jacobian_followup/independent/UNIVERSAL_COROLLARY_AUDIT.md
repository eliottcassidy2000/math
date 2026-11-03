# Independent audit: universal one-positive-weight boundary tower

**Verdict: accepted as stated.** Audited artifact: [universal boundary corollary](../universal_equivariant_boundary_corollary.md) from the Session 144 producer. This is a proof audit, not a new formal build.

The six claims are valid for the same-action equivariant Keller category with weights `(1,-w_2,...,-w_n)`, all `w_i>=1`, over an algebraically closed field of characteristic zero.

1. The unique positive-weight block in the derivative at zero proves `l(0)!=0`. Negative-weight coordinates vanish on the positive-weight axis, so their invariant numerators vanish at the quotient origin and `H(0)=0`.
2. Given any quotient point on `l=0`, choose a source point in the chart `x=1`. In orbit coordinates, the differential of the first output is `d(xl)=x dl` there. A zero `dl` would give a zero Jacobian row, contradicting Keller. Hence the orbit-zero hypersurface is smooth.
3. Divisibility of all coordinates of `H^(i+1)` by `f_i`, together with `H(0)=0`, proves `f_j-l(0)` is divisible by `f_i` for `i<j`. Since `l(0)` is a unit, levels are pairwise comaximal.
4. At a point of `f_j=0`, all earlier levels are nonzero. The quotient Jacobian product is therefore nonzero, and pullback of the nonzero gradient of `l` proves smoothness of each level.
5. The determinant formula follows by the chain rule. Dominance of `H` follows from its nonzero Jacobian and makes pullback injective on polynomial rings. Thus a nonconstant `l` remains nonconstant under every iterate. Its zero sets are nonempty over the algebraically closed field. This yields at least one distinct irreducible component per level, not an irreducibility claim.
6. For equivariant polynomial automorphisms, composing the orbit scale with that of the inverse gives a polynomial unit. Thus the scale is constant. The inverse descends to the quotient, so equivariant two-sided changes preserve the number of orbit-zero components.

The scope exclusions are necessary and adequate: no universal exact degree recurrence, no universal description of the zero fiber beyond the proved collapse, and no invariance under arbitrary non-equivariant source/target automorphisms. Smoothness of the quotient orbit-zero hypersurfaces does not make the quotient map etale on them. The lifted Keller map remains etale everywhere.

The quartic boundary correction is compatible with this corollary: its invariant target hyperplane can lose sheets even though its source orbit-zero hypersurface is smooth. These are different predicates and require different data.
