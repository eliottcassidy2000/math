# Curve inputs for the non-rational formal Hamiltonian flow theorem

**Status: CITED PRIMARY INPUTS; exact statements read September 6, 2026.**
These are the external dependencies of
[the fixed-time non-rationality theorem](../results/planar_jc_long_20260906_nonrational.md).
The hyperelliptic models, their squarefreeness, the computed genera, formal
convergence, and the source-coordinate comparison are proved locally.

- [Stacks Project, Theorem 53.2.6, tag 0BY1](https://stacks.math.columbia.edu/tag/0BY1):
  the equivalence between function fields of transcendence degree one and
  nonsingular proper curves converts a field embedding over `K(c)` into a
  nonconstant selfmap of its smooth proper generic curve. The direction
  is contravariant. Rationality of the inverse is not assumed.
- [Stacks Project, Section 53.12, tag 0C1B](https://stacks.math.columbia.edu/tag/0C1B):
  Riemann--Hurwitz applies in characteristic zero. For genus greater than
  one a nonconstant selfmap must have degree one, because its ramification
  divisor is effective. The same formula computes the genus from the
  squarefree odd-degree double-cover model.
- [Stacks Project, Section 109.7, tag 0DST](https://stacks.math.columbia.edu/tag/0DST):
  the automorphism group scheme is of finite type; absence of global
  derivations implies dimension zero. The local bridge is
  `deg(T_C)=2-2g<0`, giving no nonzero global tangent sections for `g>=2`.
  Its geometric automorphism group is therefore finite, and the group
  over `K(c)` injects into it.

The consumer retains an explicit invariant `c=I_a`, rather than replacing
it by `f(I_a)` and assuming geometric integrality of that different generic
fibre. It retains characteristic zero, geometric genus at least two, a
nonconstant selfmap, and a legitimate scalar-time group law. It proves
neither non-rationality for every source carrier nor a statement about all
planar Keller maps. JC(2) remains open.
