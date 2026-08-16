# Keller level four: the old boundary cancels again

**Status: PROMOTED as THM-3498 after independent exact/proof audit.**

The independent audit reconstructed the exposed face with a separate SymPy
factorization and nested-Horner finite-sheet evaluation.  It also replaced
the submitted finite-field engine by FLINT regular-representation matrices
and multiplicative Fourier inversion on `F_101^*`; the held-out `X=0` value,
degree `81`, derivative gcd, and coefficient hash all agree.  A literal
parity ledger retained the constant class `[2]`.  The promotion proves the
old-`L` boundary, localization, generic separability gate, and
`[Delta_4]=[2G]`; it does not promote the open image/factorization claims
listed below.

**Subsequent closure.** THM-3504 supplies the lawful squarefree specialization
requested below and proves `G` is the fourth prime image with
`S_(F^4)=V(LHJ G)`.  The remainder of this reflection preserves the exact
THM-3498 boundary before that later input.

## Inheritance pass

- **Closest proved mechanism.** THM-2582 gives the odd-block discriminant
  identity, and THM-3495 gives the fixed-map level-three normalization
  `N(H)=J/(2^35 L^7)` and square class `[Delta_3]=[-2J]`.
- **Canonical hostile example.** At the generic divisor `(L)`, one inverse
  sheet stays finite while two have half-integral pole order.  Ignoring the
  finite sheet or allowing the two divergent leading terms to cancel gives a
  false norm valuation.
- **Corrected near miss.** The older heuristic “only the newest divisor is
  odd” is not an induction theorem.  At every rung one still needs the new
  norm boundary and an actual separability/coprimality witness.
- **Least-used sidecars.** The Newton face `max(i-k)` and a modular regular-
  representation determinant avoid forming either the enormous global
  `N(J)` or a degree-81 discriminant.

## 1. Exact old-boundary valuation

Let `J` be THM-3495's primitive 66,146-term polynomial and let `N` be the
cubic function-field norm induced by the fixed Keller map `F`.  At the
generic DVR of `(L)`, write `u=1/w` on either divergent inverse branch.  The
already proved inverse-chart expansion is

```text
x=w~u^-1,       y~D/S,       z~-3(D/S)u,
```

where `D/S` is a unit.  Exact extraction from the frozen coefficient ledger
of `J` gives

```text
max(i-k)=43,
in(J)=-2^58 3^51 13^8 79^4 313^2
      *x^43(3xz-2y)^15.                                  (1)
```

There are exactly 16 face terms.  Since

```text
3xz-2y=-11D/S+O(u),                                      (2)
```

the parenthesized factor is a unit and each divergent sheet contributes

```text
v_L(J(q))=-43/2.                                         (3)
```

The hostile target `(a,b,c)=(2/27,1,1)` simultaneously has
`(c,T,S,D)=(1,1,1,1/3)`.  Its finite inverse sheet is

```text
q=(2,5/6,-7/8),
```

and exact substitution into `J` is nonzero.  Hence that sheet contributes
zero, not an unobserved positive valuation.  Therefore

```text
v_L(N(J))=-43.                                           (4)
```

This is an equality: (1)--(2) prevent divergent-sheet cancellation and the
finite hostile value proves generic unit behavior on the third sheet.

## 2. The denominator-cleared fourth norm

Put `U=A^3\V(L)`.  THM-2473 identifies `F^{-1}(U)->U` as finite etale of
degree three.  Because `J` is regular, its norm belongs to
`Q[a,b,c,L^{-1}]`.  Equation (4) therefore implies

```text
G:=L^43 N(J) in Q[a,b,c],       gcd(G,L)=1.               (5)
```

No claim is made here that `G` is irreducible, primitive over `Z`, or the
equation of a new image component.  Those require an image-multiplicity and
factorization audit analogous to THM-3495.

## 3. A full degree-81 separability witness

Over `F_101` at target `(1,1,1)`, build the first three inverse cubic
algebras.  Their dimensions are `3`, `9`, and `27`; every inverse-chart
denominator and every defining cubic derivative is a unit.  In the
27-dimensional algebra, evaluate the norm of the fourth cubic core at the
82 points `X=0,...,81`, interpolate its degree-at-most-81 polynomial, and
check an unused point `X=82` by a fresh determinant.

The resulting polynomial has exact degree 81 and derivative gcd one.  Its
ascending coefficient ledger has SHA256

```text
1c05c0fd5ee48fc2dd030aebdb9ad6ddd8185fb933eb91e7e39ff553424ef5a7. (6)
```

Thus the characteristic-zero fourth `x`-eliminant generically has full
degree 81; its 27 cubic blocks are separable and pairwise coprime.  This is a
good-reduction existence certificate, not a claim that the fibre splits
over `F_101`.

## 4. The level-four square class

Apply THM-2582's norm-product discriminant identity with outer index degree
three and odd block degree 27.  The witness in Section 3 supplies its nonzero
gate.  Using THM-3495's `[Delta_3]=[-2J]` and the first cubic class
`[Delta_1]=[-L]`, one gets

```text
[Delta_4]
 = [N(Delta_3)] [Delta_1]
 = [(-2)^3 N(J)] [-L]
 = [8 L N(J)]
 = [8 G/L^42]
 = [2G].                                                  (7)
```

Consequently `L` has even valuation in the level-four discriminant square
class.  This is the second exact old-boundary cancellation after the
level-two and level-three rungs, but it is still not an all-level induction:
the exponent 43 came from the new face (1), not from a proved recurrence.

## 5. What remains

The next decisive computation is a tractable specialization or modular
representation of `G=L^43N(J)`.  It should test:

1. whether `G` is squarefree or a proper power;
2. the generic degree of `V(J)->V(G)`;
3. whether `G` is coprime to `L,H,J`;
4. whether `S_(F^4)=V(LHJ G)` with four irreducible components.

Nothing here proves a counterexample to the Jacobian conjecture, settles
LRC(14), or identifies a physical/cohomological current.  It is a fixed-map
level-four discriminant and boundary result.

## Reproduction

```bash
python3 04-computation/keller_level_four_old_L_boundary_norm_probe_20260816.py
python3 -O 04-computation/keller_level_four_old_L_boundary_norm_probe_20260816.py
python3 04-computation/keller_level_four_degree81_finite_field_probe_20260816.py
python3 -O 04-computation/keller_level_four_degree81_finite_field_probe_20260816.py
```
