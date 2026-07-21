# The rank-two Poisson suspension leaves two rigidity gates: DC(2) and planar JC

## 1. What the supplied abstract was hiding

The simple coordinate `R=x(2-3xq)` is not an isolated guess. In the THM-1300
counterexample, shear the source by `q=y+xz/3` and reorder the target. The third
output becomes exactly `R`, independent of the sheared third coordinate.
Replacing that coordinate by the tangent momentum

```text
L=3x^2p+(2-6xq)z
```

turns the three-dimensional determinant identity into the Nambu formula

```text
{A,B}=-det d(R,A,B)/d(x,q,L).
```

This gives the second canonical pair. A single cubic coefficient equation
chooses the shift `ell=L+g`; an explicit polynomial connection correction then
makes the first Hamiltonian mate commute with the second pair. THM-2044 records
the complete formulas and proof.

This is a **symplectic suspension with only one new variable**, not the usual
cotangent lift that doubles three variables to six. It is why the Poisson
counterexample drops from three canonical pairs to two.

## 2. The exact status after reconstruction

- `PC(2)` is false: the four-variable map is symplectic and has exactly three
  points in one fibre.
- `DC(2)` remains open. The four commutative symbols do not automatically obey
  the Weyl relations after ordering.
- Planar `JC(2)` remains open. THM-2045 proves something stronger about this
  particular route: no member `R=x(a-bxq)`, `ab!=0`, has a polynomial planar
  Jacobian mate.

The last statement is a clean de-stabilization obstruction. In `s=xq`
coordinates, the planar Hamiltonian derivation raises the Laurent `x`-sector by
one. Only sector `-1` can create a constant. Its coefficient must solve

```text
((a-bs)f(s))'=1,
```

forcing a constant `f`, while polynomiality of the original `x,q` expression
requires `s|f`. The two requirements are incompatible.

## 3. The direct DC(2) experiment

Weyl-symmetric ordering was challenged rather than assumed. Its cubic Moyal
term is nonzero for five of the six canonical relations; the term-count
fingerprint is

```text
(D,R):42, (S,T):42, (R,T):0,
(R,S):3,  (D,T):165, (D,S):273.
```

So the naive quantization fails decisively. It does not follow that every
ordering fails. For `(D,R)`, the anomaly lies in the `R`-centralizer. Multiplying
it by the elementary mate `D0` and iterating once gives a 332-term corrected
symbol satisfying the exact normalized Moyal relation `M(Dq,R)=1`. This is the
first constructive foothold toward an `A_2` endomorphism.

The hard gate is simultaneous correction. Altering `D` disturbs its relations
with `T,S`; correcting those can create higher odd Moyal orders. Formal
deformation theory predicts order-by-order corrections on affine symplectic
space, but DC(2) needs a terminating polynomial operator at `hbar=1`, not an
infinite formal series.

Two attacks now look concrete:

1. Use a PBW ordering adapted to `(R,D0,L)` rather than symmetric ordering. The
   suspension is triangular in those source coordinates, so the apparent
   infinite Moyal tower may be an artifact of the wrong polarization.
2. Solve the simultaneous `hbar^2` cochain equation exactly and measure whether
   its support stays in a finite weight cone. A strict descending weight would
   prove termination; an expanding ray would be evidence against this direct
   lift.

## 4. The planar positive program

Lee--Li's current Magnus/Newton-polygon program packages a hypothetical planar
counterexample through principal, inner, and innermost polynomials. THM-2045
suggests a local test on every exposed edge: isolate the unique Laurent sector
capable of producing the constant Jacobian, then ask whether its coefficient
ODE is compatible with the semigroup of actual polynomial exponents.

For the suspension coordinate it is not. The next theorem target is to prove
the same incompatibility whenever an inner edge has two factorized endpoints
and only one resonant constant-producing diagonal. This would close a genuine
Newton-polygon stratum rather than run another bounded coefficient search.

For the Weyl floor, the Joseph-normal-form work gives a complementary sieve:
possible `A_1` counterexamples can be organized by pairs `(k,n)` and are already
excluded when the entries are coprime or one is prime. That is relevant to
`DC(1)`/planar JC, but it should not be conflated with the new `A_2` ordering
problem.

## 5. Literature and implication audit

Primary sources checked:

- Belov--Kanel and Kontsevich, stable `JC(2n) -> DC(n)`:
  <https://arxiv.org/abs/math/0512171>.
- Adjamagbo--van den Essen, the Jacobian/Dixmier/Poisson equivalence framework:
  <https://arxiv.org/abs/math/0608009>.
- Lee--Li, the current inner-polynomial/Newton-polygon program for planar JC:
  <https://arxiv.org/abs/2408.01279>.
- Han--Pan--Chen, Joseph normal forms and remaining `A_1` cases:
  <https://arxiv.org/abs/2407.11291>.

No public source containing the exact phrase or polynomial from the supplied
abstract was found in the web search. Provenance is therefore recorded as
owner-supplied; THM-2044 is an independent reconstruction and exact audit.

## 6. Incoming HYP-8810: useful analogy, not a common reduction

The newly pulled HYP-8810 proposes that planar JC and LRC(14) reduce to one
`n=12` arithmetic-progression rigidity. A source audit does not support that
strong formulation. The cited klein-S329 result is the conditional
Euler--Zariski cover-degree-three bootstrap: it reduces that stratum to pushing
the universal cubic root cover's ramification parabola to infinity. It contains
no arithmetic-progression reduction. The separate mac-mini-S137 calculation
finds Fibonacci degree pairs at the longest *Euclidean proxy* chains and
explicitly labels this a frame, not a theorem in the Abhyankar--Moh polygon
calculus. Its integer `12` is not the twelve-element AP uniqueness wall in LRC.

The honest connection is methodological. Both problems use semigroups of
admissible exponents and continued-fraction heuristics, and THM-2045 now gives
an exact JC-side model of the relevant question: the unique Laurent sector that
could create the constant Jacobian is excluded by the polynomial exponent
semigroup. Extending that sector/semigroup incompatibility to Lee--Li inner
edges could turn the current CF analogy into a theorem. Until such a map is
proved, HYP-8810 should be read as a research analogy, not a shared reduced
crux or an implication between JC(2) and LRC(14).

## 7. Assumption challenge and Tournament Analysis

Candidate tournament vertices considered were the four outputs, six bracket
obligations, Moyal anomaly supports, Laurent sectors, Newton edges, and
quantization correction stages. None preserves the operative invariant: the
symplectic condition is one antisymmetric matrix equation, while a tournament
orientation would select directions and lose the cancellations between paired
entries. Tournament Analysis is therefore not imposed on the theorem.

The challenged assumptions that mattered were instead:

- rank-two Poisson does not mean planar Jacobian;
- four classical generators do not mean a proved `A_2` quantization;
- Weyl-symmetric ordering is a test, not the only possible polarization;
- the complicated four-variable map is not independent of THM-1300 but a
  one-variable symplectic suspension of it.
