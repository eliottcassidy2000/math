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
  route: no member `R=x(a-b x^r q^s)`, `ab!=0`, `r,s>=1`, has a polynomial
  planar Jacobian mate.

The last statement is a clean de-stabilization obstruction. Grade `x^i q^j`
by `s*i-r*j`; the planar Hamiltonian derivation raises this grade by `r`.
Only one sector can create a constant, and it consists of
`q f(x^(r/g)q^(s/g))`. Its top coefficient after applying the derivation is
`-b(r+1+(s/g)deg(f))lc(f)`, which cannot vanish. This covers a full
two-exponent Newton-edge family, not only the supplied coordinate.

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
symbol satisfying `M(Dq,R)=1`. The `(S,R)` anomaly is lower order and dies after
one analogous correction, giving an 85-term `Sq` with `M(Sq,R)=0`; meanwhile
`M(T,R)=0` already. Thus the entire `R`-column is exact. The next residual
`M(Sq,T)-1` has only 59 terms and momentum degree two, but its correction must
preserve all three repaired relations.

The adapted-centralizer calculation pushes further. A correction `f star ell`
with explicit `f=x^6F(xq)` kills the tangent-degree-two layer, and two solvable
weight ODEs kill tangent degree one. The surviving constant occupies weights
six and twelve. The homogeneous mode at weight `m` then creates weights `m`
and `m+6`; exact rungs `6,12,18` produce `{6,12}`, `{12,18}`, `{18,24}` and no
linear combination closes the original constant. This is the first concrete
sign that the formal quantization may be forced into an infinite tower. It is
not a no-quantization theorem: simultaneously moving `T` or allowing higher
tangent order could escape the fixed-`T`, tangent-linear ansatz.

The hard gate is simultaneous correction of the coupled `T`-column and
`M(Dq,Sq)=0`. Altering one dual coordinate disturbs the others and can create
higher odd Moyal orders. Formal
deformation theory predicts order-by-order corrections on affine symplectic
space, but DC(2) needs a terminating polynomial operator at `hbar=1`, not an
infinite formal series.

Three attacks now look concrete:

1. Use a PBW ordering adapted to `(R,D0,L)` rather than symmetric ordering. The
   suspension is triangular in those source coordinates, so the apparent
   infinite Moyal tower may be an artifact of the wrong polarization.
2. Prove or refute the observed six-weight propagation for all homogeneous
   modes. A nonzero `m+6` coefficient at every rung would rigorously exclude a
   finite fixed-`T`, tangent-linear repair; a cancellation identifies the first
   escape rung.
3. Solve the simultaneous `hbar^2` cochain equation with `T` allowed to move.
   The current expanding ray may be a gauge artifact of freezing the quantum
   position pair rather than a true obstruction to `A_2` quantization.

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
