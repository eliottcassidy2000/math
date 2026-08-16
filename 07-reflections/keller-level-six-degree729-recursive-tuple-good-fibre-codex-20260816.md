# The 729 gate is computationally smaller in the polynomial norm coordinate

**Current status: PROMOTED AS THM-3526 AFTER AN INDEPENDENT LITERAL-MATRIX
FOURIER AUDIT.**  For the fixed sporadic Keller map, one lawful fibre at target
`(1,1,1)` over `F_733` has inverse-tower dimensions
`3,9,27,81,243`.  Its sixth x-core norm is an exact squarefree polynomial of
degree `729`.  THM-3526 promotes this to the generic sixth-eliminant
degree/separability gate and proves `[Delta_6]=[2R_6]`; it remains neither an
image, irreducibility, all-level, arbitrary-map, nor general-JC result.

## Inheritance pass

- **Closest proved mechanism:**
  [THM-3525](../01-canon/theorems/THM-3525-level-five-degree243-separability-and-discriminant-square-class.md)
  turns a lawful degree-`243`, squarefree fifth-core norm over `F_251` into the
  generic fifth block-separability gate.
- **Canonical hostile:** THM-3523 already has a `243`-dimensional evaluation
  algebra, but a tower dimension does not imply that a fresh-variable
  eliminant has degree `729` or nonzero discriminant.
- **Corrected near miss:** extending THM-3525's pointwise interpolation route
  literally would spend most of its time recomputing the same transitive norm
  at 730 scalar values.  It would be exact, but it hides the polynomial norm
  structure and makes the sixth step look artificially expensive.
- **Least-used sidecar:** norm transitivity also holds over the polynomial
  base.  Keeping the fresh x-variable symbolic turns every tower contraction
  into a polynomial-valued `3 by 3` determinant.

## Feasibility scout and route change

The first timing scout used the inherited flattened inverse solver.  The fifth
inverse stage took about `81.7 s`; one scalar absolute norm took about
`0.61 s`, projecting roughly `445 s` for a 730-value interpolation bank.  This
was feasible, but it exposed two avoidable costs.

The repaired route uses the adjugate of each relative `3 by 3` multiplication
matrix and recursively inverts its determinant in the preceding tuple layer.
The five stage times fell to approximately

```text
0.0003, 0.004, 0.045, 0.57, 7.42 seconds.
```

Then the sixth cubic is normed as a polynomial, not sampled and interpolated.
The symbolic contraction itself took about `5.3 s` in the scout.  One
deliberately flattened `243 by 243` determinant remains as an expensive hostile
and took about `25 s`.  These timings are scout observations only and are not
part of the deterministic result transcript.

## Exact construction

Let `A_5/F_733` be the five-stage inverse algebra at the fixed target, and let
`q_5=(x_5,y_5,z_5)` be its universal inverse point.  The tested polynomial is

```text
E_6(X)=L(q_5) X^3 + (4-3y_5z_5) X - 2z_5,
P_6(X)=Norm_(A_5/F_733)(E_6(X)).
```

The script imports neither THM-3525 computation.  It imports only the audited
THM-3498 generic tuple definitions, then replaces every flattened inversion by
recursive cubic adjugates.  At each inverse level it records the absolute norms
of `L`, the cubic derivative, the y-coordinate chart denominator, and `x^3`:

```text
(1,3,25,593,330,455)
(2,9,663,178,54,347)
(3,27,172,118,287,557)
(4,81,192,42,634,44)
(5,243,511,85,465,723).
```

Every entry is nonzero in `F_733`, every recursive inverse multiplies back to
one, and all five graph identities `F(q_i)=q_(i-1)` hold.  The terminal leading
coefficient also has nonzero absolute norm `364`.

The polynomial-valued relative norms have exact degree ledger

```text
3 -> 9 -> 27 -> 81 -> 243 -> 729.
```

The final leading coefficient is independently checked to equal
`Norm(L(q_5))`, and the constant coefficient equals `Norm(-2z_5)`.

## Hostiles and what they establish

FLINT gives `gcd(P_6,P_6')=1`.  Its factorization has 42 factors, all exponent
one, with degree/exponent ledger frozen in the output.  The fibre is reducible;
that is not a defect and must not be promoted to a generic factorization claim.

Seven independently evaluated scalar recursive norms agree with Horner
evaluation, including all three points `730,731,732` beyond the additive
`0,...,729` interpolation window.  At `X=730`, the literal flattened
`243 by 243` regular determinant is `74`, agreeing with the recursive route.
Two deterministic hostile controls fire:

- changing the constant term from `-2z_5` to `+2z_5` negates the `X=0` norm
  because the rank `243` is odd (`438 -> 295` modulo `733`);
- deleting the degree-`729` coefficient changes the value at `X=730` (the
  truncated value is `549`, versus the true value `74`).

Ordinary and optimized replays are byte-identical.  The ascending coefficient
digest is

```text
7aba23e306b00b14b8c60c34f9762ba8b35aecac111065058dfe9d4b3f1ecd51,
```

and the semantic digest is

```text
8009c86f1c8f290829df2ba8332dc2b09929b08cbd55376f48a13acd8c2c427c.
```

## Connection contract and exact boundary

The map is

```text
243 labelled cubic blocks over the split A_5 fibre
    -> their product P_6(X) in F_733[X].
```

It preserves the total root multiset, degree, and discriminant nonvanishing,
but forgets block labels, monodromy, and image-divisor identity.  Full degree
and squarefreeness restore exactly the predicates needed here: every one of the
243 blocks remains cubic and separable, and distinct blocks are pairwise
coprime.  They do not restore irreducibility of `P_6`, identify `R_6` as an
image prime, or prove a sixth nonproperness component.

The independent audit reconstructs all `730` coefficients from `732` literal
`243 by 243` determinants by multiplicative Fourier inversion, agrees with the
digest here, and verifies that the only two possible high Fourier coefficients
vanish.  THM-3526 therefore applies the same good-reduction openness argument
as THM-3525 and makes the next THM-2582 odd-block calculation lawful:

```text
[Delta_6]
 = [N(-2R_5)][-L]
 = [8 L N(R_5)]
 = [2 R_6],
```

because `R_6=L^1699 N(R_5)` and `L^1698` is a square.  This square-class
calculation is not an `R_6` image or irreducibility statement.

## Reproduction

Run from the repository root:

```text
python -B 04-computation/keller_level_six_degree729_recursive_tuple_probe_20260816.py
python -B -O 04-computation/keller_level_six_degree729_recursive_tuple_probe_20260816.py
```

The matching script/output LF-normalized SHA-256 values are

```text
ad5d57b124f37b23fe9541c61bc4d919106b2db6c87a38632fabf3e7c076b8de
7d152cdbf720012f1dd162bfbe603ae43691717b8fb7550a2770e3d67016eba1.
```

No META-PATTERNS card is added: the recursive polynomial-norm route is a
strong method, but the evidence still consists of two representations of one
frontier extension rather than distinct mathematical threads.
