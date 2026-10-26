# Independent audit of definedness and finite comparison certificates

**Status: PASS after three bounded repairs, independently verified 2026-10-26.**
Reviewer: tournament-taxonomy lane. Scope: the sibling certificate README,
`ZERO_AND_FINITE_CERTIFICATES.md`, `GALOIS_COMPARISON_GATES.md`, their exact
input, PARI replay, and `check_certificates.py`. This is a proof and source
audit in addition to a computational replay. It does not audit the sibling
quadratic-fibre-atlas report or prove a geometric interpretation of the TT7
raw functions.

## Repairs requested and verified

The first review found three narrow omissions. Their repairs are present in
the frozen bytes listed below; no outstanding mathematical correction remains
within this audit's scope.

1. The late-perturbation hostile now requires two computably specified,
   distinguishable codomain values. A singleton codomain would invalidate the
   unrestricted assertion. This does not change the integer-valued TT7 case.
2. The Collatz program now starts at `n+1`, so its input domain is all natural
   numbers while the conjecture concerns positive starting integers. Starting
   the usual map at zero would introduce an unintended divergent input.
3. The checker now separately verifies `row[5] == q*row[4]`. Its other modular
   coefficient checks did not explicitly discharge this derived column.
   All thirteen input rows already had the correct values.

The refreshed normal and optimized executions are byte-identical, both with
**1,494 explicit gates**, and equal the frozen receipt. Fresh PARI output also
matches its frozen table. These counts describe executable checks, not a
replacement for the quantified proofs below.

## 1. Computability: carriers, promises and quantifiers

The five index-set classifications are correct for an acceptable enumeration
of partial computable functions. A bounded halt predicate gives the stated
upper bounds. Triggering a constant output after a specified computation
halts reduces nonhalting to both empty domain and absence of a nonzero output;
choosing output zero or one gives the respective reductions. For an arbitrary
predicate `forall x exists s R(e,x,s)`, searching for such an `s` on input `x`
and returning zero reduces that predicate to both totality and total-zero
behavior. Thus the displayed completeness claims follow from explicit
many-one reductions, rather than from Rice's theorem alone.

Blank-input divergence is a different index set from empty domain: the
program can diverge only at zero. The graph-order argument is also sound.
A missing *fixed* input can be assigned a value by a finite program change;
this proves existence of a computable proper extension, without pretending
that a missing input can be discovered uniformly from arbitrary code.
Consequently the maximal elements are exactly the total computable functions,
while the incompatible constant-zero and constant-one functions preclude a
greatest element.

Under the promise that a program computes either total zero or bottom,
running input zero semidecides the first case. A bottom semidecider on this
promise would semidecide nonhalting through the same trigger reduction.
Removing the promise restores the termination quantifier: one observed zero
does not decide total-zero behavior. The discussion correctly permits Rice's
witness to be any function on the opposite side from bottom. This agrees
with [Simpson's Theorem 8.4](https://alexksimpson.github.io/Teaching/tinotes.pdf);
the citation to Theorem 8.5 is used only for its finite-positive-information
condition.

For *uniformly total*, exact trace evaluators with decidable equality, each
comparison at an input is a decidable predicate. Equality is therefore
Pi-zero-one, and bounded simulation supplies the claimed sharpness example.
This is a different carrier from arbitrary possibly partial evaluators.
The repaired Collatz paragraph establishes only a Pi-zero-two form for one
fixed conjecture. It claims neither completeness nor undecidability. A
nontrivial cycle is a finite counterexample certificate; failure to find
one does not exclude an unbounded trajectory.

## 2. Finite candidate separation is conditional identification

Theorem C1 is correct under all three promises: a finite list, total
evaluators, and pairwise extensional distinctness. Every pair eventually
produces a disagreement; finitely many pair searches therefore terminate.
Their union separates the list. No effective test for the promises is supplied
or needed by this conditional procedure.

Identification of an unknown target requires membership in the list. The
repaired perturbation argument disproves identification without membership:
choose an input outside the finite witness set and change its value using
the two specified codomain values. Because the witness set already separates
the list, the perturbed function agrees with one listed member on the bank
and cannot equal any other listed member globally. Thus this example is
outside the whole candidate list, not merely different from its first member.
The singleton-codomain exception is now expressly retained.

## 3. The thirteen-prime determinant has exactly the advertised scope

I compared the four raw-sum formulas and their exclusions with the inherited
`audit_hyp9391/ysums.c`: the signs, factors, and conditions
`s,t not in {0,1}, s != t` agree. All four inherited source-file hashes in the
frozen input still match their files. The producer checker rebuilds these
raw values directly rather than importing that C implementation.

The modular reconstruction is exact at the thirteen displayed primes. A fresh
PARI replay agrees with the table; independently asking PARI for the level-32,
weight-three coefficient field returns `y^2 + 1`. Thus its norm calculation
as the sum of the squares of the constant and linear coefficients uses the
correct quadratic field. The separately checked eta products and fourth
moment identity are additional controls, not assumptions needed to turn a
numerically fitted matrix into exact data.

An independent Gaussian elimination over Python `Fraction`, using ordinary
rational row normalization rather than the producer's Bareiss method, gives

    295608050652838265614106593306487371138496200704.

This equals the frozen nonzero determinant. A hypothetical global rational
linear relation among the thirteen exact functions would restrict to a
null vector of this matrix and must therefore be zero. In particular, a
relation among the four raw Y classes modulo the nine defined inherited
functions would be such a relation. The global independence consequence is
valid; it is not merely finite numerical evidence for independence.

The scope exclusions are essential and correct: this certifies no motivic
rank, purity, irreducibility, known-part subtraction or conductor. Thirteen
rows are minimal for *this thirteen-column full-rank certificate*; no
smallest-prime assertion follows. The four different values at prime 19
separate four specified raw-function names only under candidate membership.
Agreement on the thirteen primes alone cannot exclude a change at 107.

## 4. Frobenius comparison theorem and the rank-eight hostile

The all-prime statement has the right target: semisimplification, via
Chebotarev and Brauer--Nesbitt. The abelian-variety upgrade additionally uses
semisimplicity and Faltings to obtain isogeny, not integral-lattice or variety
isomorphism. These boundaries agree with
[Kedlaya, Proposition 9.4.5 and the adjacent Faltings statements](https://kskedlaya.org/weil-cohom/chapter-9.html).

The precise finite theorem was checked against
[Dieulefait--Guerberoff--Pacetti, Theorem 4.1, printed page 8](https://arxiv.org/pdf/0804.2302).
The packet retains imaginary quadratic base field, dimension two, a finite
extension of Q_2, controlled ramification, identically zero residual trace
functions, equal residual determinants, a noncubic Frobenius set, and exact
trace **and determinant** comparisons there. The conclusion is equality of
semisimplifications. The packet neither substitutes a higher dimension into
this statement nor applies it to unidentified TT7 objects. A continuous
representation of the profinite group has a stable lattice; this is where
the residual comparisons live.

The counterexample to an unrestricted dimension-free noncubic shortcut is
an actual arithmetic construction. The squareclasses `-1,2,3,5` are independent:
valuations at 2, 3 and 5 fix the last three exponents and the sign fixes the
first. Quadratic extensions of Q unramified outside these finite primes have
no other squareclass generators. Their full multiquadratic compositum thus
has group F_2^4. Each direct sum consists of eight genuine one-dimensional
finite-image Q_2 representations, so it is continuous and semisimple.

For each coordinate character, exactly four summands in either parity class
contain it. Both determinants are consequently one. Expanding the product
of four character differences gives the displayed trace-difference formula;
it is supported only at `(1,1,1,1)`, with value 16. Reduction modulo two makes
both integral representations trivial, an even stronger residual agreement
than equality of residual traces.

The noncubic proof is complete: as functions on F_2^4, homogeneous cubics span
the fourteen nonempty squarefree monomials of degrees at most three. A function
in this span vanishing away from zero and `(1,1,1,1)` also vanishes at zero.
Its only possible nonzero support is that last point, whose Boolean polynomial
has degree four. Uniqueness of the squarefree Boolean expansion excludes it.
An independent bitset elimination confirms evaluation rank 14. Independent
quadratic-residue-set enumeration, rather than Euler's criterion, confirms
all fourteen stated Frobenius masks and the omitted prime 43; its trace pair
is exactly `(8,-8)`.

This directly refutes the dimension-free reading of the abbreviated lecture
statement over Q. It does not contradict the precisely cited rank-two theorem
over imaginary quadratic fields. The critique of the lecture's
dimension-at-most-two remark is also justified by its displayed definition:
the whole vector space is automatically noncubic, and the nonzero vectors
already suffice in dimensions one and two. See
[Kedlaya, Definition 10.2.3 and Remarks 10.1.3, 10.2.5](https://kskedlaya.org/weil-cohom/chapter-10.html).

## 5. Information losses and reconstruction hypotheses

The three loss examples are proved independently of their finite controls.
The procyclic unipotent representation has invariant-space dimension one,
whereas the split representation has dimension two; their characteristic
polynomials agree on every group element. It is expressly a local/procyclic
example, so no unsupported global geometric claim is used. The two integral
involutions are rationally conjugate by the stated eigenbasis but cannot be
integrally conjugate, since their reductions modulo two have different
identity behavior.

For the elliptic example, substituting the proposed rational map gives
`Y^2 = X^3 + 4X`. The map extends to smooth projective curves, carries origin
to origin, and has degree two with the stated kernel. It is an isogeny, so
its rational Tate-module consequence holds beyond the finite point-count
tests. The rational roots of the two cubic polynomials give respectively
four and two rational 2-torsion points. Trace equality has genuinely lost
that arithmetic coordinate.

Newton's recurrence has the correct signs for `det(1-FT)`. Appending all
m-th roots of unity preserves the first m-1 power traces while changing
rank, including within weight-zero pure spectra. In rank four the paired
spectrum `alpha,mu/alpha,beta,mu/beta` gives exactly the displayed reciprocal
formula. Rank and the symplectic similitude pairing are hypotheses, not
conclusions of an integrality check. The TT7 warning retains both.

## Reproduction and frozen identity

From the certificate directory, the following were rerun independently:

```sh
gp -q replay_modular_rows.gp > /tmp/certificates.audit.modular_rows.out
cmp modular_rows.out /tmp/certificates.audit.modular_rows.out
python3 -B check_certificates.py --output /tmp/certificates.audit.normal.json
python3 -B -O check_certificates.py --output /tmp/certificates.audit.optimized.json
cmp /tmp/certificates.audit.normal.json /tmp/certificates.audit.optimized.json
cmp results.json /tmp/certificates.audit.normal.json
```

All exit successfully. Separate reviewer computations used rational Gaussian
elimination, tuple-indexed character sums, binary bitset elimination, direct
quadratic-residue sets, and explicit source-file hashing. No implementation
from the producer was imported into those computations. The mathematical
arguments above are the general certificates; the exhaustive arithmetic
controls are finite checks of their displayed instances.

| Audited file | SHA-256 |
|---|---|
| `README.md` | `1af84a3cdbde514527d994bc07fa6ca790af2c414a4db442c972ae67f0205147` |
| `ZERO_AND_FINITE_CERTIFICATES.md` | `707654e8f7202cfe7feec2c48fada85c701e51276d1c3f0d3d84703c2ef14a8f` |
| `GALOIS_COMPARISON_GATES.md` | `5249204df093a866bbf6df75160b1c5a268e4bb9e937b043bfc98f46c604af09` |
| `check_certificates.py` | `4e377f86cadc1df2000b30bdb6e8ced3de71dfdb750c3e7cd3ba13034f283f19` |
| `tt7_rank13_input.json` | `da9d5262dde64ae62a5eaaeef2cdbb6e7e20ff1e92d5de32e682faea57fdbab6` |
| `replay_modular_rows.gp` | `dfdeb47aa43b6b159888b4093bf9d4d5b50f573b3a0cc8b4cd5fbb31ea524fac` |
| `modular_rows.out` | `175b5f68e3042453168786064b27a4653b2a9bddf9e6cdcec5bc21b77a99376b` |
| `results.json` and `results.optimized.json` | `e35d5a5b0c0e0ba3aac8b842858c6df82988f110832ec028a8b3c89bf289fb13` |

The README's then-current “audit pending” phrase is a status marker for the
owner to update; it is not an outstanding proof objection. No certificate
producer file was edited by this reviewer.
