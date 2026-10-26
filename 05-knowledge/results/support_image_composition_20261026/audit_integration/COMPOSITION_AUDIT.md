# Independent review of the composition criteria

**PASS for the stated elementary models and the four scoped applications.**
This is a proof and type review of [COMPOSITION.md](../COMPOSITION.md), not
a claim that all repository theories or open problems are equivalent. The
reviewed bytes and application-source versions are bound in MANIFEST.json.

The review requested four bounded repairs before freeze:

1. Before conjugating operations through `E:X/~→Z`, explicitly require
   their values **and legality predicates** to descend to `X/~`, and require
   `d∘E=id` on that quotient. Without this assumption the formula is not
   typed: for example, a raw operation legal at only one of two equivalent
   states has no well-defined quotient domain.
2. Describe the Hamilton path-cover splitting as a witness-to-output
   bijection. It does not alone prove associativity of an unspecified
   binary convolution. The separate associativity hypothesis remains
   necessary.
3. Identify the eight-point observer's character as the Legendre
   determinant character on `PGL2(F7)`. It is not a determinant character
   over `F8`; the inherited `F8` Frobenius joint-degree data are distinct.
4. Say that unlabelled equinumerosity alone does not specify an
   attachment-compatible encoding. It does not rule out the inherited
   SCC/Lyndon, primitive, or arbitrary ranked type bijections. The unequal
   labelled counts at order two rule out a labelled bijection, not an
   unlabelled set bijection.

The final reviewed copy incorporates these repairs.

## C1: quotient descent

The iff criterion is correct. Operation arities and all boundary transport
choices are part of the declared input. Constancy of legality on each
coordinatewise observation fibre is necessary to define the quotient
domain; constancy of the observed result on legal fibres is necessary to
define its value. Together they are sufficient by representative choice.
The operation lives on `q(X)`, so no surjectivity onto an unused ambient Y
is being assumed. Nullary operations are included. A bare equality of
object statistics proves neither condition.

For finite weighted transitions, constancy of pushed-forward transition
rows is the corresponding aggregation condition. It can hold while
individual named transitions fail to descend; the note does not conflate
those two requirements. The inherited reset/complement example is scoped
to its declared alphabet and its loss of coloured joint alignment.

## Encoding and transported operations

With the repaired quotient-domain assumption, the decoder is the inverse
of E on its image, and the displayed conjugation is well-defined.
Induction on expression trees preserves both legality and values, so
partial identities transfer with their domains. The theorem makes no
claim about a native target operation that has not been conjugated.

The sandwich application satisfies these hypotheses by S1–S3: it has an
exact image guard and a decoded joint type independent of the exceptional
triangle choices. The automorphism kernel `C3^epsilon` is correctly retained.
The support application uses the existing equivariant reconstruction map,
the new exact image guard, and decode–delete–resupport. Ordinary induced
deletion is correctly excluded. Equivariant injectivity does imply equality
of automorphism groups for the support encoding; isomorphism-class injection
alone does not imply that for arbitrary encodings.

## Weighted witnesses and realization

The pushforward description counts an explicit finite witness fibre.
Over a commutative semiring, the stated weight-preserving bijections
between **both** parenthesizations and one admissible multi-input witness
set are sufficient for associativity, including equality of outputs.
Finite distributivity proves the result. The note correctly does not
assume that arbitrary pairwise guards have this property, or that weighted
witness counts automatically equal numbers of distinct unrooted objects.

The applications retain their different acceptance semantics:

* The Hamilton path-cover expansion I2 is bijective after retaining vertex
  disjointness, exhaustive coverage, block order and actual missing
  connector arcs. Its nonnegative integer coefficients make
  `F(1)=F(0)` equivalent to absence of positive-degree witnesses. The
  directed-four-cycle control has precisely the cooperating two-chord
  monomial described in the note.
* The E1–E4 root-card operator uses neighbourhood-marked isomorphisms that
  really extend to automorphisms. Its character is multiplicative and
  intentionally cancels odd children. Zero signed weight does not mean
  there were no witnesses. Its `1/N` orbit correction is part of the
  specific enumeration target, not a generic semiring-counting rule.
* The Legendre determinant character on `PGL2(F7)` is multiplicative in
  its own action on the eight-point observer. No identification with
  root-card signs, `F8` phase, or Collatz parity is inferred.

## One common realization

The distinction between `∀L ∃x_L` and `∃x ∀L` is exact. The all-odd clock
branch has positive sources at every finite depth but its sole common
integer source is −1; it therefore supplies the claimed failure in the
positive domain. The growing-wheel boundary is a separate arithmetic
carrier, not a conjugate dynamical system.

For nested nonempty subsets of one fixed finite `X_H`, the intersection
argument is valid: choosing a disappearance stage for every element and
taking their finite maximum contradicts nonemptiness. The bound of at
most `|X_H|` strict deletions is valid, while no waiting-time bound follows.
Replacing the fixed finite source set by growing height bounds changes
the hypothesis. The fixed-height sieve cutoff and exact clock-cylinder
counts therefore support finite realization tests and not an unbounded
prime or Collatz conclusion.

No further proof correction was identified. These elementary arguments
do not need another empirical checker. The mathematical audits and exact
receipts for the concrete lanes remain separate evidence with their own
declared finite universes; this integration review does not enlarge them.
