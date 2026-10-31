# Session 131 independent adversarial audit

**Verdict: PASS.** General statements below were checked by proof, not inferred
from the bounded computations. All five independent checkers use only the
Python standard library, import no producer code, retain explicit exception
gates under `-O`, and produce byte-identical normal/optimized receipts.

## 1. Integral reflection module

The blue-first forest construction gives the stated **integral**, not merely
rational, module decomposition. A blue chord closes entirely in the blue
forest. Each black forest twin supplies one sign cycle. For a black quotient
chord the two usual fundamental cycles differ from the proposed `z,Rz` pair
only by black-forest twin differences. Replacing the second vector therefore
has a triangular integer change-of-basis matrix with diagonal units.

The resulting multiplicities are

    trivial = beta_1(B),  sign = c(B)-c(Q),
    exchange = P+c(Q)-c(B).

The eigensum index is `2^exchange`. A black bridge is a sign block with index
one; a black loop is an exchange block with index two. Isolated vertices,
disconnected graphs, loops, parallel edges, arbitrary coherent orientations,
and the empty graph cause no exception. Reversing fixed edges or allowing
reflection to permute vertices is outside the theorem's hypotheses.

`check_module.py` exhausts **4,165** graphs on 0–3 vertices, allowing zero/one
blue edge and zero/one black pair at every unordered pair including loops.
It adds 360 seeded larger multigraphs and explicit boundary controls. Every
case checks cycle boundaries, the complete integral basis determinant, the
reflection blocks, and the exact eigensum index. The order-five pinned atlas
has 1 trivial, 2 sign, 10 exchange blocks, even/odd ranks 11/12 and index 1024.

The order-four illustration was independently reconstructed from all eight
fixed-path tilings, all 24 vertex relabelings, and converse. It has three
merged vertices, a blue bridge, a blue loop, and a black bridge pair. Its
module is 1 trivial plus 1 sign: index one does not imply a zero odd sector.

Coefficientwise odd CG growth in the integral block basis is injective,
boundary-preserving and equivariant, but not additive or basis-independent.
Neither arbitrary edgewise growth nor independently shifted eigen-coordinates
has those properties. The `(s,d)=(3,1)` hostile is valid.

## 2. Nonconsecutive numeration and carries

The exact `s,d` parity condition, binary gluing code, transported carry law,
and reflection formula are sound, including negative coordinates with
remainders in `{0,1}`. Contracting blue components identifies the allowed
black parity vectors with the binary cycle space of the contracted graph.

The CG and Zeckendorf automata use different digit-injection formulas but
the same homogeneous zero-input matrix. `check_numeration_certificates.py`
checks the certificates against **direct Fibonacci evaluations and forbidden
substrings**, rather than importing the producer transition recurrences:

- Both machines have 17 states; every access word and all 136 distinguishing
  suffixes per machine pass direct semantic verification.
- Exhaustive paired-word comparisons: 7,381 CG words and 21,845 Zeckendorf
  words in the declared padded finite universes.
- All 2,584 legal eight-place CG words encode exactly `[0,2583]` and have
  the stated exact Beatty shift.
- All 2,584 legal sixteen-place Zeckendorf words round-trip through signed
  even-Fibonacci coefficients. The inverse is independently recovered by
  solving the local binary constraints.

The bottom-carry theorem is valid: every positive rewrite lowers coefficient
sum, and only rewrites touching index zero lower the raw shifted value,
by exactly one. The ceiling identity makes the total strategy-independent.
The signed two-input correction belongs to `{-1,0,1}`.

The new three-step theorem and Zeckendorf corollary also pass. Direct numeral
lookups verify both signed recurrence identities on `[-100,100]`, and the
precise return classification on all 3,281 equal-parity pairs in `[-40,40]^2`.
The sharp trajectories are

    CG: (3,1) -> (8,3) -> (21,8) -> (55,21),
    Z:  (3,1) -> (5,2) -> (8,3) -> (13,5).

In both languages compatibility is either permanent or returns exactly at
multiples of three. This is a parity-gluing theorem in chosen coordinates,
not a graph-wide edge map or an additive endomorphism.

## 3. Hamiltonian packing coefficients versus binary digits

The modulo-16 derivation is sound. Distinguish a literal odd cycle orbit from
an orbit of disjoint two-cycle packings. Since all reflection-fixed odd cycles
meet the unique fixed vertex, fixed two-packings are precisely disjoint
mirror pairs. Consequently

    alpha_1=q+2i+2j, alpha_2=j+2p, alpha_3=t (mod 2),
    H=1+2q+4i+8(j+p+t) (mod 16).

The actual third digit is `(floor(q/2)+i) mod 2`; the actual fourth includes
`floor((q+2i)/4)`. Reducing counts to parities before adding loses carries.
The exact identity through order eight and the even-order specializations
follow from the declared packing-size bounds.

`check_tournament_carries.py` enumerates literal Hamiltonian permutations and
literal rooted cycle permutations independently for all **599** standard-
reflection-compatible tournaments through order six. It also reconstructs
all eight distinct named hostile fixtures, including the order-nine triangle
stack. All coefficients, reflection charges, digits and congruences match.
This is not an exhaustive labelled order-seven census.

The alpha-to-H noninjectivity and old path-formula correction were inherited,
not newly discovered. Exact H=7 is not the same predicate as H=7 modulo 8.

## 4. Square-sum seams, path ports, and root words

The eight-residue lower bound `(4,7,8,7,8,7,8,7)` is valid for every existing
Hamiltonian cycle on `[1,N]`, `N>=8`. Endpoint incidence counting remains
correct at singleton paths. The residue-four strengthening is essential:
`t=4` forces one path in each block and leaves the residue-2 endpoints of
block B unable to join another block. They would self-close, contradicting
connectedness. Thus the next possible value is eight.

The endpoint-aware Hamiltonian-path equations preserve the two original
ports with multiplicity. The proof does not assume that a local path forest
closes into one cycle. The stated affine-preserver classification is sound:
integer square-root differences force a square slope, then the asymptotic
integer offset forces zero translation.

`check_square_sum.py` independently verifies:

- all eight N=80,...,87 sharp witnesses, with no optimizer trust;
- all 668 single-edge-deletion paths from those witnesses;
- the inherited G32 cycle (16 seams), singleton endpoint incidences, and its
  32 additional deletion-path controls;
- all nine long witness root words against direct walking;
- 27,200 bounded root-word cases, checking the complete closure/bounds/
  injectivity compiler against direct walking;
- the exact eight H40 path components, four components at every N=44,...,256,
  and earlier odd-square neighbors for every arrival 13,...,2000.

Sharpness is **residuewise uniform**, witnessed at eight orders. It is not
attainment at every N in those residues, nor a new all-order Hamiltonian
existence theorem. Root reflection is partial and not generally an
endomorphism of the entire square-sum graph.

## 5. Moonshine arithmetic: an independent algebraic path

The group-order/index arithmetic, supplied-D distinction, and register-role
versus prime-label boundary are correct. The actual 2A stabilizer is `2.B`.
The recovered Sylow argument does not assert that every order-p element
fixes a point. The first-six support, Baby/Monster filtration and level-two
parallel-edge pairs were already in the repository.

Primary sources independently inspected:

- [ATLAS Monster](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/M): order and
  the `2.B` subgroup index.
- [Duncan–Ono](https://arxiv.org/html/1411.5354v2): Ogg's all-supersingular-
  j rationality predicate and its fifteen-prime list.
- [Aricheta](https://sigma-journal.com/2019/007/sigma19-007.pdf): Theorem 1.1/
  3.1, Table 1 and Appendix A. Preserve the specified Atkin–Lehner quotient;
  the level-six full quotient is not merely the Fricke quotient. Bad-level
  characteristics are excluded, not silently proved rational.

`check_moonshine_polynomials.py` uses neither class numbers nor Fp2 point
enumeration. For every prime 5<=p<=199 it tests the coefficient of `x^(p-1)`
in `(x^3+a*x+b)^((p-1)/2)` for each canonical j-model, matching all 44 producer
point-count sets. The coefficient equals trace modulo p; Hasse's bound makes
trace-zero equivalent to this criterion in the tested characteristic range.

For each of the thirteen Ogg primes at least five it then independently
eliminates u from

    (u+256)^3-j*u^2=0,   u^2-s*u+4096=0.

Reduction gives `A(s)u+B(s)`, where

    A=s^2+(768-j)s+192512,
    B=4096*(j+3328-s).

The norm `B^2+s*A*B+4096*A^2` is the elimination polynomial. Squarefree lcm
across all supersingular j gives the complete finite s-locus. Gcd with
`s^p-s` gives its rational roots; gcd with `s^(p^2)-s` certifies every root
is already in Fp2. All total/rational counts match the producer, and the
only irreducible quadratic remnants occur at 29,41,59,71, each once.
This independently closes the finite field-enumeration completeness gate.

Frobenius permutation modules have no isolated sign block and even rank at
least odd rank. Hence the n5 metagraph module cannot be such a permutation
module even over Q. The exchange-block mechanism transfers; a Monster
representation or global dynamical conjugacy does not.

## Reproduction and final boundary

Run `python3 verify_independent.py`. It replays each independent checker in
normal and optimized Python, compares their saved receipts, and verifies
that every intentional failure control exits one. Inputs are sibling packet
receipts or data files used only as comparison targets, not executable
producer imports. `INPUT_PINS.json` pins the reviewed proof and data bytes;
`SHA256SUMS` pins the independent packet itself.

No mathematical defect remains in the reviewed final packets. No repository
publication was performed by this audit agent; publication belongs to the
parent session. No global LRC, arbitrary square-sum existence, moonshine
explanation, or canonical graph-growth claim is inferred.
