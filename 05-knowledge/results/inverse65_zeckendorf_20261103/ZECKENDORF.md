# Numeral inverses must compile into legal source operations

Status: **PROVED** restrictions on the precisely specified compiler models,
with independent **FINITE-EXACT** spectral certificates. Existence of a larger
compiler, or of a non-finite-state uniform surjectivity proof, is **OPEN**.

## Z1. What survives the older bijections

The older golden SCC word count has

    sum a_n z^n = 1/(1-z-z^3-z^4)
                 =1/((1+z^2)(1-z-z^2)),
    F_(n+1)=a_n+a_(n-2),

with negative-index a equal to0. Its actual monomer/domino parsing uses
`1 <-> a, 21 <-> b, 22 <-> c`, with a retained terminal2 flag. See
[arithmetic.md](../scc_basins_20261020d/arithmetic.md). This enumerates source
words by vertex order. It does not invert their Hamilton-path counts; multiple
marked Q4 charts and the boundary port have been forgotten in the count.

Three recovered mechanisms are useful here:

- [Native Zeckendorf register](../zeckendorf_guard_automaton_20261003.md):
  `(X,Y)->(Y+d,X+Y+2d)` with a previous-digit guard. Both registers are needed;
  the value alone does not determine the shifted value.
- [Literal golden carries](../golden_digit_carry_20261003.md): literal base-phi
  evaluation, integer Zeckendorf normalization, and orbit parity are different
  maps. A radix position/carry witness cannot be erased by calling them all
  Fibonacci representations.
- [Reflection/CG normal forms](../session131_reflection_numeration_moonshine_20261031/zeckendorf/PROOFS.md):
  parity of the sum/difference basis, carry phase, and graph boundary are
  integral coordinates. Losing them changes the image, even if a rational
  linear formula survives.

These are precedents for retaining an extra coordinate, not evidence that a
numeral digit is already a legal tournament operation.

## Z2. Prefix-correct append-only compilation is impossible

Use canonical Zeckendorf digits most significant first, with weights F2=1,
F3=2,..., no leading0 and no adjacent1. If a prefix has value n>=1, its shift
Y satisfies `n+1<=Y<=2n`: every weight strictly increases and by at most2.
At a prefix ending0, both next digits are legal and yield Y and Y+1.

Suppose a compiler's current actual source has H=2n+1 at every prefix and it
can only append bank symbols. Both child targets `2Y+1,2Y+3` lie strictly
between h=2n+1 and3h. By INVERSE(I3-I4), only the singleton a can produce a
result in that interval. One operation has one result, so it cannot provide
both children from the same retained source. Additional abstract control
states do not fix this physical branching obstruction.

This does not exclude final-only correctness, delayed output, two-ended joins,
re-marking with its proper state, or rebuilding a source after each prefix.

## Z3. Final-only finite-state compilation has a spectral obligation

Consider a deterministic finite-state transducer reading the same canonical
digits. Each transition emits a finite word of the65 bank symbols; a fixed
initial word and a terminal word depending on the final state are allowed.
Emitted words concatenate in processing order. Only final correctness
`H=2n+1` is required. There is no replacement of earlier output, no arbitrary
new source selected by an oracle, and no unbounded external counter.

On inputs `1 0^t`, the zero transitions eventually enter a cycle of length q.
Fix a phase in that cycle. Its emitted word U cannot be empty, since output
H must grow along this phase. Let P_U be its3x3 product. With fixed prefix
and terminal suffix, the H sequence is a matrix coefficient of `P_U^k` and
must equal `2F_(qk+c)+1` for a fixed c.

The latter has minimal polynomial

    (X-1)(X^2-L_q X+(-1)^q),                            (Z1)

where L_q is the Lucas number. Its three roots are distinct and its three
coefficients in the exponential expression are nonzero. Therefore its minimal
polynomial divides the degree3 characteristic polynomial of P_U and equals it.
In particular `det(P_U)=(-1)^q`. Every emitted symbol must be unimodular.

Exactly seven letters are unimodular:

    a, b, q0, r0, r20, r30, r6.

Modulo2, all their matrices are block upper triangular, with top-left1 and
lower block either I or the2x2 swap. Products have characteristic polynomial
`(X+1)^3` over F2. Comparing (Z1) implies L_q even, which is equivalent to
`3|q` by the Lucas parity recurrence. This parity/clock restriction is derived
from the actual matrices, not imposed from a Fibonacci analogy.

## Z4. Two exact exclusions of q=3,6,9

If U has ell emitted symbols, norm growth forces

    (9/5)^ell <= phi^q.                                 (Z2)

Indeed the norm of the growing actual prefix cannot exceed230 times the
final H, while final H grows at rate phi^q. For q=3,6,9, respectively,
ell<=2,4,7. Exact inequalities use phi^3<5, phi^6<18, phi^9<77 and
`9^3>5*5^3`, `9^5>18*5^5`, `9^8>77*5^8`.

The required `(trace, second characteristic coefficient, determinant)` are
`(5,3,-1)`, `(19,19,1)`, `(77,75,-1)`. The C++ verifier tests all960799
nonempty words of length at most7 over the seven letters, with each length
cutoff respected. None has the required characteristic polynomial. The only
trace candidate is q0 for q=6; its second coefficient fails.

An independent exact-rational Python certificate uses symbol-specific norm
bounds `lambda_l=min_i (P_l u)_i/u_i`:

    a:9/5, b:17/3, q0:31/3, r0:59/3,
    r20:253/6, r30:179/6, r6:155/6.

Only words with product lambda below5,18,77 can qualify. There are respectively
2,8,53 such words; their full characteristic coefficients are recorded in
`zero_cycles.json` and none qualifies. Positive lower bounds greater than1
make this pruned enumeration finite and complete.

Thus q is a positive multiple of3 exceeding9: **q>=12**. Every such final-only
compiler needs at least12 states on its eventual zero cycle. This is a lower
bound, not an existence result and not a bound on arbitrary algorithms.

## Z5. Relation to the arithmetic clocks

The classical Fibonacci shift has determinant-1 and is invertible modulo
every modulus. The quadratic lane in CHARGES.md has determinant13: precisely
at powers of13 it loses a mode and develops a transient. The correct common
object is a pair register with its legal transition, initial source, and phase.
Calling both registers a single scalar residue would erase the information
needed for either inverse. Prior norm-controlled observers are in
[Session91](../session91_packing_golden_carriers_20261029/zeckendorf/NORM_CONTROLLED_OBSERVERS.md).

This supplies a practical next test: before attempting a Zeckendorf compiler,
find a legal emitted cycle with (Z1), then test its prefix and suffix cones
and literal port interpretation. A matching characteristic polynomial alone
still would not prove a transducer or global coverage.
