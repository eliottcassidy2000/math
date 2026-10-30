# What changes when 3 joins the seven-prime sigma code?

Session 108. **Elementary proved deductions and exhaustive finite guards; no novelty or open-problem claim.** Inherited input: Session 107 `incoming/CONCURRENT_EXTENSIONS.md`, §2, already identifies the seven-prime sigma bank as the Hamming/simplex pair. The general squareclass observer and character law are inherited from `characters/CHARACTER_STRUCTURE.md`, §3. This note supplies the zero-column extension and exact symmetry distinction, not a new discovery of classical codes.

## 1. Three named numbers, three different observers

Write sigma(n) for the sum of **all** positive divisors, and s(n)=sigma(n)-n for the aliquot sum. For primes:

| label p | binary notation of p | sigma(p) | squareclass of sigma(p) | s(p) |
|---|---|---|---|---|
| 3 | 11 | 4 | 1 | 1 |
| 5 | 101 | 6 | 6 | 1 |
| 19 | 10011 | 20 | 5 | 1 |

Binary digits of an integer are not its prime-exponent parity coordinates. The direct squareclasses of (3,5,19) are independent (rank 3). The sigma squareclasses (1,6,5) have rank 2: the 3-port is always zero as a character bit, while the 5- and 19-ports are independent. All three aliquot outputs are 1 and have rank 0. Character statements exclude primes dividing the arguments; for simultaneous direct and sigma observers it suffices to exclude 2,3,5,19.

For the joint six-port observer (3,5,19; sigma(3),sigma(5),sigma(19)), rank is 4, with latent primes (2,3,5,19). Precisely, sigma(3)'s bit is zero and sigma(19)'s bit equals direct 5's bit. The direct 3-bit XOR sigma(5)-bit recovers the missing 2-bit. The four rows direct 3, direct 5, direct 19, sigma(5) are independent; hence these two relations are complete. Thus sigma alone loses information, but joint direct/sigma observation can add information relative to direct observation alone.

## 2. The inherited seven-column system

Order the bank as B=(2,5,7,19,29,59,89). In squareclass coordinates (2,3,5), with integer masks 1,2,4 for these basis vectors, its columns are

    M = (2,3,1,4,7,6,5).

They are all seven nonzero vectors of F₂³. The selection code K=ker M has parameters [7,4,3]; the character-response code C=im Mᵀ has parameters [7,3,4]. These are different spaces: four free selection coordinates versus three free character coordinates. Their orthogonality records the statement that every square-product selection has even character parity.

For squarefree N=product over selected p in B, multiplicativity of sigma gives sigma(N)=product(p+1). Thus K exactly describes the selections giving square sigma(N), not only a formal squareclass analogy.

## 3. Adding the actual prime 3: five selection bits, still three response bits

Because sigma(3)=4 is a square, prepend a zero column to get M₀=(0|M). Then

    ker M₀ = F₂ × K,       dimension 5, minimum distance 1;
    im M₀ᵀ = {0} × C,     dimension 3, minimum distance 4.

There are 32 square-sigma selections and 8 possible character responses. This is an exact setting in which **five selection bits and three response bits coexist**. They encode different questions. The kernel is [8,5,1], not the [8,4,4] extended Hamming code. Its weight enumerator is

    1 + t + 7t³ + 14t⁴ + 7t⁵ + t⁷ + t⁸.

For any positive N coprime to 3, sigma(3N)=4 sigma(N), so square-sigma status is unchanged on toggling a single 3-factor. This coprimality condition matters: sigma(9)=13, whereas sigma(3)²=16. It is not a free toggle at arbitrary powers of 3.

## 4. Exact permutation symmetry, including pointed 19

Every invertible linear map of F₂³ permutes the seven nonzero columns and preserves K and C. Its order is (8-1)(8-2)(8-4)=168. Conversely, a coordinate permutation preserving K preserves its weight-three supports, precisely the triples {x,y,x+y}. It therefore induces an additive bijection of F₂³ after setting 0 to 0, hence is linear. Orthogonal duality gives the same group for C. These are the complete coordinate-permutation groups, not merely subgroups.

For the zero-column extension, the zero response coordinate is the unique coordinate identically zero on all responses. Equivalently the corresponding coordinate unit vector is the unique weight-one kernel word. Any permutation automorphism fixes it, and restricts to one of the above 168 symmetries. Merely adjoining the prime 3 does **not** increase the permutation group to an affine group.

The 19-column is v=(0,0,1), mask 4. Its stabilizer has order 168/7=24. There are seven Fano lines {x,y,x+y}, three through v, and four not through v. The three through v have the concrete prime-label pairs

    {19,2,59},   {19,5,29},   {19,7,89}.

Each support has square product of sigma outputs. In the first case the product is 20·3·60=3600; similarly 20·6·30=3600 and 20·8·90=14400.

The stabilizer is isomorphic to S₄, not just equal in order: the four lines avoiding v are kernels of the four linear functionals f with f(v)=1. The stabilizer acts on these four functionals. If it fixes them all, it fixes their pairwise differences and their linear span, which is the full dual space; it is therefore identity. The action is faithful, and its order 24 equals |S₄|.

Here is the justified **7/21** relationship: seven lines each contain three points, so there are exactly 21 point-line incidences (flags); equivalently seven points each lie on three lines. GL₃(2) acts transitively on flags, with flag stabilizer order 168/21=8. This does not identify the numerical label 7, the integer 21, or any unrelated 21-object construction with these incidence sets without an explicit map.

These code symmetries relabel this bank while retaining square-product constraints. They do not assert permutations of all integers preserving sigma, addition, primality, or the original named-prime identities. If every port name is fixed, the coordinate-permutation symmetry is trivial.

## 5. The separate affine completion

Index eight coordinates by every x in F₂³. The zero-extended image consists of linear truth tables x↦z·x. Translating x↦x+t sends this to z·x+z·t, usually with a nonzero constant term; it need not stay in the linear code. For example x↦x₁ becomes x↦x₁+1 under translation by (1,0,0), and its value at 0 is now 1.

To obtain affine symmetry, change the construction: add a constant-one generating row to the response matrix. The resulting code is

    A = { (a+z·x)_(x∈F₂³) : a∈F₂, z∈F₂³ }.

It has dimension 4; nonconstant functions have weight 4, the nonzero constant has weight 8. Thus its parameters are [8,4,4], with weight enumerator 1+14t⁴+t⁸. Two affine functions have intersection support of size 0,2,4, or 8, so all codewords are orthogonal; dimension 4 then makes A self-dual. This is the usual first-order Reed–Muller realization of the extended Hamming code. Puncturing its zero coordinate yields a dimension-four length-seven code orthogonal to C, hence K; the deleted coordinate is its parity extension.

Its full permutation automorphism group is AGL₃(2), order 8·168=1344. Affine bijections clearly preserve affine truth tables. Conversely any coordinate permutation preserving A pulls back each of the three coordinate functions to an affine function; its three components therefore define an affine map. Being a permutation forces its linear part to be invertible. This proves completeness.

The added constant bit is an additional mathematical degree of freedom, not an independence supplied by sigma(3)=4. The constant character of a square is +1 (bit 0); it is not a freely chosen constant bit. Even adjoining the fixed base -1 while restricting q≡3 mod4 supplies a fixed bit 1, not a new independent bit.

## 6. Arithmetic realization: add prime-count parity, not a fictitious square bit

**Parent-session deduction, verified here.** There is a concrete arithmetic meaning for the constant row in §5. For a squarefree selection N of the eight primes B₀=(3,2,5,7,19,29,59,89), observe not only the squareclass of sigma(N), but also

    epsilon(N) = Ω(N) mod 2,
    (-1)^epsilon(N) = lambda(N),

where Ω counts prime factors with multiplicity and lambda is the Liouville sign. In this squarefree bank Ω(N) is simply the number of selected ports. Set s(N) to the three-bit sigma squareclass syndrome. For the selection vector w the joint observer is exactly

    Jw = (epsilon(N),s(N)),      J = [ all-ones row ; M₀ ].

The four rows are independent: every row of M₀ is zero at the prime3 coordinate, whereas the all-ones row is 1 there. Thus J has rank4, has 16 syndromes with 16 selections each, and

    ker J = {w in ker M₀ : weight(w) is even}.

The row space of J is the affine truth-table code A of §5; self-duality gives ker J=A. Equivalently, the **squarefree products with square sigma and an even number of prime factors** form the [8,4,4] extended Hamming selection code. The dimensions change from five selection/fiber bits and three syndrome bits to four selection/fiber bits and four syndrome bits. Here a fiber's bits count the choices remaining after fixing its observer value; they are not extra independent observed quantities.

Every square-sigma selection on the original seven primes has exactly one extension with even prime-count parity: adjoin 3 precisely when its original parity is odd. Multiplication by 3 changes parity and leaves sigma squareclass unchanged, because these original products are coprime to 3. Thus the actual prime3 coordinate supplies the usual parity-extension coordinate. The result has 16 selections, with weights 0,4,8 and multiplicities 1,14,1.

### Why the larger group is not a monotonicity paradox

Let an affine permutation of the eight squareclass-labelled ports be x↦Lx+t, with L invertible. Relabel the selected prime ports accordingly, producing N'. Directly summing the selected columns gives

    epsilon(N') = epsilon(N),
    s(N') = L s(N) + epsilon(N)t.

Therefore the joint syndrome changes by the invertible linear rule

    (epsilon,s) ↦ (epsilon, Ls + epsilon t).

In particular every affine port permutation preserves the joint **zero fiber**, hence the extended Hamming code. But a translation does not preserve the old sigma zero fiber on odd selections: the single prime3 has s=0 and odd parity; translation by the 2-squareclass mask sends its port to prime7, whose sigma is 8 with nonzero squareclass2. Thus prime3 and prime7 can be interchanged as part of an affine symmetry of the new zero fiber, although prime3 was fixed by every old sigma-zero-fiber symmetry.

There is no theorem that the automorphism group of a smaller selected subset must be a subgroup of the automorphism group of its ambient subset. Here the allowed observable equivalences changed: the new group is allowed to mix parity into the sigma syndrome. It is **not** the group that preserves both the old sigma observer and parity as separately fixed observables. It also does not preserve every named joint syndrome fiber individually; the displayed rule specifies exactly how fibers are transported. If preserving the old sigma zero fiber as well is required, the translation part must vanish and the 168-element group returns.

### Character-side caveat

The dual characters of the joint syndrome are

    (epsilon,s) ↦ (-1)^(a epsilon + z·s),  (a,z) in F₂ × F₂³.

On a single prime port x, these yield the affine truth tables a+z·x, providing a natural interpretation of the constant generator. The new a is a chosen dual-character index: it is **not** an independent random bit automatically supplied by a growing observation prime q. With a fixed q and a fixed choice a, the arithmetic character is lambda(N)^a chi_q(sigma(N)); merely varying q in the inherited q≡3 mod4 model does not randomize the choice of a. In particular this arithmetic completion does not revise the earlier statement that the sigma-only prime3 character is identically +1.

## Verification

`python3 codes/verify_codes.py` and `python3 -O codes/verify_codes.py` both pass. The script exhaustively checks 5,040 seven-coordinate permutations and 40,320 eight-coordinate permutations, all GL/AGL maps, full codeword sets, weight enumerators, the pointed S₄ action, arithmetic examples, and hostile distinctions. The arithmetic completion additionally checks all 256 squarefree selections by summing explicitly constructed divisors, all 16 joint fibers, and 344,064 affine-map/selection syndrome identities. It uses explicit exceptions rather than removable assertions. These computations guard the above elementary proofs; bounded examples alone do not prove general sigma multiplicativity or asymptotic character frequencies.
