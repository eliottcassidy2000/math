# Arithmetic braids II: exact bridges and the coordinates they retain

**Status: PROVED scoped elementary statements; CITED background;
FINITE-EXACT enumerations; independent audits.** The original Collatz
conjecture, completeness of signed cycle catalogs, twin primes, Goldbach,
and LRC(14) remain **OPEN**. No claim of literature priority is made.

This continues the [first arithmetic-braids session](arithmetic_braids_20260917_synthesis.md).
The user confirmed floors in all three sums and retained the triangle family
`(k^2-1,2k,k^2+1)`. The earlier triangle, primitive-divisor and quadratic
three-cycle results remain in that synthesis; this pass follows the new
signals into inverse reachability, residue rings, and Fano multiplication.

## Inheritance and portfolio

Anchor: the user's reverse Collatz reachability question, especially the
three known `3n-1` basins. Niche: the three floor sums and squarefree
density. Wildcard: whether the divisor family supplies a genuine octonion
and tournament construction.

The closest proved mechanisms are the prior inverse-fibre odometer and
ordered affine carry, the divisor exponent-box classification, and
[THM-2422, paired operation fibres](../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md).
Canonical hostiles are long exponent-one growth words, squarefree composite6,
and distinct cycles sharing length and total exponent. Corrected near misses
are residue coverage mistaken for all-integer coverage, and a scalar H-count
mistaken for its prime-factor decomposition into SCC counts. Least-used
sidecars are ordinary height, zero-residue boundary, parameter gcd, and XOR
output labels.

The six live concepts throughout were: finite halving word and carry;
height versus local residue address; squarefree boundary; parameter ideal;
prime squareclass; and oriented multiplication with its XOR labels.

## 1. Any finite halving word can be completed to any admissible target

Let `T_b(n)=(3n+b)/2^v2(|3n+b|)`, with odd b coprime to3. Given any
finite word `k1,...,kL` and any nonzero odd target u coprime to3, there are
infinitely many same-sign source integers whose first L halving exponents
are that word and whose next odd step lands on u.

Write `K=sum ki`, `Ki=sum_(j<=i) kj`, and
`B=sum_(i=0)^(L-1)3^(L-1-i)2^Ki`. Choose `k0=1 or2` with `2^k0 u=b mod3`.
The sources have the explicit form

```text
n_t=[2^(K+k0)4^t u-b(2^K+3B)]/3^(L+1).
```

Exactly one t-class modulo `3^L` makes this integral. Lifting t through
that class visits every source residue modulo every `3^s` while retaining
the binary word cylinder. The proof is the exact identity
`v3(4^m-1)=1+v3(m)` plus backwards odd integrality.

For `b=-1`, the basins of `(1)`, `(5,7)`, and
`(17,25,37,55,41,61,91)` are disjoint, yet each meets **every odd class
modulo `2^H3^s`** and realizes every finite halving word. This proves a
sharp limitation on any fixed residue/prefix classifier. It also applies
to the plus basin of1 without assuming the Collatz conjecture.

The family grows exponentially in t. Its 2-adic limit is a preimage of
`-b/3`, the zero-numerator singularity, and has 3-adic valuation `-(L+1)`.
One sequence can therefore escape in ordinary height, converge to a 2-adic
singularity, and cover every ternary residue in its parameter. Keeping
those three notions of size separate is essential.

[Proof and exact constructor](arithmetic_braids2_20260917_inverse_completion.md).
An independent brute-exponent referee checked 108,000 trajectories and
16,000 singular-limit identities without importing the constructor.
For historical counting context, [Krasikov--Lagarias](https://arxiv.org/pdf/math/0205002)
prove an `X^0.84` lower bound for ancestors of each fixed positive target
coprime to3 at sufficiently large X. Their height-bearing inequalities are
stronger than the logarithmically thin completion subfamily here.

## 2. The first two floor identities characterize squarefreeness

For every `n>=2`, let `R3(n)` count residues `x modn` with `x^3=0`,
including zero, and `Z2(n)` count pairs `1<=i,j<n` with `ij=0 modn`.
The three exact identities are

```text
sum_(k=1)^(n-1) floor(k^3/n)
  = (n-2)(n-1)(n+1)/4 + (R3(n)-1)/2;

sum_(k=1)^((n-1)(n-2)) floor((nk)^(1/3))
  = (3n-5)(n-2)(n-1)/4 + (R3(n)-1)/2;

sum_(i,j=1)^(n-1) floor(ij/n)
  = (n-2)(n-1)^2/4 + Z2(n)/2.
```

The first two corrections vanish iff n is squarefree; the third vanishes
iff n is prime. The mechanisms are antipodal residue pairing and counting
the common boundary of a lattice rectangle. No cube-permutation assumption
is needed: cubes modulo7 do not permute the six nonzero residues.

This is a direct algebraic bridge. The cube identities test whether `Z/nZ`
is **reduced**, meaning it has no nonzero nilpotents. The product identity
tests whether it is a **field**. The squarefree composite6 separates them.
The first two predicates consequently hold on a set of density `6/pi^2`.

For `n=prod p^a`, `R3(n)=prod p^(a-ceil(a/3))`. At the three divisor-balance
patterns `p,p^3,p^2qr`, the cubic corrections are respectively
`0,(p^2-1)/2,(p-1)/2`. The common link is the exponent profile; these are
not three copies of the same zero-defect condition.

[All-modulus proofs and odd-degree extensions](arithmetic_braids2_20260917_floor_reciprocity.md).

## 3. Every finite halving word survives the squarefree restriction

Relative squarefree densities in the odd rows1,3,5 modulo6 are
`9/pi^2,6/pi^2,9/pi^2`. Conditioning on any specified finite halving word
leaves these constants unchanged: the word fixes only a dyadic congruence,
where oddness already removes the square obstruction at2.

There is a stronger joint theorem. For every fixed L, a positive-density
set of q makes all `L+1` nodes

```text
n_j=3^j 2^(L+1-j)q-1,         0<=j<=L,
```

squarefree. Each step grows, with exact halving exponent1. The local
excluded-root counts are zero at2, one at3, and
`min(L+1,ord_(p^2)(3/2))` at primes at least5. CRT followed by a proved
tail estimate gives the positive infinite Euler product.

Hence both cubic floor equalities can hold at **every node of an
arbitrarily long growing Collatz segment**. Squarefree status alone cannot
force descent in any fixed number of steps. This does not rule out a richer
joint invariant, nor does it prove a squarefree terminating completion:
the two separate existence constructions cannot simply be intersected.

[Density proofs, exact bounds, and joint witnesses](arithmetic_braids2_20260917_squarefree_symmetry.md).

An explicit exceptional lifting example is
`3^11-2^11=175099=23^2*331`. Since11 is prime and `3/2` is not1,
its order is exactly11 both modulo23 and modulo529. Therefore the excluded
square roots at23 stop increasing after eleven prefix nodes. Equivalently,
the affine map `x->(3x+1)/2` modulo529 has order11, after translating
`x+1` to expose multiplication by `3/2`. This concerns the exponent-one
branch; it is not a period claim about unrestricted integer Collatz orbits.

## 4. The parameter controls cycle content; order controls the carry

Negating state and parameter conjugates `3n+b` to `3n-b` while preserving
forward arrows. Reversing arrows is a separate, multivalued operation.

For `3∤b`, `gcd(n,b)` is conserved. On a cycle, even for general odd b,
this gcd equals the common divisor of all cycle nodes. Removing it gives
a primitive parameter. For a cyclic word with `D=2^K-3^L`, its reduced
rational denominator is

```text
q=|D|/gcd(B,|D|).
```

The word occurs as a signed integer cycle at parameter b **iff `q|b`**,
and its node content is `|b|/q`. Thus parameter-five cycles split exactly
into scaled parameter-one cycles and primitive denominator-five cycles.

There is a ring-theoretic connection to the floor defects. Modulo `|b|`,
a step multiplies n by the unit `3*2^(-k)`, so its principal ideal is
preserved. Closing a word requires `nD=0`: D lies in the annihilator of n.
The gcd and rational denominator encode this condition. The carry B is
still needed for exact integer closure; annihilation alone is insufficient.
An even finer sidecar is the coset of `<2,3>` occupied by `n/d` modulo
`|b|/d`. At b=23, this subgroup is the quadratic residues, of order11.
Thus residue-square class is conserved even among states of equal gcd.
The 23-adic growth coincidence above is a different, compatible statement:
the ratio `3/2` has exceptional order already modulo23 squared.

At `b=-5` we exhibit nine cycles. Two primitive negative cycles have
period17, total halving exponent27, and denominator5, but different carries
`189900931` and `352383011`; their least absolute nodes are -187 and -347.
The first word census missed them because its declared length bound was10.
A separate forward census exposed them and their complete exact witnesses
are retained. No list is advertised as globally complete.

The sole positive-to-negative crossing is `1 -> -1`, and -1 is fixed.
Thus the other negative cycles cannot be reached from positive starts.
The three known positive cycles are fivefold copies of the three known
positive `3n-1` cycles. By the completion theorem, the four exhibited
basins meeting the positive integers each cover all odd binary/ternary
residue classes; the modulus5 gcd distinguishes the portal basin from
the three scaled basins.

[Parameter theorem, sign geometry, and frozen cycles](arithmetic_braids2_20260917_signed_cycles.md).

## 5. The third divisor family contains an actual Fano geometry

The seven proper nontrivial squarefree divisors of `p^2qr` are
`p,q,r,pq,pr,qr,pqr`. Multiplication modulo rational squares identifies
them with the seven nonzero vectors of `F2^3`. Its seven triples
`{a,b,a xor b}` partition the 21 unordered pairs.

A declared Cayley--Dickson octonion basis supplies a sign to each product
`e_a e_b=+-e_(a xor b)`. Orient `a->b` when the sign is positive. This is
an actual tournament, each Fano line becoming a directed triangle.
The map from squareclasses supplies the grades; it does not turn the
commutative integer product into octonion multiplication. Signs and the
nonassociative multiplication law are additional data.

The 128 basis sign choices give sixteen orientations, since eight sign
characters change no edge. All sixteen have 189 Hamiltonian paths.
Independently orienting the seven lines gives128 orientations: the other
112 have171 paths. These exhaustive counts have three checks: dynamic
programming, permutation enumeration, and the full odd-cycle formula.
For a gauge tournament there are14 triangles,42 five-cycles,24 seven-cycles,
and seven disjoint odd-cycle pairs, giving `1+2*80+4*7=189`.

The seven planted Fano triangles are only a small part of the odd-cycle
structure; using only them would incorrectly give15. XOR-preserving
relabelings split the sixteen gauge tournaments into two classes of eight,
exchanged by complete reversal, while unrestricted relabeling merges them.
The lost XOR labels therefore erase an exact two-way distinction.

The numbers seven and21 here count points and pairs. Their being forbidden
Hamiltonian counts elsewhere is a different statement, governed by
[THM-1370-h-spectrum-omits-7-21-all-n](../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
The stale all-`7*3^k` corollary in THM-200 was withdrawn with its H63
counterexample and preserved correction lineage.

[Full Fano construction and declared gauges](arithmetic_braids2_20260917_squarefree_symmetry.md).

The gauge calculation extends to an exact coding and lattice bridge.
Label each Fano line by its nonzero normal `a in F2^3`, and let `g(x)`
record a changed sign at point x. Its line reversal pattern is

```text
(Bg)(a)=sum_x g(x) + a dot sum_x x g(x)  (mod2).
```

These are precisely the restrictions of affine binary functions to the
seven nonzero points: the Hamming `[7,4,3]` code, with sixteen words.
Every one of the112 nongauge line patterns has a unique one-line repair
into this code coset. Appending the parity coordinate gives the doubly-even
self-dual `[8,4,4]` code. Construction A gives the even unimodular lattice

```text
Lambda = {z/sqrt2 : z in Z^8, z mod2 in the extended code},
```

explicitly isometric to `D8+`, the E8 lattice. Its240 norm-two vectors
split into16 signed coordinate roots and224 roots on fourteen weight-four
code supports. The source of this eight-dimensional object is now an
explicit construction, not a numerical resemblance.

[Proof, unique-repair algorithm, and E8 isometry](arithmetic_braids2_20260917_fano_code.md).
The240 roots and the earlier240 labelled copies of the tournament are
different carriers; no bijection follows from equal counts. With coordinate0
identified as the scalar unit, this normalized lattice is not closed under
octonion multiplication: it contains `sqrt2*1` but not its square `2*1`.
Bott periodicity
requires further representation/topological data, and no global Collatz
invariant has been transported through this chain.

The role of three bits is itself exact. In binary dimension r, the same
incidence/parity construction gives `RM(1,r)`, with length `2^r` and dimension
`r+1`. Self-duality requires `2(r+1)=2^r`, which among r>=2 holds only at3.
At r4 the resulting length16 code has dimension5 and its Construction-A
lattice has covolume8, rather than being unimodular. This provides a
decisive boundary for the eight-dimensional construction.

## 6. Transfer ledger and remaining research obligations

| Source -> target | Map / preserved predicate | Information that must be retained |
|---|---|---|
| Halving word -> terminating path | Affine carry plus principal-unit logarithm | Original integer and ordinary height |
| Cube graph -> root graph | Transpose finite lattice rectangle | Equality-boundary points |
| Floor defects -> residue ring | Count cube-zero or product-zero fibres | Reducedness differs from fieldness |
| Finite word -> squarefree inputs | Fixed dyadic cylinder plus squarefree sieve | Intermediate nodes require joint conditions |
| Signed cycle -> rational cycle | Divide state by parameter | Reduced denominator and ordered carry |
| Proper squarefree divisors -> Fano points | Prime exponents modulo two | Exponent magnitudes and properness |
| Octonion table -> tournament | Positive product sign | XOR output, identity, squares, sign gauge |

The next productive Collatz target is a **height-controlled reverse cover**,
not another proof of local residue fullness. A proposed invariant should be
tested against the three dense minus basins and against squarefree growth
segments of arbitrary finite length. A successful connection must explain
how its retained coordinate behaves on a fixed integer's actual trajectory.

The reproducible [manifest](arithmetic_braids2_20260917_manifest.json) records
all run universes and LF-normalized hashes. Each finite table is evidence
only within its declared universe; universal statements have separate proofs.
