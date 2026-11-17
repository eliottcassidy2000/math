# Counting the transverse packing lattice and spending one bit

**Proved.** At every order n the full integer difference lattice
is explicit, and its rank at transverse depth k is p_3(n−3k),
the number of partitions into three positive parts. The complete
census through order nine separates the lattice from the actual
realizable profiles. Actual tournament pairs first disclose a
difference at every depth, with quadratic-size constructions.
At depths one and two a single additional bit distinguishes
prime versus composite in a fixed marked experiment. Separately,
a literal primitive three-color ×7 switch removes the common-factor
obstruction from the earlier 137-vertex switch. None of these is
a universal prime selector.

## 1. The exact lattice per order

Use G_T(x,y) with weight x per nontrivial directed odd cycle
and y per singleton; P_T(q)=G_T(2q,q) is the inherited
[Rédei–Berge specialization](https://arxiv.org/abs/2307.05569).
Write z=x−2y. The signed **integer lattice** L_n is the span of
all differences G_T−G_U with |T|=|U|=n.

If a pack uses c nontrivial cycles and s singletons, its order
has the form n=3c+2d+s with d≥0. Put j=c+d.
The complete basis of possible difference monomials is

    x^c y^(n−2j−c),   j≥1,   1≤c≤min(j,n−2j).

No linear lattice congruence is missing. To obtain a
unit-triangular integer proof, use the transitive block on
odd ell≥3 vertices with just its first-to-last arc reversed.
A cycle must use the reversed edge, so there is at most one
per pack and the directed full-ell cycle is unique. Subtract
the transitive block polynomial y^ell. Now take an ordinal
sum of c such blocks, with lengths 3,...,3,3+2(j−c), and
pad with s singletons. The product of the differences is a
signed sum of literal order-n profiles (choose each block
cyclic or transitive independently). It has the requested
monomial with coefficient one, and all other terms have
smaller j, since they use shorter cycles in some block.
Induction isolates every monomial over Z.

For degree D=n−2j and M=min(j,D), an equivalent lattice
basis is

    x z^h y^(D−1−h),    h=0,...,M−1.

The change is unitriangular in x at fixed total degree.
Multiples of z^k in this block are exactly the span of the
basis elements with h≥k: one may divide out x and z^k
over Z[x,y] since they are coprime and monic in x, or test
the first coefficients successively. So the rank is

    rank(L_n intersect z^k Z[x,y])
        = sum_j max(0,min(j,n−2j)−k) = p_3(n−3k).

For the last identity reindex j by j−k. The original
monomials correspond to the partitions (c,c+d,c+d+s)
of n. For nonnegative n, p_3(n) is the nearest integer
to n²/12 (and is zero below 3); the defining finite sum
above also gives the exact count without rounding.
The kth graded rank is max(0,floor((n−3k−1)/2)).
In the first unrevealed layer the coefficient per allowed
degree can be any even integer in the signed lattice:
x z^k y^(D−1−k) leads with 2z^k y^(D−k).
In the aggregate-deletion convention G(2q,q+t), the same
first coefficient is a multiple of 2^(k+1), since z=−2t.

This is a statement about the **signed lattice**. An
arbitrary lattice vector need not be the difference of
two actual tournaments. Any actual nonzero difference
first seen at depth k requires n≥3k+3. The next sections
count the physical profiles and construct actual pairs.

## 2. Actual order-by-order counts and the bit budget

In this table depth 0 retains P, depth 1 adds the coefficient
of z, and depth 2 adds the coefficient of z². Trailing
constant counts are omitted. All classes, not just strong
ones, are counted. Converse pairs are not identified at
input; their polynomials already agree.

| n | Isomorphism classes | Distinct P | Through first transverse read | Distinct G | Largest P bucket (distinct G) | Extra bits |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 1 | 1 | 0 |
| 2 | 1 | 1 | 1 | 1 | 1 | 0 |
| 3 | 2 | 2 | 2 | 2 | 1 | 0 |
| 4 | 4 | 3 | 3 | 3 | 1 | 0 |
| 5 | 12 | 8 | 8 | 8 | 1 | 0 |
| 6 | 56 | 24 | 29 | 29 | 2 | 1 |
| 7 | 456 | 156 | 203 | 203 | 3 | 2 |
| 8 | 6,880 | 1,728 | 3,007 | 3,007 | 10 | 4 |
| 9 | 191,536 | 50,511 | 91,288 | 91,797 | 16 | 4 |

The bits are the worst-case fixed-length information needed
to specify **G given P** in this finite universe: the lower
bound is the size of the largest bucket, and the upper bound
is indexing each bucket. This is not a claim that four
arbitrary local queries, or four bits for every order, suffice.
It is also not the number of bits needed to recover the
isomorphism class; distinct classes can share all of G.
Strong-only counts are in [summary.out](summary.out).

The finite counts use the complete [n≤8 catalogue](../actual_hamilton_multiplicities_20261021/session24-multiplicities/README.md)
and independently audited [order-nine class ordering](../blue_quotient_20261020b/README.md).
[census.cpp](census.cpp) enumerates directed simple odd cycles
anchored at their least vertex and combines disjoint supports,
counting unordered sets of cycles once. Through order nine
there are at most three cycles; three forces three triangles
covering all nine vertices. [summary.py](summary.py) changes
basis and counts profiles (not just source classes), compares
every n8 H to the archive, and checks the weighted n9
Hamiltonian mass against n!·2^binom(n−1,2). The lattice proof
is independent of these finite data.

A minimal depth-one actual pair on six vertices is

    A=(38,8,2,5,15,30),    B=(58,4,1,6,14,30),
    G_A−G_B=xz.

The census finds a minimal depth-two pair on nine vertices:

    C=(260,9,50,53,195,83,15,111,254),
    D=(36,65,266,131,15,94,157,55,251),
    G_C−G_D=xz².

Both are strong, with the same P and same first aggregate
deletion polynomial. Their difference is a replayable
literal proof, not an inference from a class count.
For n=6m, ordinal products using i copies of A and m−i of B
give m+1 distinct G's all with the same P. Distinct i's
are separated already by the first transverse read, whose
coefficient depends linearly on i. Thus the worst-case extra
bit budget is unbounded even when only the first transverse
read is needed.

## 3. Actual hidden depth grows by a physical determinant

Use the bivariate cyclic transfer of the
[previous session](../packing_hadamard_flatness_20261117/README.md).
In (z,q) coordinates write L(F,G) for the polynomial of
C3[A,B,1] when F=G_A and G=G_B after substitution x=2q+z,y=q.
It is bilinear over Z[z] (or over Q[z] as written):

    L(F,G)=qFG+(q+z/2)(Delta_+F Delta_+G+Delta_-F Delta_-G),

where Delta_± shifts q by ±1 at fixed z and subtracts the
original. At z=0 it is q times the three-term diagonal read.
Now take equal-order physical pairs (A,B) of depth k and
(C,D) of depth l. Write F_A−F_B=z^k U, F_C−F_D=z^l V.
Make two literal tournaments:

    C3[A,C,1] -> C3[B,D,1],
    C3[A,D,1] -> C3[B,C,1].

Ordinal composition multiplies the bivariate polynomials.
The difference by bilinearity is exactly

    z^(k+l) [L(U,V)L(F_B,F_D) − L(U,F_D)L(F_B,V)].

The bracket is nonzero at z=0. Let the degrees there be
u=deg U, v=deg V, b=deg F_B and d=deg F_D.
The latter two are the orders of their tournaments, and
u<b, v<d because the leading empty-pack terms cancel
before dividing by a positive power of z. In L(f,g) the first
correction to qfg has degree two lower than qfg, with
leading multiplier 2 deg(f) deg(g). In the bracket the
pure products cancel identically. At degree u+v+b+d its
leading multiplier, times all four leading coefficients,
is 2(uv+bd−ud−bv)=2(u−b)(v−d), nonzero. Products of two
corrections are lower still. Thus the construction has
**exact depth k+l**, not just a lower bound.

The output order is 2(n+m+1). Combining near-equal depths
gives pairs of order O(k²) at every depth k, starting from
the literal orders 6 and 9. For powers of two k≥2, repeated
equal splitting from the order-nine pair gives the exact
size (29/12)k²−2/3; for arbitrary k the balanced split
from order 6 gives the bound 10k²−4 by induction.
An optional final common C3[block,1,1] makes both strong
without changing depth: at the first transverse coefficient
the operator has leading term q² times the input.
[replay.py](replay.py) independently checks the depth-four
pair polynomial at order 38.

The comparison to [AMM 12592](https://sites.math.rutgers.edu/~zeilberg/EM26/AMM326.pdf)
is about information kept by an observation, not forecasting a
random coin. The repository's [orientation split](../../../01-canon/theorems/THM-3344-orientation-splitting-saves-exactly-one-dyadic-donor-bit.md)
uses one retained bit to refine an otherwise too-coarse
composition ledger. Here an observer that only sees P
and the first k−1 additional deletion sums cannot identify
the remaining hidden direction, no matter how expensive
its computation. One more reading separates the constructed
pair. Read depth, number of bits, and cost of producing the
read are three distinct measures.

## 4. Spending one promised bit on an actual prime read

Use the [previous exact two-color factor-seven switch](../fermat_flip_stern_20261117/FACTOR_SEVEN.md).
For an input profile (u,R), its lower cyclic count is

    p=3390151449 R−3390038347 u,

and its upper count is literally 7p. Uniformly using A or
B in this reader gives the same prime 701,763,272,677,
because each has (u,R)=(17,224). This is also a useful
negative control: hidden packing data cannot change a
uniform scalar comparison.

Now run the *same fixed marked cut* on the two sides of
each promised pair. Delete the indicated vertices and
append one ordinal singleton (doubles R, leaves u).
The saved bit is bit 1, counting from zero, of
the coefficient of z^k q in G(2q+z,q).
For A,B this coefficient is 6 versus 4; for C,D it is
127 versus 125. The bit distinguishes each two-element
promise; it is not a universal primality test.

| Promise | Fixed cut (zero-based) | Prime lower read | Composite mate |
|---|---|---:|---:|
| Depth 1, A/B | delete vertex 1 | 416,989,193,737 (A) | 539,035,098,309 (B) |
| Depth 2, C/D | delete vertices 1,2 | 4,241,082,742,657 (D) | 8,634,728,973,537 (C) |

The composite mates are divisible by 3. Both sevenfold
prime outputs are in genuine dyadic gaps. The literal
prime read is certified by trial division and the cut's
cover profiles by independent subset DP. The marked cut
is specified in the promise; the aggregate transverse
moment alone does not select arbitrary vertices of an
unlabelled tournament.

## 5. The primitive three-color repair now exists

The 137-vertex [previous three-color repair](../third_cover_binary_tilings_20261117/README.md)
had gcd(H,P2)>1. The new [compressed physical pair](primitive.txt)
has orders 270,948 and 271,000, and consists of repeated
ordinal blocks with literal outgoing masks. Each line
records a signed multiplicity and the claimed P1,P2,P3
of a small block. Negative copies go to the lower arm.
There are only 67 distinct listed factors, but this is
**not** claimed to be a small or optimal construction.

Replay verifies each block by subset convolution and
checks its prime-exponent vectors. The upper exponent is
exactly one greater at prime 7 for q=1,2,3, and otherwise
the same. The lower H and P2 supports are disjoint, so
their gcd is one. Color 4 fails, checked by a nonzero
modular residue. A first common cyclic wrap with a
singleton retains both the ×7 agreement through q=2 and
primitivity of its lower pair: if the original lower
profile is (u,R,S), the new reads are u+R and 4R+2(S−u).
Their gcd is gcd(u+R,S−3u)=1 by exact computation
(the first is odd). A second wrap therefore has genuine
prime-ready input, not the old forced factor 3.

This removes the fixed common-divisor barrier without
claiming that any output of **this huge new pair** is a prime.
The concrete primes in section 4 use the smaller older
two-color reader. The remaining proof task is a physical
selector that hits a prescribed prime in a repeatable
profile family, not simply another multiplicity or information
bound.

Reproduce: `python replay.py` (normal and optimized agree).
For the finite tables, from the result folder run
`g++ -std=c++17 -O3 census.cpp -o census` and
`python summary.py /path/to/checkout ./census`. Option
`--ref 1b91742ce55a202b73c80d56561ba56d00b503a6`
pins the catalogue used here.
