# Multiplication fibres and an exact golden content monoid

**PROVED elementary statements and scoped golden-word classification.**
**FINITE-EXACT:** small modular controls and explicit integer word products.
The unrestricted Chan–Pak spectrum and a reusable two-port return controller
remain **OPEN**. No external mathematical-priority claim is made.

The useful new construction is a complete classification in the inherited
golden bridge language. An exact universal multiplier r is realizable by
a word that multiplies the gcd of **every** integer three-coordinate
input by exactly r if and only if

    r = 5^b 9^k 21^c,   b,k,c >= 0.                  (G)

This is stricter than the previous uniform-divisibility theorem, whose
admissible moduli have arbitrary prime factors in{3,5,7}. Here the normalized
integer map is unimodular, so no further modular information is lost at any
prime. Explicit words for9 and21 complete the inherited word for5.

## Inheritance and concept board

Closest mechanisms: the exact three-coordinate single-bridge transfer in
[golden RECURRENCE G1–G4](../../golden_hamilton_transfer_20261103/RECURRENCE.md),
and [G105/Gphase](../../chanpak_55_105_golden_20261110/golden/NOTE.md).
The latter constructs divisibility by every105-smooth modulus, but its
primitive contractions are individually invertible only at their own prime.
The Smith-kernel formula is already used explicitly in
[the norm-controlled golden observer, section4](../../session91_packing_golden_carriers_20261029/zeckendorf/NORM_CONTROLLED_OBSERVERS.md).
The scalar multiplication facts below are standard elementary algebra,
included with proofs to specify the information metric, not claimed as new.
Targeted repository searches for the exact W9 word, the two displayed
matrix fingerprints, and scalar/content/unimodular reset terminology
found no earlier copy of these constructors or this attainable-multiplier
classification. This is a local provenance check, not an external novelty
claim; W5, the bridge transfer, G105 and the Smith-observer method are
inherited dependencies.

Hostile: a composite modulus still has perfectly invertible unit rows;
even a prime table's central2-by2 observer can have rank one. Corrected near
miss: a raw zero residue is not the erasure of the normalized phase.
Least-used sidecar: the complete Smith partition, rather than determinant
support alone. The live board is row fibres / central projection / CRT
phase / singular-direction alignment / exact content / physical source.

## I1. The precise prime/composite split

For m>=2 and a mod m, let mu_a(x)=ax mod m on Z/m. Put d=gcd(a,m).
Then mu_a has image d(Z/m), cardinality m/d, and every nonempty fibre has
exactly d elements.

Proof. Write a=da', m=dm', with gcd(a',m')=1. The equation ax=y mod m
is solvable exactly when d divides y. After division it fixes x uniquely
mod m', leaving d lifts modulo m. QED.

Thus the nonzero row a is a permutation exactly when gcd(a,m)=1. All
nonzero rows are permutations if and only if m is prime: a proper divisor
of a composite m supplies a nonzero nonunit. Composite tables have exactly
phi(m) invertible nonzero rows. Their unit rows are just as reversible as
prime rows; zero divisors, not every composite row, create collisions.

For a uniform input X on Z/m and known row a, every output fibre is uniform.
In base-two units the exact conditional information loss is

    H(X | aX) = log2 gcd(a,m),
    H(aX) = log2(m/gcd(a,m)).                         (I1)

These formulas specify the distribution, input alphabet and observation.
They make the user's efficiency intuition precise as a no-collision
property of every nonzero row for prime moduli. They do not establish a
global minimal-size or maximal-efficiency claim across different alphabets.

For the full row alphabet a in Z/m, the average loss is

    sum_(p^e || m) log2(p) * sum_(j=1..e) p^(-j).     (I2)

Indeed v_p(gcd(a,m))=sum_(j=1..e) [p^j divides a], and each event has
probability p^(-j). Excluding the zero row changes the average to

    [m*(I2) - log2(m)]/(m-1).                        (I3)

This vanishes exactly for prime m. Examples: modulo9, row2 is a permutation
while row3 has fibres{0,3,6},{1,4,7},{2,5,8}; modulo15, row2 is again a
permutation while row5 has five-element fibres. A central unit row therefore
cannot distinguish prime and composite moduli by its own reversibility.

## I2. The central block observes a difference

The recovered central block is C_m=[[u,v],[v,u]], with u+v=m and
u=4^(-1), v=-4^(-1) modulo odd m. Both u and v are units. Nevertheless,

    (x,y) C_m = u(x-y, y-x) mod m.                   (I4)

Its kernel is the diagonal{(t,t)}, containing exactly m elements. Its image
is the anti-diagonal{(t,-t)}, also containing m elements. For a uniform
two-coordinate input, it retains log2(m) of the original2 log2(m) bits.
This remains true when m is prime. In module notation C_m=(1/2)P_-, where
P_-=(1/2)[[1,-1],[-1,1]] projects to the anti-diagonal. The block itself
is a scaled projection, not an idempotent.

Over the integers the two row sums are m rather than zero; they are the
lifted copies of m in the user's observation. They do not create two extra
entries in the four-cell block. The central labels (m-1)/2,(m+1)/2 are
units for every odd m, since any common divisor with m would divide1.
The centre is therefore a universal odd-modulus pattern, not a primality
test. Its full table and its two-coordinate observer are different maps.

## I3. Finite coordinates for every modulus

For an integer d-by-d matrix A with nonzero determinant and Smith factors
s_1,...,s_d, integer unimodular changes of basis give diagonal coordinates.
Each change is invertible modulo every m. Applying I1 to each diagonal
coordinate gives

    |ker(A mod m)| = product_i gcd(s_i,m),
    |im(A mod m)| = m^d / product_i gcd(s_i,m).        (I5)

For uniform input the lost information is the logarithm of the kernel
size. At p^j its exponent is sum_i min(j,v_p(s_i)); at an arbitrary m
the Chinese remainder theorem takes the product of these prime-power
answers. A finite Smith list therefore answers the infinite family of
modulus questions exactly. Determinant support finds the possible bad
primes; the complete Smith list records how much is lost at each precision.

For the inherited W7 matrix, the Smith factors are(7,21,1575). Its
normalized matrix U7 has factors(1,3,225). It is invertible at every7-adic
precision but not at3 or5. Conversely the normalized W3 matrix has
determinant5*7^3. Multiplying those two locally primitive words and W5
therefore gives a normalized determinant105^3, not a globally unimodular
map. Separate local primitivity cannot be silently replaced by simultaneous
primitivity. The following construction repairs that exact limitation.

## G1. Content-preserving normalization is equivalent to unimodularity

For an integer row v write cont(v)=gcd of the absolute values of its
coordinates, with cont(0)=0. Let M be a nonsingular integer matrix. For
an integer r>=1 the following are equivalent:

    cont(vM)=r cont(v) for every v in Z^d;
    M=rU for an integer unimodular matrix U.          (I6)

The reverse direction follows because both U and U^(-1) are integral,
so each side's coordinate ideal contains the other. For necessity, the
standard basis rows show that every entry of M is divisible by r; write
M=rA. If a prime p divided det A, there would be a nonzero left-kernel
vector modulo p. Scale one nonzero coordinate to1 and lift the resulting
vector to integers with that coordinate exactly1. This primitive row v
would have every entry of vA divisible by p, contradicting cont(vA)=1.
Thus det A has no prime factor and equals+1 or-1. QED.

For the golden matrices all letter determinants are positive, so the
normalized U in (I6) necessarily has determinant+1.

## G2. Exact classification of golden content multipliers

Use the same actual seven-letter single-reversed-bridge operation as the
inherited recurrence. Its row state is(H-E,E,C), with all physical vertex
names and the literal word retained separately. The determinant list is

    s,c,q0: 1;   q1,q4: 9;   q2: 5;   q5: 21.       (I7)

**Theorem Gcontent.** For an integer r>=1, there exists a golden word w
with M(w)=rU and U integer unimodular if and only if r belongs to the
monoid generated by5,9,21.
Equivalently, writing r=3^a5^b7^c, one needs a>=c and a-c even, with
no other prime factors. This is also the exact set of universal content
multipliers in (I6).

Necessity: let n_5 count q5 letters, n_2 count q2 letters, and n_9 count
q1 and q4 together. From det M=r^3 and (I7),

    n_5=3c, n_2=3b, 2n_9+n_5=3a.

Thus a-c is a nonnegative even integer and n_9=3(a-c)/2. This proves
the claimed obstruction at all word lengths, independently of any search.

For sufficiency use the inherited W5 with M(W5)=5U5, det U5=1, and the
following explicit new words, in processing order:

    W9  = q1 s s q1 s s q1,
    W21 = q5 s q0 c q0 q5 q0 c s s s q5.              (I8)

Their vertex costs are16 and34 **added** vertices. Applying either word
to a nonempty marked source produces a tournament larger by that cost;
34 is not the total order of a tournament realizing the second operation.
No minimum-cost claim is made. Exact multiplication gives M(W9)=9U9,

    U9 = [[1451,380,2358],
          [1283,336,2085],
          [3937,1031,6398]],

and M(W21)=21U21, where the integer product before division is

    M(W21) = [[8643759306,2133974409,14984227839],
              [6911718198,1706367474,11981680263],
              [12221861547,3017337582,21186980283]].

All entries of the latter are divisible by21. The determinant list gives
det M(W9)=9^3 and det M(W21)=21^3, so both U9 and U21 have determinant1.
The independent checker also evaluates their3-by3 determinants directly.
Consequently

    W5^b W9^k W21^c

has matrix5^b9^k21^c times an integer unimodular matrix, completing the
proof at arbitrary exponents. The normalization matrices need not commute;
the word order is retained. The empty word handles r=1. QED.

This gives a simultaneous primitive105 contraction W5 W21. In contrast
to the earlier product W3 W5 W7, division by105 leaves an integer
unimodular phase map. No search at modulus105 is needed for this proof.
Its added-vertex cost is56, compared with65 for that earlier product.
These are costs of the displayed words, not optimality statements.

The prime/composite distinction must consequently be kept at the right
type. Prime moduli make every nonzero table row reversible. Nevertheless
in this operation language the composite21 has a globally primitive
constructor while the primes3 and7 separately do not. The determinant
valuation of q5 couples3 and7. Reversibility of table rows does not imply
that prime-valued operations are intrinsically the most efficient.

## G3. Procedural singular alignment rather than mass target search

The words were found by rotating surviving row subspaces with the
determinant-one letters s,c,q0. For W21, apply q5 three times. Each q5
must drop rank once simultaneously over F3 and F7, exhausting exactly the
determinant valuations required for scalar21. Between q5 letters, search
the finite pair of row-subspace orbits. The displayed stages visit1,155,
and532 orbit states before finding a suitable rotation. Each intermediate
rank-one or rank-two pair has at most13*57=741 possibilities. No assertion
of globally shortest word or minimum vertex cost is made.

For W9, after q1 the row image modulo3 is a line. The unit word s s
aligns it with the kernel plane of q1. Multiplication then gives a matrix
divisible by3. Divide that integer matrix by3; its determinant is3, so
its mod3 image is a plane. The same s s aligns that plane with the
kernel of q1. A final q1 makes it divisible by3 again. Thus the total
factor is9 with no residual determinant. Each small alignment search
visits10 states. The checker replays these subspace walks, but explicit
integer products alone suffice to certify the resulting infinite family.

This is the recursive structural pattern made precise: track the small
surviving subspace, choose a legal rotation, spend a singular letter, and
record the exact valuation already removed. The finite quotient selects
an actual word. It does not certify a bounded global state space for every
other operation language.

## G4. What is preserved, and the remaining Chan–Pak obligation

For any constructed r and every precision m, normalized output

    (vM(w))/r = vU mod m

determines v mod m uniquely, because U^(-1) is integral. Its determinant1
is a certificate valid for all m at once. The raw output, before dividing
by r, has exactly gcd(r,m)^3 preimages per attained residue modulo m.
This is an exact information budget, not an entropy heuristic.

The source is an actual marked tournament and a legal bridge word; the
target is its three-coordinate bridge observer. The map preserves every
future H response in that language. The literal word, old names and full
arcs preserve the physical source separately. Scalar r and the exponent
triple do not determine the normalized matrix U, the old graph, or a
future arbitrary-cut response. Its integer inverse need not be a legal
positive bridge word, and a normalized vector is not automatically an
actual tournament state.

In particular the constructions do not multiply H by r: the observer
readout after the word is r(vU)_1+r(vU)_2, not generally rH. They are
neither H-neutral resets nor witnesses for forbidden counts7 or21.
They do not close the seventeen-coordinate two-port observer or preserve
the81/27 automorphism-order hostile by arithmetic alone. A useful next
bridge is an actual source-changing operation whose normalized state
returns to the branch-ready predicate while retaining named joint ports.
That realization and exact-height step remains OPEN.

## Reproduction and scope

Run `python3 -B check.py` in this directory and its `-O` version. The
checker pins the inherited matrices and contraction input, evaluates
the new integer products through two implementations, checks adjugate
inverses and Smith invariants, and tests the declared small finite modular
universes. Counts for a source plus34 added vertices are not enumerated:
their matrix action
uses the inherited all-order physical bridge theorem. The row/table
controls include prime5,7 and composite9,15. The central block is audited
as a difference observer at every odd modulus through101.
