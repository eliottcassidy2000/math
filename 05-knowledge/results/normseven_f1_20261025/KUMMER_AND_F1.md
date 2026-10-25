# Marked norm-seven pairs: arithmetic reconstruction and a monomial repair

2026-10-25. **PROVED**, by the explicit algebra below. The Fano-pair
inheritance is pinned in `SOURCES.json`; the exhaustive Boolean-basis census
is **FINITE-EXACT**. No integer-square hourglass point is constructed.

## 1. A universal root-algebra statement

Let R be a nonzero commutative ring, let b=(b_1,...,b_n) be units of R,
and let m>=2. For an integer n-by-n matrix C put

    A_b = R[u_1^±1,...,u_n^±1]/(u_i^m-b_i),
    s_j = product_i b_i^C_ji,
    A_s = R[w_1^±1,...,w_n^±1]/(w_j^m-s_j).

The map w_j -> product_i u_i^C_ji is an R-algebra isomorphism **if and
only if C modulo m is invertible**. Neither m nor det C need be units of R.

Proof: both algebras are free over R on the monomials with exponent vectors
in {0,...,m-1}^n. The map sends each basis vector to a unit times the basis
vector indexed by its image under C^t modulo m. It is therefore an
isomorphism exactly when that finite group map is bijective. Nonbijectivity
misses a basis vector; coefficients in a nonzero ring cannot repair it.

An explicit inverse is useful. Choose an integer matrix B with BC=I+mH.
Then

    u_i = product_j w_j^B_ij / product_l b_l^H_il.             (K1)

These formulas prove recovery on the original ring, including its constants.
For m=2, binary B gives a particularly small parity decoder, but its
correction monomials in b must not be discarded. In characteristic two this
is still an algebra isomorphism; it is not an assertion of an etale sign cover.

Equivalently, the square with horizontal torus monomial maps C and vertical
m-th power maps is Cartesian whenever C mod m is invertible. The horizontal
map need not itself be an isomorphism. This separates root-fibre recovery
from recovery on an unmarked ambient torus.

## 2. The two norm-seven pairs over the original field

Use the inherited branch order

    (L_0,...,L_6)=(E,E+x,E-x,E+z,E-z,E+x+z,E-x-z).

Let k be any field of characteristic different from two, K=k(P^2), and
choose d_0,...,d_6 in k*. On E!=0 put b_i=d_i L_i/(d_0 L_0) for 1<=i<=6.
The six classes of b_i in K*/K*^2 are independent: a nonempty product has
odd valuation on a selected distinct line L_i=0. Consequently

    L=K(u_1,...,u_6), u_i^2=b_i,

has degree 64 and deck group G=C2^6. The constants may be arbitrary twists.

The inherited norm-seven residue maps a+2b and a+4b on the Eisenstein
lattice select these four Fano label rows, in the same branch order:

    U2=(1,2,4,5,6,3,7), D2=(1,2,4,6,3,7,5),
    U4=(1,2,4,5,7,3,6), D4=(1,2,4,7,3,6,5).

Here the integers 1,...,7 mean nonzero three-bit vectors. For each row the
three bit characters give a map G -> C2^3 and its kernel H. Each U/D pair
has complementary kernels. Taking the corresponding three root monomials
generates exactly L^H: their independent characters give degree eight and
are annihilated by H. Stacking the U and D characters gives an invertible
six-by-six matrix modulo two. Section 1 therefore proves, with every twist
retained,

    L^HU tensor_K L^HD = L.                                  (K2)

On the complement V of all seven branch lines, the same identity holds for
the finite etale cover algebras over every affine open of V, hence globally
as covers. At each k-point of V, pairs of quotient points above that **same
base point** correspond bijectively to full-cover points. Two independently
soluble quotient surfaces, without a common base point, do not suffice.

Over k of characteristic zero, the inherited free-action proof remains true
after twisting, since over an algebraic closure the twists can be scaled
away. On the finite normal branched models, normalization of the fibre
product recovers the root cover. One may either take its unique dominant
component or use the inherited finite-flat quotients to see integrality:
the product is locally free and embeds in its generic field (K2).
Do not replace those finite models by smooth resolutions inside the product.

The field and unramified algebra claims work in every characteristic !=2;
the Campedelli/surface classification is only imported in characteristic zero.

## 3. An explicit integral inverse, and why the old determinant was not intrinsic

For U2 use masks (1,2,4). For D2 first use masks (1,2,4). After deleting
branch coordinate zero the ordinary integer exponent matrix is

    C0 = [0 0 1 0 1 1]
         [1 0 0 1 1 1]
         [0 1 1 1 0 1]
         [0 0 0 1 1 1]
         [1 0 1 1 1 0]
         [0 1 1 0 1 1].

Its determinant is 3 and its Smith factors are (1,1,1,1,1,3).
Thus it is a degree-three torus isogeny, although K1 already reconstructs
every square-root fibre. Its kernel as a group scheme is mu_3; in
characteristic three this is nonreduced, not three distinct geometric points.
For a finite field F_q its torus map has gcd(3,q-1) rational kernel points.

Now change just the third D mask from 4 to 5. The last row becomes
(0,1,1,1,0,0), giving a matrix C1 with determinant 1 and inverse

    C1^-1 = [ 0  1  0 -1  0  0]
            [ 1  2 -2 -1 -2  3]
            [ 0 -1  1  0  1 -1]
            [-1 -1  1  1  1 -1]
            [ 1  1 -2  0 -1  2]
            [ 0  0  1  0  0 -1].                            (K3)

If w_4,w_6 were the old D masks 1,4, the new last coordinate is

    w'_6=w_4 w_6/(b_5 b_6).                                 (K4)

The two quotient subfields and their kernels have not changed. Equation K4
does require the base radicands. It is not an invertible reparameterization
of the old unmarked six-torus. For the inherited U4/D4 bit bases the
determinant is -3; the same arithmetic reconstruction applies.

For one fixed U2/D2 pair there are 28 unordered character bases on each
side. All 784 Boolean stacks have absolute determinant distributed as

    1:462, 3:231, 5:63, 7:21, 9:7.                           (K5)

The complete check of all 120 unordered disjoint Fano pairs gives this same
distribution for each pair (94,080 matrices, with branch zero fixed).
In particular every such pair has a Boolean unimodular choice. This is a
finite statement about seven labelled branches, not an all-rank theorem
about Boolean lifts. The number seven occurs as a presentation index in
only some choices; it is not an intrinsic obstruction of the norm-seven pair.

There is nevertheless an all-rank Laurent repair for quadratic roots.
Every invertible matrix over F_2 is a product of row swaps and elementary
row additions. Lift that product to integer elementary matrices. The result
Q is unimodular and Q=C mod 2. Writing Q=C+2D, replace each root coordinate
w_j by w_j product_i b_i^D_ji. This changes its exponent to Q and retains
its **individual** binary character. Thus the repair can even preserve the
marked character rows, and hence the two prescribed quotient groups.
The price is allowing Laurent monomials and keeping the base units.

## 4. Exactly which part is over F1

Here F1 means the pointed-monoid model: a split n-torus is represented by
the group Z^n with an absorbing zero adjoined. An integer matrix is a
monoid morphism on these groups; it is an isomorphism iff its determinant
is +/-1. Base extension takes it to the same Laurent monomial map over Z.
Thus K3 is an actual F1-torus isomorphism, and K1 is an actual Cartesian
root-fibre construction after base extension. The monoid convention and
base extension are **CITED** from Deitmar; the displayed maps are proved here.

This does not supply an F1 torification of the branched surface. Two exact
boundaries remain:

* C1 has a Laurent inverse with negative entries. On affine six-space its
  positive monomial formula collapses each coordinate axis to the origin,
  since every output monomial uses at least three distinct inputs. It is
  not an automorphism of the free positive monoid or of affine six-space.
* The substitutions b_i=d_i L_i/(d_0 L_0) contain additions such as E+x.
  The multiplicative coordinate change preserves these equations as an
  external arithmetic constraint; it does not encode addition in a free
  monoid. The inherited six-root K3 and seven-root general-type obstruction
  to torification in characteristic zero still applies.

There is a small exact six/seven comparison, with a useful hostile. For one
marked m-th-root character, exponent a is a fibre automorphism iff
gcd(a,m)=1. It can be replaced, by multiplying by a power of b, with a
unimodular one-dimensional torus exponent iff a=+/-1 mod m. Every unit
modulo 6 qualifies. Only two of the six units modulo 7 qualify; exponent
2 is a counterexample. Modulo 5 exponent 2 already fails, so this is not a
phenomenon unique to seven and is not an identification with the Eisenstein
norm-seven geometry. It is an exact map between two different operation sets.

## 5. Arithmetic consequences with the source retained

For char k!=2 the map on constant twist classes

    (k*/k*^2)^6 -> (k*/k*^2)^3 x (k*/k*^2)^3

induced by either complementary character pair is an isomorphism. K1 gives
its inverse. Over a finite field of odd cardinality, at a branch-free base
point set epsilon_i=chi(b_i). The full rational fibre has size
product_i(1+epsilon_i); the two quotient fibre sizes are the two products
of three factors 1+product_i epsilon_i^C_ji. Their product equals the full
fibre size because C modulo two is invertible. This is a pointwise equality
on a common source, not a product of independently summed point counts.

These formulas can transport local square tests and point searches without
discarding twists. They do not prove existence of a simultaneous rational
base point, an integral point, or an all-source coverage theorem.

Primary convention checked 2026-10-25:
[Deitmar, Schemes over F1, sections 1-3](https://arxiv.org/html/math/0404185v7).
His original monoids are unpointed; adjoining zero and using the contracted
monoid algebra gives the pointed convention used here.
