# A linear-size legal re-marking protocol for every recursive correction

**Status:** PROVED for the recursive ternary family; exact finite controls below.
No claim of optimal marking-graph distance or external priority. The cost is
measured in the same square-root re-marking edges as the inherited bound.

## 1. The improved bound and its contract

Let W_k be the tournament on k ternary digits, oriented at the first differing
digit around C3. Let alpha negate every digit and let C_k be its centralizer
inside Aut(W_k). Put

    N=3^k,      c_k=(3^k−2k−1)/4.

Every correction c in C_k is the holonomy of a legal based re-marking loop
with at most

    3 c_k = 3(3^k−2k−1)/4

square-root transport edges. At k=0,1 the group is trivial and the empty loop
suffices. At orders9,27,81,243 the bounds are3,15,54,174.

This improves the inherited upper bound3[1+(k−1)3^k] in
`lucas_prime_reflection_towers_20261028/FULL_HOLONOMY.md`, using the individual
triangle generators proved in Session93. It does not assert a minimum loop
length. Full anti-involution marks can have large descriptions; constructing
and storing them is charged separately. No native tournament complement move
is counted as a re-marking edge.

Source snapshot: `88c969273f2b6bc4ea2a8c5d29bb36c558eb59d2`.
Dependencies: Session93 `holonomy/FULL_CENTRALIZER_GENERATION.md` and
`ports/PORT_ORBITS_AND_LCA.md` in `session93_reversal_lucas_clocks_20261029`.
The current source already distinguishes the old loop bound from optimality.

## 2. Centralizer portraits

For a prefix w of length d<k, let r_w rotate the next digit by+1 on the subtree
with prefix w, fixing all later digits. Every automorphism has a unique
portrait of local rotations a_w in F3: on input x, output digit d is
x_d+a_(x_0...x_(d-1)). This describes its actual permutation, not just a count.

Conjugation by alpha changes the portrait to -a_(-w). Therefore c centralizes
alpha iff

    a_(-w)=-a_w,

and in particular a_(0^d)=0. Choose the representative of each nonzero prefix
pair {w,-w} whose first nonzero digit is1. The free entries are arbitrary
trits; there are

    sum_(d=1)^(k-1) (3^d−1)/2 = c_k

of them. For each representative define the paired rotation

    b_w = r_w r_(-w)^(-1).

The subtrees are disjoint, so b_w is in C_k. Applying the nonzero powers
b_w^(a_w) in **descending prefix depth** gives exactly the desired portrait.
Deeper operations happen before any ancestor rotates their address; later
ancestor rotations transport the whole subtree and leave its suffix unchanged.
Within one depth the paired rotations commute.

Thus at most c_k elementary paired rotations are needed, using exponent+1
or−1 rather than repeating a generator twice.

## 3. Every paired rotation is one actual triangle

Write w=0^a 1 z with |z|=b and let q=0^a1 0^b. The Session93 construction gives
one based triangle with holonomy b_q. Explicitly, with products acting right
to left, put

    x=r_(0^a),      y=r_(0^(a+b+1)),
    s=alpha,       t=x² alpha,       u=y² alpha.

Both t and u are anti-involutions and

    g(u,s)^(-1) g(u,t) g(t,s) = b_q,
    g(v,w)=sqrt(vw).

Choose a centralizer element h which sends q to w using local rotations
strictly above depth |w| and has identity suffix sections at that depth.
Such an h is obtained by the folded-prefix transporter: on the side1 block
map0^b to z, and impose the negated portrait on the paired side2 block.
Consequently h b_q h^(-1)=b_w exactly, without an extra suffix permutation.

Conjugate all three marks of the triangle by h. Since h centralizes alpha,
the base mark stays alpha. Naturality of the odd square root gives holonomy
b_w. This is still a three-edge loop; h is used to **construct the endpoints**,
not silently executed as an uncharged transport. Reverse the triangle to
realize b_w^(-1), also in three edges.

Concatenating these triangles in portrait order proves the bound. Every loop
returns to the same base mark, while all named ports are moved by its exact
permutation. No source weight or joint alignment is reset between factors.

## 4. Exact protocol cost under a stated source law

Under the uniform distribution on C_k, its c_k free portrait entries are
independent uniform trits, because the portrait parametrization is a bijection.
If J is the number of nonzero entries, the prescribed protocol length is

    L=3J,       J~Binomial(c_k,2/3),
    E[L]=2c_k,       Var(L)=2c_k.

The variance follows from9*c_k*(2/3)*(1/3). These are the costs of this exact
normal-form protocol; an optimal loop could be shorter.

For an explicitly chosen product portrait law q_w(0),q_w(1),q_w(2), the cost
probability generating function is

    product_w [q_w(0)+(q_w(1)+q_w(2)) z³].

A correlated original law must retain its joint coefficients rather than
being replaced by this product. The deterministic construction works for
every correction regardless of how it was sampled.

The uniform-group fixed-point law elsewhere in the session is critical
ternary branching: a local zero rotation occurs with probability1/3.
Here a nonzero paired rotation costs a triangle with probability2/3 under
another explicitly identified group portrait. Equal probabilities do not
identify this correction process with a Collatz trajectory or its clock.
The dynamics lane supplies a separate source-preserving Haar factor and
states exactly what it forgets.

## 5. Computation and limits

The program constructs every elementary triangle, verifies its actual
anti-involution endpoints and square-root holonomy, and synthesizes every
centralizer element through k=3. It samples exact portraits at k=4 and checks
the binomial cost law exhaustively at the smaller levels. Full source/port
permutations are retained. Normal and optimized outputs must agree.

With explicit permutation arrays, producing and multiplying O(c_k) marks
costs O(N c_k) elementary array operations, excluding optional adjacency
validation. Compact prefix descriptions can be stored instead; neither
model gives free access to an arbitrary full mark. The result is a sharper
legal-step upper bound, not a uniform constant-size data representation or
a proof of black-core connectivity.
