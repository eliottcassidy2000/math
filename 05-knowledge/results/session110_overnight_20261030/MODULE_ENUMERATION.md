# Module counting, shared ports, and prime composition

2026-10-30. **PROVED elementary reductions, using the classical Birkhoff submodule formula where stated.** This is a synthesis and algorithmic representation change, not a claim to have invented Burnside, Gaussian coefficients, Hall algebras or modular decomposition.

## 1. Remove the palette group from tournament enumeration

Let a vertex permutation have odd cycle type lambda, with multiplicities m_l.
Put z_lambda=product l^(m_l)m_l! and

    e(lambda)=sum_l m_l(l-1)/2
              +sum_l binom(m_l,2)l
              +sum_(l<j) m_l m_j gcd(l,j).

There are R_lambda=2^e(lambda) tournaments fixed by that permutation.
An even cycle prevents a fixed tournament, so only odd cycle types enter
Burnside. The ordinary unlabelled count is sum R_lambda/z_lambda.

Now color vertices by vectors in F_q^r and identify colors under one
global GL_r(q) action. A color matrix is classified by its row space
W<=F_q^n, with dim W<=r. Thus the combined Burnside term needs the
number of permutation-invariant subspaces W, not enumeration of GL_r(q).
The matrix rank is retained; independent gauge changes in different
components are not allowed.

### Binary semisimple formula

For q=2 every permitted permutation has odd order, so its permutation
module is semisimple. For odd d put

    f_d=ord_d(2),  c_d=sum_(d|l) m_l,

with f_1=1. The invariant-subspace dimension polynomial is

    F_lambda(z)=product_d
       (sum_(j=0)^c_d Gaussian(c_d,j;2^f_d) z^(f_d j))^(phi(d)/f_d).

Each degree-f_d irreducible factor has multiplicity c_d, and its
submodules correspond to subspaces over F_(2^f_d). The factor
phi(d)/f_d counts the distinct irreducibles of that root order.
Summing coefficients through degree r gives the GL_r palette weight.
Multiply by R_lambda/z_lambda and sum over lambda.

For AGL_r(2), replace c_1 by c_1-1 for n>0. The invariant is the
augmented row space containing the all-ones vector; quotienting by that
line removes one trivial summand. The empty case is handled separately.

The [portable checker](portable_checks/binary_module_counts.js) implements
these module weights and compares them term by term against complete
GL3, AGL3 and one-port palette groups through n=12. It also recovers
A000568 through that range. At n=4 the GL3, AGL3 and one-port counts are
188,48,756. No wall-time speedup is claimed.

## 2. Joint named ports must remain joint

Fix k independent linear functionals on the palette pointwise, with
0<=k<=r. The complete invariant is

    (W; w_1,...,w_k),  w_i in W,
    dim W-rank(w_1,...,w_k)<=r-k.

For a vertex permutation, every w_i must lie in the fixed subspace of W.
If that subspace has dimension u, the number of ordered k-tuples of rank t is

    Gaussian(u,t;q) product_(i=0)^(t-1)(q^k-q^i).

Sum only over t satisfying the displayed rank bound. This includes
dependent or zero observed words. Separate one-port marginals do not
replace this joint tuple. These ports are fixed linear functionals,
not arbitrarily chosen color points.

At k=0 this is the GL quotient. At k=r it recovers named colors with no
remaining palette gauge. The latter gives the independent termwise
control q^(r times number_of_vertex_cycles).

## 3. All finite fields: retain the primary Jordan blocks

For characteristic p, write a cycle length l=p^e d with p not dividing d.
Its polynomial satisfies X^l-1=(X^d-1)^(p^e). Each relevant irreducible
factor therefore receives a block of length p^e, not p^e separate
semisimple copies.

For a primary factor of degree f and module partition lambda, let
Q=q^f. The number of submodules of type mu is the classical formula

    product_i Q^(mu'_(i+1)(lambda'_i-mu'_i))
       Gaussian(lambda'_i-mu'_(i+1), mu'_i-mu'_(i+1);Q).

Its dimension contribution is f|mu|. In the trivial primary component,
the fixed-space dimension is mu'_1, which is what the named-port
formula requires. It is not generally the total dimension |mu|.

A primary reference for the formula and its classical attribution is
[Maglione's proof of Birkhoff's formula](https://joshmaglione.com/2023/12/17/Birkhoff.html).
The same uniformizer lifting argument applies to finite modules over
F_q[t] with residue field F_(q^f). This is a cited classical input, not
a new theorem of this session.

### Affine quotient: shorten a smallest block

The all-ones vector has a nonzero socle component in every cyclic
trivial-primary block. Choose a shortest block. Module shears from it
to longer blocks remove the other socle components, so quotienting by
the all-ones line shortens that shortest block by one, deleting it if
it had length one.

For characteristic3 and vertex cycle type(3,1), the trivial module is
J3 direct-sum J1. The affine quotient is J3. Shortening the largest
block instead would give J2 direct-sum J1 and the wrong subspace
polynomial:

    correct: 1+z+z²+z³,
    wrong:   1+4z+4z²+z³.

Likewise a single characteristic3 Jordan block J3 has one invariant
subspace of each dimension. Treating it as three trivial summands
would incorrectly produce1+13z+13z²+z³.

The current portable executable publishes the binary implementation.
The general-field statements here are proved by the primary-module
formula and the explicit quotient argument; no new sharp runtime
claim is made.

## 4. Strong components and source prefixes without losing the global gauge

Let S_lambda be the number of strongly connected tournaments fixed by
a permutation of cycle type lambda. Its unique first strong component
is a union of whole permutation cycles. Therefore

    S_lambda=R_lambda
      -sum_(empty != mu proper subset of lambda)
         [product_l binom(m_l,k_l)] S_mu R_(lambda-mu).

Use S_empty=0. Multiplying this fixed-permutation statistic by the
**whole** module/palette weight preserves globally correlated colors.
Applying a nonlinear generating-function inversion after averaging over
the palette group is generally wrong. For example the averaged AGL
series already gives a false strong count at degree2 if inverted as
though each component had an independent gauge.

For source stripping, write f=m_1. The fixed source-free count is

    B_lambda=R_lambda-f R_(lambda-1).

Exactly k successive source vertices have fixed count

    (f)_k B_(lambda-1^k),

where (f)_k is a falling factorial. Multiply by the same complete palette
weight F_lambda, not a product of separately averaged component weights.
The empty residual has B_empty=1. These formulas retain source labels,
order and the single global color action.

## 5. A finite coordinate sidecar for fixed rank

For fixed binary palette rank r, irreducibles of degree greater than r
cannot occur in a subspace of dimension at most r. Hence the only needed
divisibility counts c_d have ord_d(2)<=r. There are finitely many such d,
since each divides2^f-1 for some f<=r.

For r=3 the complete extra coordinate list is just

    c_1, c_3, c_7.

These are unbounded integer coordinates, not a constant finite number
of states. Appending a vertex-permutation cycle of length l increments
c_d precisely when d divides l. They can therefore be attached to the
existing projective gcd-response recursion without reconstructing a
palette group at each step.

The identity gcd(l,j)=sum_(d|l,d|j) phi(d) recovers these counts from
appropriate gcd responses. In general characteristic, truncated primary
heights replace the semisimple counts; the affine case must retain the
extra height needed to identify a shortest relevant block.

This proves a sufficient closed augmentation. It does not prove
minimality or a uniform speedup over the best existing A000568 counter.

## 6. Prime composition lives in a Hall-type algebra

For a finite odd abelian group G with an odd-order automorphism g,
let A(G,g) count g-invariant Cayley connection sets and P(G,g) the prime
ones. Canonical modules are g-invariant subgroups. Thus

    A(G,g)=P(G,g)+sum_(0<H<G, gH=H)
                     P(G/H,gbar) A(H,g|H).

This is the abelian specialization of the transitive-action carrier
proved in [checkpoint2](CHECKPOINT_02.md). A pair (G,g) is a finite
Z[t,t^-1]-module. In Hall notation, the recurrence is A=1+P*A.
Classical Hall multiplication counts **both** submodule and quotient
types; see [Ringel's definition](https://www.math.uni-bielefeld.de/~ringel/lectures/green/def.htm).
The inherited type-count formula alone cannot discard the quotient.

### A concrete joint-type hostile

Take G=F_p³ and g=I+J2 direct-sum J1, p odd. Put

    t=2^((p-1)/2),  B=2^((p²-1)/2),
    C=2^((2p²-p-1)/2).

There are C invariant banks, and

    P=C-2tB+t³.

To check the quotient information explicitly: the p+1 invariant lines
lie in ker(g-I). The image line has identity quotient F_p², whose prime
count is B-(p+1)t². The other p lines have a J2 quotient and prime
count zero. Among invariant planes, one is the identity plane with B
states and p have J2 with t² states. The first-prime-cut sum is therefore

    t[B-(p+1)t²] + t[B+p t²]=2tB-t³.

At p=3 the result is72 prime states of128; at p=5 it is4161600 of4194304.
In the p=5 line example, all six subgroups have identical type J1, but
their quotient prime counts are4000,0,0,0,0,0. The embedding/cotype is
the missing information.

## 7. A generating function and an honest formal q=1 limit

For G=F_p^d with identity clock, all subspaces are invariant. Write

    A_d=2^((p^d-1)/2),
    P_d=A_d-sum_(h=1)^(d-1) Gaussian(d,h;p) P_(d-h) A_h.

Using the Gaussian factorial [d]_p!, define formal power series

    Ahat(z)=sum_(d>=0) A_d z^d/[d]_p!,
    Phat(z)=sum_(d>=1) P_d z^d/[d]_p!.

Then Phat=1-1/Ahat. These are labelled Cayley-bank prime counts, not
unrestricted tournament isomorphism counts.

If one formally lets q approach1 in these expressions, Ahat becomes
exp(z) and Phat becomes1-exp(-z). Thus the resulting formal P_d(1)
alternates as(-1)^(d+1). Negative even-degree values block a literal
object-count interpretation of that naive specialization. This is a
formal analytic/algebraic limit, not a constructed field-with-one-element
category or a theorem about its geometry.

## What this unifies

The common carrier is a module or transitive action together with its
actual subobject/quotient data. This connects palette symmetry,
named-port observations, canonical tournament composition and arithmetic
ring ideals. It does not identify every application or make all open
problems mutually reducible. Ordinary counts, graph isomorphism classes,
stabilizer-weighted counts and fixed-action banks remain different
observations.
