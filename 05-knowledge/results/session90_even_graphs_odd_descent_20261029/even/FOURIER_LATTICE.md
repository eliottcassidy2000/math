# The integral cost of the decorated Fourier correspondence

**Status:** PROVED from the decorated Walsh theorem; FINITE-EXACT tests stated below. This extends the undecorated lattice statement already present in THM-4721 to arbitrary finite decorations and the paired/self-converse carrier. No external-priority claim is made.

Dependency: [DECORATED_EVEN_GRAPH_DUALITY.md](DECORATED_EVEN_GRAPH_DUALITY.md). We retain its Royle sign, finite G≤S_n, finite G-set D, and relatively even graph-decoration definition. Put m=binom(n,2). Let A index G-orbits on tournaments × D and O index relatively even G-orbits on graphs × D. Their common number is r.

## 1. Two integral lattices, one dyadic change of coordinates

Choose one signed orbit sum Phi_O for each surviving dual orbit. Its overall sign is a free basis convention. Form the square integer matrix

    M[A,O] = Phi_O(A).

The columns are invariant integer-valued functions on the source objects. Thus M maps the free integer lattice on surviving dual orbits into the free integer lattice on source orbit indicators.

**Theorem.** M is nonsingular, its determinant has absolute value 2^k for an integer k≥0, and its integral cokernel is killed by 2^m. In particular every Smith invariant is a power of two dividing 2^m.

**Proof without hidden odd divisions.** Inverting the Walsh transform on the labelled tournament coordinate gives, for any integer-valued invariant function f,

    Fourier coefficient at (X,d)
       = 2^(-m) sum_{labelled tournaments T} f(T,d) chi_X(T).

The numerator is an integer. On a relatively even dual orbit the coefficients differ only by its prescribed signs, so reading the coefficient of its chosen representative gives its signed orbit-sum coordinate. On a nonsurviving orbit the coefficient is zero by invariance. Therefore all coordinates of M^(-1)f lie in 2^(-m) Z, and N=2^m M^(-1) is an integer matrix. This proves nonsingularity and both integer identities

    MN = NM = 2^m I.

Consequently the cokernel is annihilated by 2^m. Its finite order is |det M| and has no odd prime factor, proving the determinant statement. QED.

If D is empty, r=0 and the empty determinant is 1; the assertions hold. When m=0 the matrices are unimodular. The theorem does not require division by |G|.

## 2. Exact determinant and stabilizer budget

Let D_A and D_O denote diagonal matrices with entries |A| and |O|. Partial Walsh orthogonality says

    M^t D_A M = 2^m D_O.

Taking determinants and substituting orbit-stabilizer gives

    (det M)^2
      = 2^(m r) product_O |O| / product_A |A|
      = 2^(m r) product_A |Stab(A)| / product_O |Stab(O)|.

All source stabilizers are subgroups of tournament automorphism groups and therefore have odd order. Hence

    2k = m r − sum_O v_2(|Stab(O)|).

In particular the right side is a nonnegative even integer. The total 2-adic stabilizer contribution on the relatively even graph side is bounded above by m r and has the same parity as m r. This is a global identity, not a matching of the stabilizers of individual objects.

Taking odd parts in the determinant identity yields

    oddpart(product_A |Stab(A)|)
       = oddpart(product_O |Stab(O)|).

The products agree at every odd prime even when the multisets of stabilizer groups differ drastically. There is no claimed classwise stabilizer correspondence.

The weighted formula for the inverse is

    N[O,A] = |A| M[A,O] / |O|.

This expression may visibly contain odd factors in its denominator. The labelled inverse-transform proof above establishes that they cancel **entry by entry**. The cancellation must not be assumed merely from equality of the global products.

## 3. Odd primes are uniformly unexceptional here

For every odd prime p, reduction of M modulo p remains invertible. Thus the relatively even signed orbit sums form an invariant-function basis in characteristic p even if p divides |G|. No Reynolds averaging by |G| is used. The exceptional characteristic for this construction is 2.

This provides a precise connection to the repository's oddness and denominator work: the obstruction to moving between these two integer coordinate lattices is entirely 2-primary. It does **not** single out 7, 21, Monster primes, or any particular odd-prime valuation behavior in a different operator. To transfer such a claim one would need an explicit compatible operator, not merely the shared appearance of odd primes.

Nor does the cokernel measure entropy or a number of bits required to store an object. It measures the index of one integer lattice of invariant functions in another. Its order, its exponent and the number of its generators are different quantities.

## 4. Exact paired-tournament computations

Use D=Graph_n × {0,1}^n as in the preceding note. The checker constructs complete source and target orbit lists, builds every signed Fourier column through n=3, computes its determinant using exact Bareiss elimination, verifies both MN and NM, and checks the entrywise inverse divisibility. Through n=4 it also compares the prime factorizations of the products of all source and target stabilizer orders.

Results:

- n=0: r=1, m=0, |det|=1, cokernel exponent 1.
- n=1: r=2, m=0, |det|=1, cokernel exponent 1.
- n=2: r=8, m=1, |det|=2^2, cokernel exponent 2. Source stabilizer product is 1; target product is 2^4.
- n=3: r=88, m=3, |det|=2^116, cokernel exponent 8. Source stabilizer product is 3^4; target product is 2^32 3^4.
- n=4: r=2752, m=6. Source stabilizer product is 3^32; target product is 2^632 3^32. The theorem therefore gives |det|=2^7940. The large n=4 matrix was **not** constructed or its determinant independently eliminated.

The n=3 example separates the two costs sharply: lattice index 2^116, but every missing integer coordinate is repaired by multiplying by at most 8. Nothing requires a division by 3 despite both symmetry products containing 3^4.

Run from this directory:

    python3 check_lattice.py
    python3 -O check_lattice.py

Expected: 249,754 explicit runtime checks in each mode; PASS. The checker imports only the explicit permutation setup from `check_decorated_walsh.py`; its integral matrix, determinant and divisibility tests are new, but it is not an independent implementation of that permutation setup. The session's separate audit should be used for independent action conventions.

## 5. What changes for a general affine binary action

The annihilator theorem and weighted determinant identity remain valid for any finite affine binary G-space of dimension m with finite decorations, using its signed dual action. Source stabilizers need not then have odd order. The correct valuation formula is

    2k = m r + sum_A v_2(|Stab(A)|)
                  − sum_O v_2(|Stab(O)|).

Dropping the source term outside the tournament setting would be an error. Odd-part product equality still follows because M is invertible over Z[1/2].
