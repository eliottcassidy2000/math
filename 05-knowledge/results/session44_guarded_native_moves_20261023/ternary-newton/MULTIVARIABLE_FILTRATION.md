# The exact multivariable ternary augmentation filtration

**PROVED** on the regular table representation of `(Z/3^k)[F3^r]`, for integers `r,k>=1`. The all-ring Newton compiler in THEOREM.md is the starting mechanism; this supplement identifies the full multiaxis nilpotency index, not merely a sum of separately nilpotent coordinate bounds. All original native source, shared-port and electrical guard obligations remain external.

## 1. Carrier and statement

Let `M=(Z/3^k)^G` be the module of all literal tables on `G=F3^r`. Let

    (T_i f)(a)=f(a+e_i),       D_i=T_i-I.

The `T_i` commute and have order three. Let `A` be their generated operator algebra and `J=(D_1,...,D_r)` its augmentation ideal. This is the regular representation: the map from the group algebra `(Z/3^k)[G]` into `End(M)` is faithful, because applying an operator to the delta table at zero recovers its group coefficients, up to reversing their indices.

**Sharp theorem.**

    J^(2(k+r)-1)=0,       J^(2(k+r)-2) != 0.

Thus the exact ideal nilpotency index is `2(k+r)-1`. At `r=1` this is the previously proved single-axis `2k+1`; at `k=1` it becomes the characteristic-three truncated-polynomial index `2r+1`.

The last nonzero ideal is particularly simple:

    J^(2(k+r)-2) = {0, z, 2z},
    z = 3^(k-1) sum_(g in G) T_g.

On tables, `z` sends `f` to the constant table `3^(k-1) sum_a f(a)` modulo `3^k`. The surviving final datum is exactly total table mass modulo three. It is not a claim that the whole table, its input support, or its native phase is determined by that mass.

## 2. Exact integral valuation of a monomial

First work on integer-valued tables, before reducing modulo `3^k`. Write `content_3(U)` for the minimum3-adic valuation among the nonzero entries of an integer operator matrix. The identity operator has content zero.

For a single coordinate, THEOREM.md proves

    D^(2s+1)=(-3)^s T^s D,
    D^(2s+2)=(-3)^s T^s D^2          (s>=0).

Both `D` and `D^2` have some integer entries prime to three. Since `T` is a permutation matrix, multiplication by its powers only permutes rows. Consequently

    content_3(D^j)=floor((j-1)/2)     for j>=1,
    content_3(D^0)=0.

For a multiexponent `j=(j_1,...,j_r)`, the literal monomial matrix is a tensor product of the one-axis matrices. Every entry is a product, so every entry is divisible by the sum of the coordinate valuations. Choose an entry of minimal valuation independently in each factor to attain that sum. Therefore the exact, not merely lower, valuation is

    nu(j) = sum_(i:j_i>0) floor((j_i-1)/2).

The monomial is nonzero modulo `3^k` iff `nu(j)<k`. When nonzero, its additive order as an operator is exactly `3^(k-nu(j))`.

This equality is a statement about one monomial in the faithful regular representation. In a sum of monomials, cancellation can increase the valuation; one must not replace the valuation of a general sum by the minimum of its summands without a noncancellation proof.

## 3. Sharp nilpotency index

Let `s` be the number of positive exponents in a monomial. For each positive exponent,

    j_i <= 2 floor((j_i-1)/2)+2.

If the monomial survives modulo `3^k`, then `nu(j)<=k-1`, whence

    |j| <= 2nu(j)+2s <= 2(k-1)+2r.

Every monomial of total degree `2(k+r)-1` therefore vanishes. Such monomials generate `J^(2(k+r)-1)`, proving the upper bound for the entire ideal, including arbitrary coefficient combinations.

For the lower bound, take exponents

    j_1=2k,      j_2=...=j_r=2.

Their total degree is `2(k+r)-2` and valuation is exactly `k-1`. Their operator is consequently nonzero modulo `3^k`. This proves sharpness.

## 4. Identifying the last surviving ideal, not just its degree

At the maximal surviving total degree `2(k+r)-2`, equality must hold in every inequality above. Thus each coordinate is active, each exponent is even, and

    j_i=2(v_i+1),       v_i>=0,       sum_i v_i=k-1.

There are `binom(k+r-2,r-1)` such monomials. Their number is not the dimension of the top ideal: they all agree up to the same unit scalar.

Indeed each such monomial is

    (-3)^(k-1) product_i (T_i^v_i D_i^2).

Modulo three,

    D_i^2=(T_i-I)^2 = I+T_i+T_i^2,

and this cyclic norm is invariant under multiplication by `T_i`. Since the expression already contains `3^(k-1)`, reduction of the remaining factor modulo three is exact modulo `3^k`. Hence every surviving monomial at that degree equals

    (-1)^(k-1) 3^(k-1) product_i (I+T_i+T_i^2)
      = (-1)^(k-1) z.

All other monomials of that total degree vanish. Multiplying `z` by an arbitrary group-algebra element multiplies it by the coefficient sum, because every translation fixes the full norm. Thus the whole ideal is exactly the three-element group `{0,z,2z}`. Its generator is nonzero in the regular representation and is killed by three.

## 5. Scope boundary and relationship to the actual guard problem

On every `(Z/3^k)[G]`-module the same ideal-power upper bound holds, because it is already an identity in the group algebra. Sharpness is asserted only on the regular representation or a module known to retain its top norm element. On the trivial representation every `D_i` is zero, so the nilpotency index drops immediately.

Tensoring coordinate bounds alone would give the weaker bound `2kr+1`; the exact improvement to `2k+2r-1` comes from sharing the same factor of three across coordinates. This is a genuinely coupled coefficient effect, not independent bookkeeping of coordinate nilpotency indices.

The statement applies to supplied ternary input tables. It does not reassign the original electrical guard from modulus8 to modulus9, prove invariance under native moves, identify observer equivalence with connectivity, or infer all-size metagraph wiring. If different authenticated observables use moduli8 and9, they remain distinct output components with different algebraic decompositions.

## 6. Independent exact checks

The separate audit packet contains multivariable_check.py and normal/optimized receipts. For every `1<=r,k<=4`, it enumerates every monomial at degrees `N-1,N,N+1`, with `N=2(k+r)-1`. Two independent paths construct its integral regular delta column: tensoring one-coordinate difference columns, and literally applying each coordinate's translations and differences on all `3^r` table entries.

It checks exact matrix-content valuations, every survival/vanishing decision, all boundary survivor counts, equality of every top survivor to the same signed norm column, translation-invariance of that norm, and its exact additive order three. Every matrix column is a translate of the delta column, so these are full regular-operator tests, not evaluations on an arbitrary test vector. The proof above, not the finite range, establishes the all-`r,k` statement.

    python3 -B multivariable_check.py
    python3 -O -B multivariable_check.py

All arithmetic is integral and all tests survive optimized Python.
