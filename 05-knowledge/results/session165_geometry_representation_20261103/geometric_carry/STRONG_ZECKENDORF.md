# Strong Zeckendorf codes and their exact 2-adic geometry

**Status:** proved synthesis and finite independent controls. This upgrades the S158 generally nonstrong Zeckendorf product code using S164's already proved physical guard returns. The return operations and the elementary 2-adic decoding are inherited; the strong implementation, exact distance formula, and constrained limit-set description are the conclusions here. No literature-priority or full-spectrum claim is made.

## 1. A physical, source-closed code

A guarded tournament has a vertex q with unique in-neighbor r. Its endpoint state is `(s,u,d)`, where s counts Hamilton paths beginning at q, u counts all remaining Hamilton paths, and d counts Hamilton paths of the deletion X−q beginning at r. Start with a directed triangle and any guard, with state `v0=(1,2,0)`.

A one-vertex attachment adds a with `q→a` and `a→X−q`. Marking a as new guard is C; retaining q is D. S164 proves, by a path decomposition,

    C(s,u,d)=(u,s+u+d,s),
    D(s,u,d)=(s+d,2u,0).

Each operation preserves strong connectivity: the old tournament is strong and the new directed cycle `r→q→a→r` joins a to it. Each operation costs one vertex. Words are read temporally from left to right.

For j≥1 put `R_j=D^(j−1) C D`. Direct substitution gives

    R_j[t v0]=(1+2^j)t v0,           cost j+1 vertices.       (1)

Thus these are reusable physical returns, not numerical multipliers without an output interface.

Let `F_2=1,F_3=2` and use the usual Fibonacci recurrence. For N≥0 write its unique Zeckendorf expansion as `N=Σ_(i∈I) F_i`, where `I⊂{2,3,...}` is finite and has no adjacent indices. Apply the blocks `R_(i+2)` in increasing i. The resulting guarded tournament Z_N is strong and canonical, with

    H(Z_N)=3 ∏_(i∈I)(1+2^(i+2)),
    |V(Z_N)|=3+Σ_(i∈I)(i+3).                              (2)

In particular N=0,1,2,3,4 yield `(H,order)` equal to `(3,3)`, `(51,8)`, `(99,9)`, `(195,10)`, `(3315,15)`. The vertex count is `O((log(N+2))²)`: the largest index is `O(log(N+2))`, and the sum of all indices through it is quadratic. An ordering of the blocks is fixed to specify actual graphs; commuting scalar multipliers does not establish graph isomorphism between differently ordered blocks.

For any nonempty product of distinct factors `(1+2^j)`, its smallest exponent is the 2-adic valuation of the product minus one. Divide by that factor and repeat. Applying this to `H(Z_N)/3` recovers I and then N. Hence (2) is injective. The source states close exactly after every block.

This is the same digit-factor code as S158, with its scalar core 265 replaced by the canonical triangle factor 3 and with strong, reusable source realizations. It does **not** preserve S158's square-sum vertex labels or reversed-interval geometry. It does not supply affine `3h+2` or `3h+4` returns: every count in (2) is divisible by 3.

## 2. Exact prefix separation

Write `H_I=3∏_(i∈I)(1+2^(i+2))`, initially for arbitrary finite sets, without the no-adjacency restriction. If I≠J then

    v_2(H_I−H_J)=2+min(I symmetric-difference J).            (3)

**Proof.** Cancel the common factors, which are odd. Put k=min(I△J). In the two remaining products, exactly one contains the factor `(1+2^(k+2))`; every other factor is 1 modulo `2^(k+3)`. Their difference is therefore `±2^(k+2)` modulo `2^(k+3)`. Multiplication by 3 and the common odd factors does not alter the valuation. ∎

Thus low 2-adic precision recovers exactly the initial digit block. There is no collision hidden by large later digits. At modulus `2^K`, K≥4, only indices `2≤i≤K−3` are visible. There are exactly

    F_(K−2)                                                (4)

residues attained by the Zeckendorf family: there are K−4 available positions, and binary strings of length n without adjacent ones number `F_(n+2)`. Every such prefix occurs as a finite Zeckendorf support, and (3) distinguishes them. The occupied fraction of all `2^K` residues therefore tends to zero.

## 3. Infinite addresses and a 2-adic Cantor set

For an arbitrary infinite support I⊂{2,3,...}, the product in (2) converges in Z_2, because the factors tend to 1. Formula (3) continues to hold: only finitely many factors affect any specified residue. It gives a homeomorphism from binary supports, with distance `2^(−2−min(I△J))`, onto the whole ball `3+16 Z_2`. To prove onto, solve one binary digit at a time. At each new exponent j≥4, multiplying by `(1+2^j)` changes exactly the next unresolved digit modulo `2^(j+1)`. This also proves uniqueness.

Restrict to no-adjacency supports. Their image is a compact, nowhere-dense set K. Finite strong codes are dense in K, and K has Hausdorff dimension

    dim_H K = log_2(phi),          phi=(1+sqrt(5))/2.         (5)

For an elementary justification, after n visible positions there are `F_(n+2)` cylinders of diameter at most a fixed constant times `2^(−n)`. This gives the upper bound. For the lower bound use the adjacency matrix `[[1,1],[1,0]]` on the last digit. Its positive right eigenvector is `(phi,1)`, with eigenvalue phi. Give an allowed transition a→b probability `A_ab r_b/(phi r_a)`, starting in state 0. Every n-digit cylinder then has probability `r_last/(phi^n r_0)`, between two positive fixed multiples of `phi^(−n)`. A 2-adic ball of sufficiently small radius fixes an initial digit block by (3), so its measure is bounded by a constant times `radius^(log_2 phi)`. The standard covering inequality (sum of ball diameters to that power dominates the total mass up to a constant) proves the lower bound. No numerical dimension estimate is used as evidence.

The full-support ball and K consist of 2-adic limits. An infinite support is **not** a finite tournament, and no Hamilton-path count is assigned to an infinite graph. Only finite supports in (2) have the literal graph interpretation.

There is a useful separation of geometries. The finite code values are infinite and hence Zariski dense in the one-dimensional multiplicative group over Q: a nonzero one-variable polynomial has only finitely many roots. Nevertheless their no-adjacency limit is nowhere dense 2-adically, has dimension strictly less than one, and occupies only the Fibonacci number of residue classes in (4). Zariski density alone cannot certify modular universality or exact odd-count coverage.

## Sources and checks

- [S158 Zeckendorf square-sum code](https://github.com/eliottcassidy2000/math/blob/89416555c4d95c6fa20c9ebafca1a180ad3fcaf3/05-knowledge/results/session158_creative_20261103/square_overlap/ZECKENDORF.md): the digit-factor choice, injective valuation decoding, and generally nonstrong original implementation.
- [S164 guard algebra, Sections 5–6](https://github.com/eliottcassidy2000/math/blob/89416555c4d95c6fa20c9ebafca1a180ad3fcaf3/05-knowledge/results/session164_emergent_20261103/profile_algebras/GUARD_ALGEBRA.md): literal one-vertex transfers, strongness, and the already proved return (1).

Run `python3 check.py` and `python3 -O check.py` for independent controls. The script builds actual tournament adjacency matrices, independently counts Hamilton paths and endpoint/deletion states on selected small codes, verifies strongness, tests decoding and exact 2-adic separation, and counts every allowed residue prefix in a specified finite range. The all-size statements rest on the proofs above.

**Verification receipt.** Both normal and optimized Python runs pass all 729,799 active checks and produce byte-identical JSON. The code controls decode N=0,…,9999, test all 523,776 unordered pairs of ten-position supports, count residue prefixes through modulus 2^19, and independently count eight actual strong tournaments through order 15.
