# A positive profile algebra and a nine-state obstruction

**Status.** Exact identities and a finite induction certificate for a specified construction grammar. This is not an obstruction for arbitrary tournaments, and no novelty claim against the literature is made.

## Sources and scope

The input is the full path-cover calculus in [S159's positive profile operator](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session159_constructions_20261103/profile_operator/NOTE.md). We retain its notation: p_k(X) counts unordered covers of V(X) by k nonempty directed paths, P_X(z)=sum p_k(X) z^k, and

    F_X(r) = sum_k p_k(X) (r)_k.

Here (r)_k is a falling factorial; empty bins are allowed in F. The construction grammar C starts with a singleton, allows ordinal sums X⊕Y, and allows triangle wrappers G(X)=X→a→b→X. All objects in C are nonempty. A root G is strong. We do not assume a root ordinal sum is strong.

## 1. The matching product

Define on polynomials with nonnegative integer coefficients

    P ⋆ Q = sum_(t>=0) z^t/t! · P^(t)(z) Q^(t)(z).                 (1)

Despite the displayed denominator, every coefficient is an integer: the contribution from z^k and z^l is

    sum_(t=0)^min(k,l) binom(k,t) binom(l,t) t! z^(k+l−t).        (2)

Then

    P_(X⊕Y) = P_X ⋆ P_Y.                                        (3)

**Proof.** Every path in an ordinal sum meets X first and Y second. Restrict a cover to the two blocks. If the restrictions have k and l paths, exactly t pairs of an X-path and a Y-path have been concatenated. Select their endpoints as whole paths and match them in binom(k,t)binom(l,t)t! ways. All joining arcs point forward. Conversely every such matching gives a unique cover with k+l−t paths. ∎

This makes the polynomial carrier a commutative associative semiring with ordinary addition, ⋆, and unit 1. One proof of associativity is the injective basis change

    Φ(z^k)=(r)_k,       Φ(P ⋆ Q)=Φ(P)Φ(Q).

Ordinary addition represents a formal alternative between cover classes; it is **not** asserted to be a tournament operation. Nor does commutativity mean X⊕Y and Y⊕X are isomorphic.

A useful formal-series identity is

    exp(az) ⋆ exp(bz) = exp((a+b+ab)z).                           (4)

Thus the algebra naturally carries the multiplicative formal-group law a+b+ab. This is a precise algebraic bridge, not a claim that a new tournament operation has been found. Multiplication by z in this algebra is the differential operator z(1+D), where D=d/dz.

S159's triangle operator is

    L = (z+z²)D² + (3z+2z²)D + (z+z²),    P_G(X)=L(P_X).       (5)

Together (1) and (5) give a positive differential calculus for all construction trees in C. The matching proof also works for marked internal arcs with coefficients in N[t_1,…,t_s], when all cross-block arcs have weight 1. S159 supplies the corresponding coefficientwise triangle identity.

## 2. Two exact finite quotients

For any integer M>=2, the vector

    f_M(X)=(F_X(0),…,F_X(M−1)) mod M

is a closed finite state. Since F_X(r) is a polynomial over Z, its values mod M are periodic in r with period M. Ordinal sum is coordinatewise multiplication, and triangle substitution is

    (L_M f)_r = r f_(r+1) + r² f_r − r f_(r−1) mod M,           (6)

with cyclic subscripts mod M. Formula (6) follows from S159's polynomial identity. For nonempty sources the zeroth coordinate is zero. This gives a deterministic finite bottom-up tree automaton for any fixed modular question about C. Saturation from the singleton terminates and is exact for this grammar, though its worst-case state space can be large.

There is also a coefficient quotient. Truncate p_0,…,p_N modulo M. Formula (2) shows that an output coefficient of degree at most N never uses an input degree above N. The only omitted input to (5) is N(N+1)p_(N+1), so this truncation is closed whenever

    M divides N(N+1).                                          (7)

The standard choice N=M−1 works, but need not be the smallest positive choice: for M=6, N=2 suffices. For a composite M, the smallest positive N can be found by assigning each prime-power factor of M to either N or N+1 and using CRT. This quotient is not claimed to determine the entire integral profile or tournament realizability.

## 3. Exact mod-12 obstruction for the entire grammar

**Theorem.** Every X in C has

    H(X) mod 12 in {1,3,5,9}.                                  (8)

In particular C contains no source with H congruent to 7 or 11 mod 12. The statement applies to every construction depth and every parenthesization, not only to one iterated family.

Here is a finite certificate. The nine possible f_12 states are listed below as evaluation vectors at r=0,…,11. All nine are reachable. In the witness column v is a singleton and multiplication denotes ordinal sum at the profile level.

| State | f_12 | Witness | L_12 state |
|---|---|---|---|
| 0 | (0,1,2,3,4,5,6,7,8,9,10,11) | v | 1 |
| 1 | (0,3,0,9,0,3,0,9,0,3,0,9) | G(v) | 1 |
| 2 | (0,1,4,9,4,1,0,1,4,9,4,1) | v*v | 6 |
| 3 | (0,3,0,3,0,3,0,3,0,3,0,3) | v*G(v) | 3 |
| 4 | (0,9,0,9,0,9,0,9,0,9,0,9) | G(v)*G(v) | 4 |
| 5 | (0,9,0,3,0,9,0,3,0,9,0,3) | v*G(v)*G(v) | 5 |
| 6 | (0,5,8,9,8,5,0,5,8,9,8,5) | G(v*v) | 2 |
| 7 | (0,1,8,3,4,5,0,7,8,9,4,11) | v*v*v | 5 |
| 8 | (0,5,4,3,8,1,0,11,4,9,8,7) | v*G(v*v) | 5 |

The complete coordinatewise product table, with rows and columns ordered 0,…,8, is

    2 3 7 1 5 4 8 2 6
    3 4 1 5 1 3 1 3 3
    7 1 2 3 4 5 6 7 8
    1 5 3 4 3 1 3 1 1
    5 1 4 3 4 5 4 5 5
    4 3 5 1 5 4 5 4 4
    8 1 6 3 4 5 2 8 7
    2 3 7 1 5 4 8 2 6
    6 3 8 1 5 4 7 6 2

The table contains the leaf and is closed under both operations, proving inclusion by structural induction. The witnesses prove equality. Reading coordinate 1 gives (8). Strong roots already attain all four residues: G(v) has residue 3; G(v*v) has residue 5; G(G(v*v)) has residue 1; G(v*v*v) has residue 9.

The projection onto mod 3 hits every residue, and the projection onto mod 4 hits both odd residues. These separate successful projections therefore miss a genuine joint mod-12 obstruction. Modular universality in individually selected moduli is not a substitute for checking their coupled state.

## 4. Why the guard interface matters

[GUARD_ALGEBRA.md](GUARD_ALGEBRA.md) gives a different literal operation, obtained by allowing both choices of the next moving guard in S163. Starting from a directed triangle, two b-guard steps give a strong seven-vertex tournament with H=31, which is 7 mod 12. It therefore lies outside C. The point is not that 31 was unknown, but that the extra physical endpoint interface provably escapes this entire profile construction grammar.

## Reproduction

Run `python3 check.py` and `python3 -O check.py`. Independent subset/endpoint DP verifies the triangle formula on all 75 labeled tournaments of orders 1 through 4 and the matching product on all 121 ordered pairs of tournaments of orders 1 through 3. The code separately checks every entry of the finite mod-12 closure certificate. Its conditions remain active under optimized Python. The all-depth obstruction is structural induction on that exact finite table.

## Endpoint module extension

[POINTED_PROFILE_ALGEBRA.md](POINTED_PROFILE_ALGEBRA.md) combines the matching product with moving guards in a single positive typed algebra. It proves full endpoint-refined cover transfers and a positive lift from an ordinary profile to the guarded triangle. The same cutoff condition M|K(K+1) gives a finite quotient of this larger grammar. The mod-12 obstruction above remains scoped to the smaller grammar C.
