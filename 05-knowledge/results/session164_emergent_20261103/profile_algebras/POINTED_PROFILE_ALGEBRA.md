# A positive endpoint module closes profiles and moving guards

**Status.** Exact general identities, with independent finite controls. This unifies the S159 full-profile operator and the S163 moving-guard mechanism, including the one-vertex factorization developed in this session. It does not prove exact odd-spectrum coverage or a closed ternary observer bank. No literature-priority claim is made.

Read [PROFILE_ALGEBRA.md](PROFILE_ALGEBRA.md) for the matching product ⋆ and [GUARD_ALGEBRA.md](GUARD_ALGEBRA.md) for the scalar endpoint calculus and explicit S159 graph isomorphism. The source identities are S159's positive profile operator and S163's moving-guard/first-step notes, linked there.

## 1. The full pointed profile

Let q have a unique in-neighbor rho in a tournament X. For each integer r>=1, define three nonnegative counts:

- S_r: ordered r-slot path covers of X in which q starts its path;
- U_r: ordered r-slot path covers of X in which q does not start its path;
- d_r: ordered r-slot path covers of X−q in which rho starts its path.

Empty slots are allowed, and F_X(r)=S_r+U_r. Starting a path means being its first vertex, including a singleton path. Every count is a polynomial in r with an expansion in falling factorials and nonnegative integer coefficients: its coefficient of (r)_k counts the corresponding unordered k-path covers.

Write their unordered-cover polynomials as S(z), U(z), and d(z). All have zero constant term, because X has at least two vertices and X−q is nonempty. We use the same letters for sequences only when their argument or subscript makes the meaning clear.

Add a vertex a with q→a and a→X−q. In move C mark a as the new guard; in move D retain q. Both are valid guards. Then the full ordered-slot identities are

    C_r(S_r,U_r,d_r)
      = ((r−1)S_r+rU_r, S_r+U_r+d_r, S_r),                  (1)

    D_r(S_r,U_r,d_r)
      = (rS_r+d_r, (r+1)U_r, (r−1)d_r).                    (2)

At r=1 these specialize exactly to the Hamilton-path C/D calculus.

### Proof

For a fixed old r-slot cover, inserting a as the initial vertex of a slot gives r possible slots, except that the slot beginning at q is forbidden. Empty slots are included: in such a slot a is a singleton. There is one additional insertion immediately after q, valid whether q ends its old path or has an old successor. These account for (r+1)F_X(r)−S_r output covers.

The only outputs whose deletion of a does not leave an old cover have the segment q,a,rho. Since rho is the only old in-neighbor of q and occurs after it, q must begin that slot. Deleting q,a leaves exactly a cover of X−q in which rho begins its slot, giving d_r additional covers. Thus

    F_new(r)=(r+1)F_X(r)−S_r+d_r.                            (3)

For C, the covers in which a starts its slot are the first insertion class, namely rF_X(r)−S_r. Subtracting from (3) gives the new U_r=F_X(r)+d_r. Deleting a restores X and the new guard's unique in-neighbor is q, so its new defect is S_r. This proves (1).

For D, an old cover with q starting its slot has r−1 allowed insertions at the starts of other slots and one insertion immediately after q, hence r insertions. The exceptional q,a,rho covers contribute d_r. Therefore its new S_r is rS_r+d_r; subtracting from (3) gives (r+1)U_r. After deleting q, a is a source. To keep rho at the beginning of its path, a must be placed at the start of one of the other r−1 slots. Thus the new defect is (r−1)d_r, proving (2). ∎

## 2. Positive operators on unordered covers

Let ∂=d/dz and let

    N=z(1+∂).

Under the falling-factorial basis change Φ(z^k)=(r)_k, N is multiplication by r. Hence (1) and (2) become

    C(S,U,d)=((N−1)S+NU, S+U+d, S),
    D(S,U,d)=(NS+d, (N+1)U, (N−1)d).                         (4)

These are coefficientwise nonnegative on the stated carrier. In particular N−1 is nonnegative on polynomials with zero constant term:

    (N−1) z^k=(k−1)z^k+z^(k+1),       k>=1.

More explicitly, writing s_k,u_k,d_k for coefficients and setting negative-index terms to zero, C sends

    s'_k=(k−1)s_k+s_(k−1)+ku_k+u_(k−1),
    u'_k=s_k+u_k+d_k,
    d'_k=s_k,

while D sends

    s'_k=ks_k+s_(k−1)+d_k,
    u'_k=(k+1)u_k+u_(k−1),
    d'_k=(k−1)d_k+d_(k−1).                                 (5)

No coefficient at a lower degree depends on a higher-degree input. Consequently these two pointed operators admit exact coefficient truncation at **any** degree, even before reducing modulo an integer.

## 3. The module action and the triangle lift

For any tournament Y, the ordinal graft X⊕Y preserves the unique in-neighbor of q. Every slot splits uniquely at the block boundary. Whether q starts its slot is decided wholly inside X, and the analogous statement holds for rho in X−q. Therefore

    (S,U,d)_(X⊕Y,q)=(S_X ⋆ P_Y, U_X ⋆ P_Y, d_X ⋆ P_Y).       (6)

This is a positive module action over the profile semiring. It is a right ordinal graft: putting a nonempty block before X would generally destroy the unique-in-neighbor condition. An ordinal graft of two nonempty blocks is not strong, and no preservation of strength is claimed for this operation.

Conversely the triangle wrapper gives a canonical way to lift an ordinary profile to a pointed one. For nonempty X let G(X)=X→a→b→X and mark b. Its unique in-neighbor is a. In ordered slots,

    d_r = r F_X(r−1),
    U_r = r F_X(r+1),
    S_r = r²F_X(r)−rF_X(r−1).                               (7)

To prove the defect formula, remove b. The remaining graph is X⊕a, and a begins its path only when it is a singleton: choose its slot in r ways and cover X in the remaining r−1 slots. To count covers where b starts, first cover X⊕a, then prepend b in any slot except a slot beginning at a. This gives r F_(X⊕a)(r)−rF_X(r−1)=r²F_X(r)−rF_X(r−1). Subtracting from S159's full triangle identity gives U_r in (7).

All subtraction disappears in unordered-cover coordinates. If P=P_X, then

    d=zP,
    U=z(1+∂)²P,
    S=z²P+(z+2z²)P'+z²P''.                                  (8)

For example P=z gives

    S=z+2z²+z³,       U=2z+z²,       d=z²,

and S+U=3z+3z²+z³ is the directed triangle's path-cover polynomial. The Hamilton endpoint state is (S_1,U_1,d_1)=(1,2,0).

Equations (4), (6), and (8), together with forgetting the guard by P=S+U, form one positive algebra for ordinary profiles, guarded profiles, ordinal sums, triangle resets, and moving guards. The realizable objects are the ones obtained by literal operations; arbitrary nonnegative coefficient vectors are not claimed to be realizable by tournaments.

## 4. Uniform finite quotients and a recurrence at every slot count

Fix a modulus M>=2 and a positive cutoff K with M dividing K(K+1). Retain coefficients through degree K modulo M, for either an ordinary profile or a pointed triple. This is a closed finite typed algebra:

- ⋆ and the module action never bring higher degrees down;
- C,D and forgetting the guard never bring higher degrees down;
- in the triangle lift (8), the sole omitted contribution is K(K+1)p_(K+1) to the degree-K term of U, which vanishes modulo M.

Thus all these operations can be combined in one finite exact modular reachability problem. This supplies a decision procedure for the specified finite-seed grammar and modular output condition, not for arbitrary tournament realizability or exact integer coverage. It is particularly useful for detecting a coupled modular obstruction before pursuing an infinite source bank.

At fixed r the matrix C_r in (1) is

    [r−1 r 0]
    [ 1  1 1]
    [ 1  0 0],

with characteristic polynomial

    λ³−rλ²−λ−r.                                             (9)

Therefore every coordinate, and F=S+U, along repeated C satisfies the same order-three recurrence

    F_(n+3)(r)=rF_(n+2)(r)+F_(n+1)(r)+rF_n(r).              (10)

This holds for every valid guarded seed. At r=1 it is tribonacci; at general r it describes the **entire cover profile**, not only the scalar Hamilton count. Since it holds at every positive integer r, it is also a polynomial identity in the slot variable.

### Full-profile generating function for the S159 interval family

Let T_m be S159's tournament on 0,…,2m+1 with backward edges (2i+3)→(2i), i=0,…,m−1. The literal isomorphism in GUARD_ALGEBRA.md identifies T_m with every second C-step from a directed edge whose sink is the guard. That edge has pointed slot state

    (S_r,U_r,d_r)=(r²−r,r,r).

The matrix C_r² has characteristic polynomial

    μ³−(r²+2)μ²+(1−2r²)μ−r².

Its first profile values are

    F_(T_0)(r)=r²,
    F_(T_1)(r)=r⁴+4r²,
    F_(T_2)(r)=r⁶+8r⁴+8r².

Consequently the entire interval-family path-cover profile has the rational generating function

    sum_(m>=0) F_(T_m)(r) z^m
      = r²(1+z)² / [1−(r²+2)z+(1−2r²)z²−r²z³].             (11)

At r=1 this recovers S159's Hamilton-path generating function exactly. Equation (11) is its full-profile refinement. It follows from the proved graph isomorphism and transfer, not from fitting finite values.

## 5. Useful next targets

1. Saturate the enlarged pointed grammar in selected composite moduli and look for the next obstruction after the old mod-12 obstruction disappears. Preserve the typed state and do not infer exact coverage from modular coverage.
2. Extend the fixed observer queries to three pointed triples and search for a genuine all-branch return with the required affine observer correction. This is the physical-interface target; a scalar residue witness is weaker.
3. Seek an odd-modulus reachability theorem, or a counterexample modulus, for the small C/D seed bank. Existing finite controls and the proved power-of-two theorem do not settle this.

## Reproduction

Run `python3 check.py` and `python3 -O check.py`. The checker independently counts unordered covers and endpoint-restricted covers by choosing the distinguished vertex's path block and applying subset/endpoint Hamilton-path DP. It verifies (4) on all 110 qualifying guarded labeled tournaments through order 4, (6) on each such source grafted to each of the three labeled tournaments of orders 1 and 2, and (8) on all 75 labeled sources of orders 1 through 4. It also independently controls the interval-profile recurrence in (11) for m=0,…,4 and r=1,…,4. These are finite controls. The insertion bijections, module factorization, and positive lift prove the general formulas.
