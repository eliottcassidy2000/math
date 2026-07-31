---
id: THM-2766
title: "Quadratic-cubic pullback, even-sign Kummer plane, and the Weyl D3=S4 coincidence"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Let R have three
  distinct nonzero roots t_i and let H(V)=R(V^2).
  Its ambient signed-root symmetry is B3=C2^3 semidirect S3.  If t1*t2*t3
  is a square in the base field, the six-root discriminant is a square and
  monodromy lies in the even-sign subgroup W(D3)=V4 semidirect S3, which
  acts on four even sign vectors as S4.  If the Kummer squareclass rank over
  the cubic splitting field is two, then a cyclic/full cubic base gives
  A4/S4 exactly.  THM-2758's quartic pair-sum sextic has this product
  constraint, and THM-2762 supplies its nonzero Keller branch.  This is an
  exact binary-over-ternary monodromy theorem, not PSL2(Z), an affine-cover
  extension, JC(2), DC(2), or LRC(14).
source: root/even-sign-binary-ternary-pullback-2026-07-28
depends_on:
  - THM-2758-quartic-pair-sum-sextic-resolvent-pullback-and-discriminant-square
  - THM-2762-quartic-opposite-sum-wall-imprimitive-d4-and-keller-exclusion
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
  - THM-2743-c2-c3-off-diagonal-projector-rank-and-s3-s4-compatibility-defect
  - THM-2756-opposite-edge-projectors-parity-cancellation-and-integral-clutch
script: 04-computation/quadratic_cubic_even_sign_weyl_d3_thm2766.py
output: 05-knowledge/results/quadratic_cubic_even_sign_weyl_d3_thm2766.out
script_sha256: a1c3f31fd0eed8ed154fbe79b6aa7992e5f61292824d2f3728891b94bf307c39
output_sha256: d4085141bfb8fbf72b6e6209e9be2659cfac290d279c8db7f7fd481954de0714
hash_basis: LF-normalized bytes
---

# THM-2766 -- the binary pullback over a ternary resolvent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The recurring primes two and three have an exact common carrier here.  A
quadratic pullback of a cubic has three binary root pairs.  Its unconstrained
symmetry is the rank-three signed permutation group.  A single product
orientation removes the odd binary direction and leaves the Klein four
plane under the ternary permutation action.  The resulting order-24 group is
simultaneously the Weyl group `W(D3)` and `S4=W(A3)`.

This is the group hidden in the pair-sum identity of THM-2758.  It is not the
free product `C2*C3`: the product constraint is precisely the relation which
turns the free binary choices into an even-sign plane.

## 1. A quadratic pullback of a cubic

Let `K` be a field of characteristic different from two, and let

```text
R(Y)=product_(i=1)^3(Y-t_i)                              (1)
```

have three distinct nonzero roots in a splitting field `E`.  Put

```text
H(V)=R(V^2)=product_i(V^2-t_i),                          (2)
Q=Gal(E/K)<=S3.
```

Choose `u_i` with `u_i^2=t_i`, and let

```text
L=E(u_1,u_2,u_3).                                        (3)
```

Then `L` is the splitting field of `H`.  Every automorphism permutes the
three unordered pairs `{+u_i,-u_i}` and may flip signs inside them, so

```text
Gal(L/K) <= B3:=C2^3 semidirect S3.                      (4)
```

The projection to `S3` has image `Q`.  Its kernel is the Kummer sign group

```text
N=Gal(L/E) <= C2^3.                                      (5)
```

If `W` is the span of `[t_1],[t_2],[t_3]` in `E*/E*2`, then Kummer theory
gives

```text
dim_F2 W=r,                       |N|=2^r,               (6)
```

and `N` is the annihilator of every linear relation among those three
squareclasses.

## 2. The product square is exactly the even-sign cut

Assume

```text
t_1t_2t_3=w^2                    for some w in K*.        (7)
```

Choose the signs of the `u_i` so that `u_1u_2u_3=w`.  Because `w` is fixed
by every `K`-automorphism, every signed permutation in `(4)` flips an even
number of the three signs.  Hence

```text
Gal(L/K) <= W(D3)
          :={(epsilon,pi) in C2^3 semidirect S3:
             epsilon_1+epsilon_2+epsilon_3=0}.            (8)
```

The same conclusion is visible from the discriminant.  The three within-pair
differences contribute `4^3 product_i t_i`; the four cross differences for
each pair `{i,j}` contribute `(t_i-t_j)^4`.  Therefore

```text
disc(H)=4^3(t_1t_2t_3) disc(R)^2.                        (9)
```

Under `(7)` this is a square.  On the six roots, a block permutation occurs
twice and is even, while one sign flip is one transposition.  Thus the
six-root sign character is exactly sign-flip parity, and `(8)` is the
permutation-theoretic content of `(9)`.

The square must lie in `K`, not merely in `E`: only then is the chosen product
`u_1u_2u_3` fixed by the whole group.  Likewise `(7)` alone does not force
rank two.  For example, `(t_1,t_2,t_3)=(1,2,8)` has square product but Kummer
rank one, while `(1,2,3)` has rank two but nonsquare product and permits odd
sign flips.

## 3. `W(D3)` is `S4`, and the cyclic subcase is `A4`

Let

```text
Omega={delta in {+1,-1}^3:delta_1delta_2delta_3=+1}.     (10)
```

It has four elements.  The even sign plane acts on `Omega` by coordinatewise
multiplication, and `S3` permutes coordinates.  This gives a faithful action

```text
W(D3)=V4 semidirect S3 -> Sym(Omega).                    (11)
```

Both groups have order 24, so

```text
W(D3) isomorphic_to S4=W(A3).                            (12)
```

Now suppose the Kummer rank in `(6)` is exactly two.  Relation `(7)` makes
`N` a subgroup of the even sign plane, and its order is four, so

```text
N=V4.                                                    (13)
```

If `Q=S3`, then `Gal(L/K)` has order `4*6=24` inside the order-24 group
`W(D3)` and is all of `S4`.  If `Q=C3`, it is the full inverse image of the
coordinate three-cycle, namely

```text
V4 semidirect C3=A4.                                     (14)
```

Thus the exact dictionary is

```text
rank-two binary Kummer plane + cyclic ternary base = A4;
rank-two binary Kummer plane + full S3 base       = S4.  (15)
```

This is a semidirect product constrained by `(7)`, not a faithful action of
the modular free product `C2*C3`.  Static binary opposition and ternary block
rotation may commute on a quotient even when the full group is `A4/S4`.

## 4. Quartic pair-sum specialization

For the quartic notation of THM-2758, let

```text
d_1=s_12-s_34,       d_2=s_13-s_24,       d_3=s_14-s_23,
t_i=d_i^2/4,         T=d_1d_2d_3.                         (16)
```

Its translated pair-sum sextic is exactly

```text
G(e1/2+V)=R(V^2),                                      (17)
```

and

```text
t_1t_2t_3=(T/8)^2.                                     (18)
```

For quartic monodromy `A4` or `S4`, the matching quotient has group `C3` or
`S3`, and adjoining the three signed differences recovers the `V4` kernel;
therefore the Kummer rank is two.  Sections 1--3 recover `A4/S4` from the
binary-over-ternary pullback itself.

For a generic degree-four Keller extension, THM-2762 proves `T!=0`, so the
nonzero/separable hypotheses above hold for every fixed primitive-element
presentation.  What remains open is geometric: neither `(17)` nor its Weyl
group supplies an affine integral model across the Jelonek boundary or the
unit/class-group Kummer realization required by THM-2655.

## 5. Exact verification and boundary ledger

Run

```bash
python 04-computation/quadratic_cubic_even_sign_weyl_d3_thm2766.py
python -O 04-computation/quadratic_cubic_even_sign_weyl_d3_thm2766.py
```

The exact companion uses explicit exceptions and no truth-bearing Python
assertions.  It enumerates all `48` signed permutations of three root pairs,
checks the six-root parity character, identifies the `24` even-sign elements,
proves their four-point action is all of `S4`, and identifies the `12`-element
cyclic-base preimage with `A4`.  It exhausts all `16` squareclass-relation
subspaces of `F2^3`, verifies the rank-two/product relation has annihilator
`V4`, checks `(9)` on `120` distinct nonzero integer triples, and verifies the
quartic specialization on all `210` distinct root quadruples in `[-4,5]`.

```text
PROVED HERE:              signed-wreath ambient group for R(V^2);
                          Kummer rank/sign-kernel dictionary;
                          product-square even-sign cut;
                          exact discriminant formula;
                          W(D3)=V4 semidirect S3=S4;
                          cyclic-base preimage A4;
                          quartic pair-sum specialization.

NOT PROVED:               a free-product PSL2(Z) action;
                          that square product alone forces rank two;
                          affine/integral extension of the sextic;
                          the THM-2655 unit/class-group carrier;
                          exclusion of A4/S4;
                          JC(2), DC(2), or LRC(14).                       (19)
```

QED.
