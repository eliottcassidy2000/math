# Complementary ears: repeated scalar offsets and two-sided composition

**Status.** This gives a constructive escape from S163's exclusion of consecutive small offsets for its *same two-vertex ear*. A different second attachment realizes two consecutive equal scalar offsets, universally under the stated guard hypothesis. A commuting two-sided construction gives the literal strong chain `5→17→53→161`, with all three steps of the form `h→3h+2`, and an infinite explicit two-parameter family. These are scalar constructions, not a closed bank for the marked observer `t=(2,−1,0)`, and not a tournament-spectrum classification. No novelty against the literature is asserted.

Source credit: the initial ear, defect, and two example sources come from [S163 FIRST_STEP](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session163_ternary_triples_20261103/strong_returns/FIRST_STEP.md). The complementary one-vertex guard calculus is proved in this session’s [guard algebra](../profile_algebras/GUARD_ALGEBRA.md); its relevant formulas are stated below so that the present constructions can be read independently.

## 1. Change the second attachment, not just the guard

Let `X` be a strong tournament, and let `q` have unique in-neighbor `r`. Write

    h=H(X),   u=H(X)−H_start(X,q),
    d=H_start(X−q,r).

A one-vertex left attachment adds `a` with `q→a` and `a→X−q`. Both the old `q` and the new `a` are again unique-in-neighbor guards. The two choices of marked guard have exact profiles

    C(h,u,d)=(h+u+d, h+d, h−u)       [mark a],
    D(h,u,d)=(h+u+d, 2u, 0)         [retain q].

Here and below words are in **temporal order**. Thus `CD` means first C, then D. Every attachment preserves strongness: the old strong graph reaches the new vertex through q, while the new vertex reaches an old vertex other than q and hence all of X.

For completeness, a path of X starting q is precisely a path of `X−q` whose first vertex is not r, so `H(X−q)=h−u+d`. Deleting the added vertex gives the one-vertex count `h+u+d`. Paths starting the new vertex number u. Paths starting the old q number `H(X−q)`. These facts give the first two coordinates above. Deleting the new guard recovers X, giving defect `h−u`; deleting the retained old guard leaves the new vertex a source, giving defect zero.

### Two equal offsets

For every integer `d≥0`, temporal `CDDC` yields

    (h,u,d) → (3h+2d, 2h+2d, 0)
            → (9h+8d, 5h+4d, h),

where the first arrow is CD and the second is DC. In particular the two scalar steps are

    h → 3h+2d → 3(3h+2d)+2d.

The variant `CDDD` has the same final count, with final profile

    (9h+8d, 8h+8d, 0).

**Proof.** Substitute the displayed C and D formulas. For the second step, starting from `(h1,u1,0)`, temporal DC gives `(h1+3u1,h1+u1,h1−u1)`, while DD gives `(h1+3u1,4u1,0)`. Taking `h1=3h+2d` and `u1=2h+2d` gives the claims. ∎

This is not S163's prohibited repetition: its next two-vertex ear is CC or CD, whereas the successful second attachment is DC or DD. Concretely, after the first triangle `q→a→b→q`, add two further vertices, each losing to a and beating every other vertex present at its insertion. The intermediate graph is the original two-vertex ear; the final graph has a transitive three-vertex wing beyond a.

### Literal +2 and +4 witnesses

Adjacency row masks use bit j of row i to mean `i→j`.

* From `X=(8,1,3,6)`, guard `q=2`, one has `(h,u,d)=(5,3,1)`. Temporal CDDC gives the strong chain `5→17→53` and final graph

      (8,1,19,6,235,15,47,111),  new guard 7,
      final guard profile (53,29,5).

* From `X=(20,1,2,7,14)`, guard `q=3`, one has `(h,u,d)=(9,6,2)`. The same word gives the strong chain `9→31→97` and final graph

      (20,1,2,39,14,471,31,95,223),  new guard 8,
      final guard profile (97,53,9).

These statements specify literal induced extensions, not merely other tournaments with matching counts.

## 2. Protected guards on opposite sides

A *protected left guard* consists of q with a unique in-neighbor r, together with a vertex b that is a source of `X−q`. Dually, a *protected right guard* consists of s with a unique out-neighbor t, together with a vertex c that is a sink of `X−s`.

Assume X is strong and the four vertices q,b,s,c are distinct. Then the left and right defects are zero. Indeed b cannot equal r, because that would make b a source of X, contradicting strongness. Hence no Hamilton path of `X−q` starts at r. The right statement is dual.

Let `L` add a vertex beating `X−q` and losing to q, retaining q. Let `R` add a vertex losing to `X−s` and beating s, retaining s. Both preserve strongness and the two protected guards:

* after L, its new vertex is the new source of the deletion of q, while the old right witness c is still a sink after s is deleted;
* after R, its new vertex is the new sink of the deletion of s, while the old left witness remains a source after q is deleted.

The four protected vertices remain distinct, using the new witness on the side just extended. In particular both defects stay zero at every stage.

L and R commute up to swapping the two new labels. In either order the new left vertex beats the new right vertex, and all other new incidences agree. Let `X(k,l)=L^k R^l X`, with the original ports q and s retained.

### Bivariate count theorem

Set

    h=H(X),
    U=H(LX)−h,   V=H(RX)−h,
    W=H(LRX)−H(LX)−H(RX)+h.

For every `k,l≥0`,

    H(X(k,l)) = h + (2^k−1)U + (2^l−1)V
                   + (2^k−1)(2^l−1)W.                   (1)

**Proof.** At every point the left defect is zero, so D sends `(h,u,0)` to `(h+u,2u,0)`. For each fixed l, the count is affine in `2^k−1`. The dual statement makes it affine in `2^l−1` for each fixed k. Their commuting actions therefore give a bilinear polynomial in these two variables, determined by the four corners `(0,0),(1,0),(0,1),(1,1)`. Those corners are exactly the coefficients in (1). ∎

This theorem is an exact composition rule for two physical interfaces. It does not assume that arbitrary scalar equalities can be reused as source states.

## 3. An infinite family and three consecutive +2 steps

Take the S163 first-step output

    T=(8,1,19,6,43,15),  H(T)=17.

The protected left pair is `(q,b)=(4,5)` and the protected right pair is `(s,c)=(0,1)`. The unique in-neighbor of q is 2, and the unique out-neighbor of s is 3. The four corner counts are

    H(T)=17,   H(LT)=29,   H(RT)=29,   H(LRT)=49.

Thus `(h,U,V,W)=(17,12,12,8)`, and (1) simplifies to

    H(T(k,l)) = 8·2^(k+l)+4·2^k+4·2^l+1
              = 2(2^(k+1)+1)(2^(l+1)+1)−1.            (2)

All these tournaments are strong and have `6+k+l` vertices. Formula (2) is proved for all k,l by the interface theorem; the checker only supplies finite controls of its seed and construction.

Combining the initial S163 ear with two left additions and then two right additions gives

    5 → H(T(0,0))=17 → H(T(2,0))=53 → H(T(2,2))=161.

All three scalar steps are `h→3h+2`, and each adds two vertices. Taking the right additions before the left ones gives the same final tournament up to new labels. The second and third steps consume different protected sides; they are not repeated uses of the same moving-guard ear.

The analogous 31-seed

    T4=(20,1,2,39,14,87,31),
    left pair (5,6), right pair (1,2),

has corner data `(31,22,22,16)`, yielding the further exact strong family

    H(T4(k,l)) = 16·2^(k+l)+6·2^k+6·2^l+3.

It gives `31→97` by two additions on either side, but `T4(2,2)` has count 307, not 295. Thus the two-sided symmetry that gives a third +2 step does **not** silently give a third +4 step.

## 4. The observer boundary

For a separate old observer port p, define the familiar state by

    H(X^a_p)=h+2x,    H(X^aa_p)=h+2(x+y).

For the two named inputs, retaining p=3 and p=4 respectively, the final CDDC states are

    (53,27,36)=9(5,3,4)+(8,0,0),
    (97,90,99)=9(9,10,11)+(16,0,0).

They are not the marked target corrections `4j(2,−1,0)`. Thus repeated **scalar** offsets have been constructed, but neither a reusable all-branch marked bank nor arbitrarily many consecutive offsets follows. The two-parameter family is a source-level construction with two fixed physical ports; its existence is not an all-count theorem.

## Reproduction

Run `python3 offset_check.py` and `python3 -O offset_check.py`. Explicit failure checks remain active in both modes. Independent endpoint/subset Hamilton-path DP verifies all 54 unique-in-neighbor guard cases among labelled strong tournaments of orders 3 and 4, both C/D transfers, both four-vertex complementary profiles, the named offset chains and fixed observers, the protected pairs and commutation, and all 18 controls with `0≤k,l≤2` in the two families. The universal formulas follow from the proofs above, not from extrapolating the finite checks.
