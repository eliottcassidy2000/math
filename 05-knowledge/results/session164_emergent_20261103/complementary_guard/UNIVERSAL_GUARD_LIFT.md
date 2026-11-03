# Every strong source admits a finite guard lift

**Status.** A general first-ear identity and an exact bridge from arbitrary tournaments to the moving-guard algebra. The initial Hamilton count changes by an explicitly defined rooted path-cover slack. This does not provide a universal exact `3H` return, ternary offset closure, or an odd-spectrum classification. No novelty against the literature is asserted.

The attachment is S163's [moving-guard ear](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session163_ternary_triples_20261103/strong_returns/MOVING_GUARD.md). Its known unique-in-neighbor special case and defect are in [FIRST_STEP](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session163_ternary_triples_20261103/strong_returns/FIRST_STEP.md). The point here is to remove that hypothesis for the first attachment and identify its exact rooted-cover correction.

## 1. An unrestricted first-ear formula

Let X be any nonempty tournament and q any vertex. Let K(X,q) count ordered pairs `(P,Q)` of disjoint directed paths partitioning all vertices of X, where P is nonempty and ends at q, and Q is allowed to be empty. Put

    h=H(X),       J(X,q)=K(X,q)−h.

Every Hamilton path of X has a unique cut just after q, yielding such a pair. Hence J is nonnegative. More precisely J counts the pairs with Q nonempty for which the joining arc from q to the first vertex of Q points in the wrong direction.

Add a,b with

    q→a→b→q,      a,b→X−q,

and call the result F(X,q). Then, **without any in-degree assumption on q**,

    H(F(X,q)) = h+2K(X,q) = 3h+2J(X,q).                    (1)

**Proof.** Partition Hamilton paths by their first vertex. If a is first, b is forced second, and the remainder is any Hamilton path of X: h possibilities. If b is first, a must follow q, so the path is exactly `b,P,a,Q` for one pair counted by K. If neither new vertex is first, the unique possible new predecessors force `P,a,b,Q`, again with `(P,Q)` counted by K. All required joining arcs hold, including empty Q. These are disjoint and exhaustive bijections. ∎

The new a has unique in-neighbor q. Deleting a leaves b as a source, so the guard defect at a is zero. Paths starting at a number h. Therefore the new exact guard profile is

    (H, H−H_start(a), d_a) = (h+2K, 2K, 0).                (2)

If X is strong, so is F: the old strong component and the cycle `q→a→b→q` join all vertices into a strong component. This also holds for the one-vertex strong source.

When q has the unique in-neighbor r, a bad join must be `q,r`. Since r then occurs in Q, P cannot have a predecessor of q; hence P is the singleton q. Thus

    J(X,q)=H_start(X−q,r),

recovering precisely S163's defect. Formula (1) is a rooted two-path-cover version of that identity, rather than a new scalar rule inferred from examples.

## 2. A four-vertex canonical normalization

Apply the same ear again, using a as the current guard and retaining the first new vertex of this second ear as its next guard. Its input defect is zero, so the canonical final profile is

    (3m,2m,0),       m=h+2K=3h+2J.                        (3)

The construction adds four vertices to X and retains X as an induced subtournament. If X is strong, both stages are strong. Thus **every strong tournament has an explicit finite induced extension into the canonical moving-guard bank**, although its normalized count is (3), not necessarily `9h`.

The canonical returns in this session’s [guard algebra](../profile_algebras/GUARD_ALGEBRA.md) then give, for every `alpha,beta≥0`, strong induced extensions with

    H = 3^(alpha+1) 5^beta (3h+2J),
    order = |X|+4+2alpha+3beta.

This is useful source normalization with a known correction. It does not assert that the original h or observer state survives unchanged. Universal access to a guarded source also follows from the known triangle wrapper `G(X)=X→a→b→X`; the contribution here is the selectable rooted correction J and its relation to the full profile, rather than the first existence of some guarded extension.

## 3. The unpointed profile records total normalization slack

Let p_2(X) be the number of unordered covers of X by two nonempty directed paths. Summing K over all possible ports counts:

- h pairs with Q empty, one for each Hamilton path and its endpoint;
- 2p_2 pairs with both paths nonempty, one for each ordering of the two cover paths.

Consequently, for n=|X|,

    sum_q J(X,q) = 2p_2(X)−(n−1)h.                        (4)

In particular the right-hand side is nonnegative, and some port satisfies

    0 ≤ J(X,q) ≤ [2p_2(X)−(n−1)h]/n.

For comparison, the triangle wrapper has `H(G(X))=3h+2p_2(X)` and its b guard has zero defect. If n≥2, the average bound above is strictly less than p_2, so some rooted ear has `J<p_2` and a strictly smaller first-stage count at the same two-vertex cost. Both routes can then use one further ear to enter the canonical bank. This comparison concerns strong X so that the rooted route preserves strongness.

This connects the ordinary positive path-cover profile to the rooted guard lift: the second cover coefficient gives the total broken-join correction across all ports, while selecting a port determines its local share. A zero total implies every port has zero correction. None of these identities assumes arbitrary cover polynomials are tournament-realizable.

## Reproduction

Place `universal_lift_check.py` beside `offset_check.py`, then run it with ordinary Python and `python3 -O`. It independently computes K by summing endpoint/subset Hamilton-path counts over all subsets containing q. It checks (1), (2), the guard-specialized defect, the four-vertex canonical normalization, strongness when applicable, and the full-profile slack identity (4) on all 75 labelled tournaments of orders 1 through 4, totaling 285 arbitrary ports. Explicit failure tests remain active under optimization. The all-order results are the bijections above, not the finite census.
