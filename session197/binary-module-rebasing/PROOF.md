# Reversible, count-preserving rebasing between two binary modules

Session197. This is an all-order physical rebuilding operation, with a strict consequence for the next one-seam extension. It does not prove arbitrary-target coverage. The operation deletes one unexposed vertex and adds one vertex, so it is expressly outside the earlier induced birth-only grammar.

## 1. Actual graphs and retained vertices

For integers `t,j,m>=1`, let `T(t,j,m)` have blocks in the following labelled order:

    a, U=Tr_t, b, X=(x0->x1), Z=Tr_j, Y=(y0->y1), R=Tr_m.

The outside block `a->U->b` is transitive. The block `I=X->Z->Y` is transitive. Orient `I->R->(a,U,b)`. Between `I` and `(a,U,b)`, every edge points from `I` to the outside block except `a->X` and `b->Y`. These rules specify every edge. Mark the **last** vertex `r` of `R`.

The cyclic order

    a, X, Z, U, b, Y, R, a

is a directed Hamilton cycle, with each transitive block read in its own order. Hence every graph in the stated positive parameter range is strong.

The count formula is

    H(T(t,j,m)) = 30 * 2^(t+j) * 3^m + 1.              (1)

For fixed `s=t+j>=2` and fixed `m>=1`, there are `s-1` distinct rooted graphs in this family with the same order and Hamilton-path count. Indeed the marked last root has outdegree `t+2`, so different values of `t` cannot be rooted-isomorphic. We do not claim that every pair is nonisomorphic after forgetting the root.

The five-parameter positive clocks and the1122 clock families are predecessors of this construction. The change here is allowing the **middle outside vertex** to become a transitive module, then physically transferring vertices between that module and the inside zero gap. The new assertion is not merely a new scalar evaluation of an old one-parameter ray.

## 2. Why the count identity is an all-order theorem

We give the bounded-run argument explicitly. For a marked vertex `v`, let `O` and `I` be its out- and in-neighbors, and let the crossing graph consist of the directed edges `O->I`, considered bipartitely. In a Hamilton word after replacing `v` by `Tr_q`, deleting its nonempty clone runs leaves internal return fragments, each starting in `O` and ending in `I`. Each such fragment contains an `O->I` edge. The fragments are vertex-disjoint, so these edges form a matching. A crossing cover of size `d` therefore bounds the number of clone runs by `d+1`.

For a fixed outside skeleton with exactly `ell` clone runs, distributing the `q` distinct, internally ordered clones among those nonempty runs gives

    sum_{i=0}^ell (-1)^(ell-i) binom(ell,i) i^q.

Consequently, for `q>=1`, the path count lies in the span of `1,2^q,...,(d+1)^q`. The same upper bound applies when the Hamilton word is required to end at a fixed outside vertex. No value at an empty module is used here.

The following crossing covers remain valid for every choice of the other positive module sizes:

- A vertex of `U`: the singleton cover `{b}`.
- A vertex of `Z`: the singleton cover `{a}`.
- A vertex of `R`: the two-vertex cover `{a,b}`.

Other vertices in the same transitive module introduce no extra crossing edges. Thus the function `H(T(t,j,m))` belongs to the tensor span

    span(1,2^t) tensor span(1,2^j) tensor span(1,2^m,3^m).

This is a12-dimensional space. Invertible Vandermonde evaluation in one variable at a time shows that the complete grid `t,j in{1,2}`, `m in{1,2,3}` determines it uniquely. The exact graph counts on that grid are:

|t|j|H at m=1,2,3|
|--|--|--|
|1|1|361,1081,3241|
|1|2|721,2161,6481|
|2|1|721,2161,6481|
|2|2|1441,4321,12961|

The right side of(1) belongs to the same tensor span and has these values. This proves(1) at every positive parameter triple. The finite grid is a complete coefficient certificate after a uniform structural bound, not an extrapolation from an observed pattern.

## 3. The replacement operations are literal and reversible

Continuously expose the ordered four old ports

    [r,b,x0,x1].

A right birth `R(mask)` makes the newborn beat just the selected old ports and lose to all other old vertices. A left birth `L(mask)` makes the newborn lose just to the selected old ports and beat all other old vertices. After either birth, discard only the newborn's temporary exposure; the four old ports stay exposed. A bit is indexed by the displayed port list.

To move from `(t,j,m)` to `(t+1,j-1,m)`, require `j>=2` and perform:

1. Delete the first vertex of `Z`. It is not an exposed port.
2. Add a right-born vertex beating only `b`: mask2, discarding the newborn exposure.

The new vertex is the last vertex of the enlarged `U`: every old vertex except `b` beats it. Every edge among retained old vertices is unchanged. Reading the newborn after the old `U` and retaining the order of all other blocks gives exactly `T(t+1,j-1,m)`. In particular, the marked old root is preserved.

The inverse move, from `(t,j,m)` to `(t-1,j+1,m)`, requires `t>=2`:

1. Delete the last vertex of `U`, again not an exposed port.
2. Add a left-born vertex losing only to `x0,x1`: mask12, discarding the newborn exposure.

This is the first vertex of the enlarged `Z`. The two moves are inverse up to the indicated graph relabeling. They can be repeated in any legal order while both module sizes stay positive. They connect all `s-1` representations at fixed `s,m`.

By(1), each deletion step reduces the count to `(H+1)/2` and the following birth restores it by `2H_intermediate-1`. The intermediate graph is still in the positive family and is strong. Thus the construction is more explicit than a statement that two formulas happen to agree.

There is no old-vertex re-exposure. The only extension to the previous port language is the expressly permitted **deletion of a specified unexposed module vertex**. This is a new rebuilding language, not an induced extension and not a counterexample to induced-supergraph obstructions.

Formula(1) is valid for every `m`, but this packet does **not** claim that the four retained ports suffice for arbitrary subsequent root-cloning births: when `U` is large, such a right birth would need exceptional incidences to its vertices. The four-port closure claimed here is precisely for the two replacement operations above.

## 4. Exact marked-root responses on the m=1 slice

Write `(h,e,R_r)` for all Hamilton paths, paths ending at the marked root, and ordered two-path covers whose second nonempty path starts at that root, with the first path allowed to be empty. On the slice `m=1`,

    h = 90 * 2^(t+j) + 1,
    e = 30 * 2^j - 3,
    R_r = (90*3^t + 24*2^t - 27)*2^j + 6*3^t + 1.  (2)

Here is an all-order certificate for the response formulas. The ending count has the same uniform size-one crossing covers in `U` and `Z`, so it lies in `span(1,2^t) tensor span(1,2^j)`. Count its four grid values directly.

Next append a vertex `w` beating only the marked root, with all other old vertices beating `w`. The elementary retained-root one-seam identity gives

    H(T+w) = h-e+R_r.                                (3)

Equivalently, partition a Hamilton path at `w`: either it ends at `w`, or the portion after `w` begins at the root. The first case contributes `h-e`, because the preceding old path cannot end at the root. The latter is exactly the reserved two-cover count `R_r`: the first old path may be empty, and the second starts at the root. Every nonempty first path ends away from the root and therefore points to `w`.

In the augmented graph, a vertex of `U` has crossing cover `{b,w}`. The new crossing edges at `U` are covered by `w`, including `w->r`. A vertex of `Z` still has crossing cover `{a}`. Therefore its count belongs to `span(1,2^t,3^t) tensor span(1,2^j)`. The complete six-point grid is:

|t|j|(h,e,R_r)|H(T+w)|
|--|--|--|--|
|1|1|(361,57,601)|905|
|1|2|(721,117,1183)|1787|
|2|1|(721,57,1813)|2477|
|2|2|(1441,117,3571)|4895|
|3|1|(1441,57,5353)|6737|
|3|2|(2881,117,10543)|13307|

It uniquely gives

    H(T+w) = (90*3^t+114*2^t-57)*2^j + 6*3^t + 5.

Combining this with `h` and `e` proves(2), including every positive `t,j`.

## 5. Rebasing changes which next targets are obtainable

The inherited ordinary one-seam gadget of index `L>=1` sends a physical rooted profile to scalar

    F(L)=h g(L)+R_r L-e(L mod2),

where `g(1)=g(2)=1` and `1<=g(L)<=L`. This is the reviewed law in `session197/pair-entry-bridge/PROOF.md` and the low-cycle `ROOTED.md`; no physical inverse closure is assumed.

For every `j>=2`, compare

    A=T(1,j,1),       B=T(2,j-1,1).

The replacement in Section3 takes A to B, preserving

    h=180*2^j+1.

Their singleton-gadget children have counts

    N_A=441*2^j+23,
    N_B=1209*2^(j-1)+59.                              (4)

These are strong physical children. Adding a vertex beating only an old root to a nontrivial strong tournament preserves strongness.

**Neither next target is obtainable by any one ordinary final gadget from the other representation.** For B, `F_B(1)=N_B>N_A`; every `L>=2` has `F_B(L)>=L R_B>=2R_B>=N_B`. For A,

    F_A(1)=N_A<N_B<F_A(2)=762*2^j+39,
    N_B<3R_A=873*2^j+57.

Thus `L=1,2` fail and every `L>=3` is too large. The inequalities hold at `j=2` and follow directly from their positive displayed powers for all larger `j`.

This proves an infinite family of strict one-generation gains from count-preserving rebuilding, rather than just different formal boundary states. It does not exclude obtaining those targets through several gadgets or through another construction.

At `j=2`, the two actual H721 graphs are:

    A=[30,4,384,1014,998,967,903,771,515,7], root9;
    B=[62,12,8,384,1006,974,911,775,519,15], root9.

They have profiles `(721,117,1183)` and `(721,57,1813)`, giving children1787 and2477. Their sorted outdegree sequences differ, so even the underlying unrooted graphs are nonisomorphic. The verifier supplies the literal rows for all grid graphs and their children.

## 6. Verification and scope

`rebase.py` gives the full adjacency construction, exact delete/birth maps, retained-port checks, structural crossing covers and Hamilton cycles. `verify.py` checks the complete count and response tensor grids, finite controls of all graph maps and covers, and the one-gadget separation inequalities. The unchanged `core.py` provides direct endpoint subset inventories; normal and optimized runs use explicit exception checks.

The proof works only with positive modules. Empty-module boundaries are not inferred from the exponential formulas. The result supplies reversible physical rebasing within an explicit smooth-amplitude family. It neither covers arbitrary odd targets nor establishes a universal returning ternary digit interface. In particular, it does not bypass the global no-induced-H13-inside-H27 theorem: vertex deletion is essential to this new operation.
