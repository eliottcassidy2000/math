# Joint observation closure: why a linear response space is not a partition

**Status:** PROVED elementary finite-state formulation; classical partition-algebra and bisimulation mechanisms are credited. The explicit controls and applications are session contributions, not claims of historical priority.

## 1. Three different minimality questions

Let S be a finite state set and A_1,...,A_k rational weighted transition matrices. On an observable f:S→Q, use the pullback convention

    (A_i f)(s) = sum_t A_i(s,t) f(t).

Fix observed functions f_1,...,f_r, such as the level, final acceptance flag, or a named-vertex test. Include the constant function1.

A partition is **exact for this menu** when all observed functions are constant on its cells and, for every i and every target cell C, the weighted row sum sum_(t in C) A_i(s,t) is constant as s ranges over a source cell. This is weighted equitability. For nonnegative stochastic weights it is the usual strong-lumping/bisimulation requirement.

There are distinct objects:

- The smallest transition-invariant **linear** space containing the observations.
- The coarsest exact **partition**, whose whole space of cell-constant functions must be retained.
- A minimal representation for a specified scalar decision, which can be smaller still.

The first space need not be closed under pointwise multiplication and therefore need not be the second. Neither dimension is automatically an information-bit lower bound for arbitrary encodings.

## 2. Exact joint-closure theorem

Write f·g for pointwise multiplication. Let U be the smallest rational vector subspace of Q^S containing1 and all f_j, closed under every A_i, and closed under pointwise products.

**Theorem.** U is exactly the space of functions constant on the cells of the coarsest exact partition for the stated menu. In particular its dimension equals the minimum number of partition cells.

**Proof.** Any unital pointwise subalgebra V⊆Q^S defines an equivalence relation s~t when f(s)=f(t) for every f∈V. It is contained in the associated cell-constant space. Conversely choose a rational linear combination of a basis of V which has distinct values on all distinct cells: only finitely many proper hyperplanes of coefficient choices must be avoided. Lagrange interpolation in this one function produces every cell indicator, so V equals the full cell-constant space.

For such a partition algebra, A_i-invariance is equivalent to equitability, because the image of each cell indicator is precisely the weighted sum into that cell. Every exact partition algebra contains the seeds and is closed under all prescribed operations, so it contains U. Conversely U itself supplies an exact partition by the first argument. Thus its partition is the coarsest, and its dimension is the number of cells. QED.

This is a standard finite algebraic form of the partition/bisimulation correspondence. A primary published precedent for unital Hadamard algebras and partitions is [Characterizing imprimitive partition designs of binary Hamming graphs](https://doi.org/10.1016/j.ejc.2003.11.001). For broader weighted modal witnesses, see [Wißmann, Milius and Schröder, Quasilinear-time Computation of Generic Modal Witnesses for Behavioural Inequivalence](https://arxiv.org/abs/2203.11175). We use the explicit finite proof above, not an unexamined transfer of a more general theorem.

### Constructive closure

Start with the span of the seeds and1. Repeatedly add their transition images and pairwise pointwise products, then take rational span. Every strict step increases dimension, which is at most |S|. The procedure terminates and proves the theorem constructively. Ordinary partition refinement supplies a separate implementation: repeatedly split cells by all weighted row sums into existing cells.

This bound counts strict dimension increases, not bit complexity or total multiplication work. Rational coefficients can grow. If the field changes, equality and interpolation are interpreted in that field; a rational proof does not authorize ignoring modular cancellations of weights.

## 3. Two separately exact quotients need not have an exact common refinement

Use six states x,y,t00,t01,t10,t11. Targets are terminal. The only transitions are

    x -> t00 and t11, each weight1,
    y -> t01 and t10, each weight1.

Observe the source/target level and, separately, the first or second target bit. At x,y define both bit observations to be0.

The first-bit exact partition is

    {x,y}, {t00,t01}, {t10,t11};

the second-bit exact partition is

    {x,y}, {t00,t10}, {t01,t11}.

Both are equitable: x and y send weight1 to each target bit value. Their common refinement has five cells, keeping x,y together and separating all four targets. It is **not** equitable: weight into {t00} is1 from x and0 from y.

If a,b are the two bit observables, the joint query repairs this:

    A(a·b)(x)=1,        A(a·b)(y)=0.

The stable joint algebra has dimension6 and separates all states. The purely linear transition closure has dimension4. Thus two complete marginal weighted responses do not provide their joint response for free; one correlation distinction remains. Dividing the source transition weights by2 gives the same distinction for probabilities. This example concerns a shared branch choice, not independent sampling of the two coordinates.

It is also the exact reason a naive product of separately minimized carriers can fail: the compatible product cells must be **stabilized again** under the native operation. A CRT argument can justify joint reconstruction in a particular arithmetic system, but no such compatibility theorem follows from marginal minimality alone.

## 4. A variance query separates equal continuation counts

Use seven states x,y,p,q,r,t1,t0. Let t1 be accepting, t0 rejecting, and put

    p -> 2 t1,       q -> 2 t0,       r -> t1+t0,
    x -> p+q,       y -> 2r.

All arrows lower the level. Every ordinary linear continuation count from x and y agrees, including total paths and acceptance counts at each remaining length. Let f indicate t1 and g=Af. Then

    g(p)=2, g(q)=0, g(r)=1,
    A(g²)(x)=4,       A(g²)(y)=2.

The missing response is a second moment of a **shared intermediate state**. At source transition mass2 the means agree but the variances differ. A linear count can legitimately discard this correlation when no later operation tests it; a compositional carrier that permits the product query cannot.

This is a small explicit counterpart of the repository's run-order warning. It does not claim that the particular tournament obstruction is a probability model or that its missing coordinate is exactly this numerical variance. The map to that setting must preserve its actual seam and order operation, not just two unequal output counts.

## 5. How this applies to the folded-path and arithmetic results

The companion observer lane reconstructs the64-state folded recurrence. Its coarsest partition for level, final flag and weighted transitions has16 cells, while the corresponding transition-invariant linear observation space has rank12. Naming current vertex3 changes the allowed observer and yields35 cells. These are exact but different minima. A count-only equivalence need not be equitable.

The Session91 Zeckendorf theorem is a positive contrasting case: a proven invariant kernel and CRT compatibility of both registers close the joint operation. The Session91 hard-guard Gram is another contrast: its almost-full rank is a lower bound for a specified linear return interface, not for arbitrary scalar counts or finite-state partitions.

Thus the reusable composition test is explicit:

1. State the observation menu and weighted native transitions.
2. Compute transition closure if only linear responses are required.
3. Add pointwise products when joint observations of the same state are allowed.
4. Recompute closure when observers or operations are enlarged.
5. Retain source weights and actual alignment when composing the resulting objects.

This is an exact finite method, not an algorithm solving tournament isomorphism, LRC, or every open arithmetic problem.

## 6. Reproduction

`check_joint_closure.py` implements rational basis closure and weighted partition refinement independently. It compares them on every binary transition matrix through order3 with every Boolean observation, all acyclic binary matrices at order4 with every Boolean observation, and deterministic signed-weight larger controls. Both explicit joint/variance witnesses are checked. Normal and optimized Python runs must agree; checks remain active under `-O`.
