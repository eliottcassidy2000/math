# Partial-star growth: exposure flags, stabilizers, and symmetry phase changes

**PROVED:** elementary all-order statements below. **FINITE-EXACT:** exhaustive
checks on complete bases of orders 0–6, and unrooted partial-star groups on
bases of orders 0–4. No historical-priority or complexity-speedup claim.

## Inheritance and precise category

Session79 `ROOT_SWITCH.md` already proves that attachment masks modulo parent
automorphisms count **pointed**, not ordinary, child tournaments; its root-switch
operator Q has Q²=NQ and tr Q=N a_N. We retain that correction rather than
rediscover it. The October20 future-boundary packet already supplies an exact
finite-horizon Burnside response recurrence; nothing here replaces it with a
bounded-state scalar recurrence.

Fix a complete tournament T on n vertices with G=Aut(T). Adjoin a named root *.
A partial star is a word x in {?,−,+}^V: root-to-old arc unknown, old beats
root, or root beats old. Isomorphisms must preserve the root. Thus exact
partial classes are G-orbits, and H_x=G_x is the automorphism group of the
rooted partial object. Orientation reversal is NOT quotiented out.

## 1. Exposure-polynomial classification and generating function

Let c_j(g) count j-cycles in g. The exposure-layer polynomial is

    L_T(z) = (1/|G|) sum_g product_j (1+2 z^j)^{c_j(g)}.

Its coefficient z^r counts partial classes with exactly r exposed arcs.
Proof: a fixed coloring is constant on each cycle; that cycle is either all
unknown, or exposed with either of two orientations. This also gives the
in/out refinement product_j(1+u^j+v^j)^{c_j(g)}.

Write Z_Tourn for the ordinary tournament species cycle index, with empty
object included. The all-base-size formal generating function is

    P(t,z)=Z_Tourn(p_j=t^j(1+2z^j)).

Equivalently, using only odd partitions lambda,

    P(t,z)=sum_lambda 2^{e(lambda)}/z_lambda * t^{|lambda|}
                          * product_{j in lambda}(1+2z^j),
    e(lambda)=sum_j (j-1)/2 + sum_{unordered pairs a,b} gcd(a,b).

This is a formal series, not an analytic convergence assertion. At z=0 its
coefficient t^n is a_n; the diagonal coefficient t^n z^n is the pointed child
count p_(n+1), not (n+1)a_(n+1). The expression is a cycle-index specialization,
not a newly discovered Burnside formula.

## 2. Exposure metagraph incidence and exact mass

For orbit representatives x of rank r and y of rank r+1, let E_xy count the
literal ways to expose one unknown coordinate of x and land in orbit y.
Let D_yx count literal erasures from y landing in x. Let s_x=|Gx|.
Then

    sum_y E_xy=2(n-r),    sum_x D_yx=r+1,
    s_x E_xy=s_y D_yx.                                      (1)

Proof of (1): count adjacent labeled pairs in Gx × Gy in the two directions.
This is an exact finite-state recurrence once orbit representatives are
retained; rank and raw orbit counts alone do not determine its entries.

Let w_x be total weighted directed paths from the all-unknown state to x,
with each step weighted E. Then, at rank r,

    w_x=r! s_x.                                             (2)

Every labeled partial word in the orbit admits exactly r! exposure orders.
Thus traversing the quotient DAG without the inverse orbit-size correction
counts labeled exposure histories, not partial isomorphism classes. Dividing
only by r! still leaves orbit mass. The correct rank count is
sum_rank-r x w_x/(r!s_x). It is a check and a state-based recurrence, not a
shortcut that removes orbit classification.

## 3. The completion stabilizer is a FLAG, not a final quotient

For a partial x, let U_x be its unknown coordinates. The number of completions
modulo automorphisms preserving the partial exposure is

    F(x)=(1/|H_x|) sum_{h in H_x} 2^{c(h|U_x)}.              (3)

This counts pairs (full attachment, retained partial exposure pattern) up to
G. Forgetting the exposure pattern can merge distinct such completions.
The exact map is from H_x-orbits of extensions of x to G-orbits of full
words whose orbit contains an extension of x. It is surjective, not generally
injective. No stabilizer-only formula above performs this forgetting.

**Minimal hostile.** T=C3 with cyclic rows (2,4,1), x=(+,?,?). Here G=C3,
H_x=1. Four completions survive modulo H_x. But (++−) and (+−+) are in the
same G-orbit, leaving only three unflagged full-attachment classes. Bases of
orders 0,1,2 have trivial G, so no such collision occurs earlier.

A second quotient remains after this flag forgetting: full attachment classes
are rooted children. Session79's exact unpointing is

    a_(n+1)=1/(n+1) sum_[(T,S)] |Orb_Aut(T+*_S)(*)|.

It cannot be replaced by division of the number of full-attachment states by
n+1. The smallest root-weight hostile is already inherited at n=2: equal
attachment Hamming weight can yield transitive or cyclic triangles.

## 4. A sharp symmetry phase change for one-star defects

Now forget that the root was named, while keeping the actual partial directed
graph. Let d be its number of missing edges, all incident with the root.

**Theorem.**
* If d≥2, the unrooted partial automorphism group equals H_x and has odd order.
* If d=1, its order is |H_x| unless the two missing-edge endpoints have
  identical orientations to every other vertex. In that twin case its order
  is 2|H_x|, and its 2-adic valuation is exactly one.
* If d=0, it is the completed tournament group, again odd, but need not fix
  the old root.

Proof: for d≥2 the missing-edge graph is a star with unique degree-d center.
Every automorphism fixes this center and restricts injectively to G, whose
order is odd. A tournament has no involution: any transposed vertex pair
would have its arc reversed; Cauchy's theorem proves odd order.
For d=1, automorphisms preserve the missing pair, so the root stabilizer H_x
has index one or two. If index two occurs, the total group is even and has
an involution. Its action on the other vertices is an automorphism of their
complete induced tournament, so it fixes them all. It therefore swaps just
the missing endpoints, forcing them to be twins. Conversely twin endpoints
can be swapped. Since H_x is odd, the valuation is exactly one. The d=0
case is the same tournament oddness argument.

This refines the general defect-budget symmetry bound in the parent session:
a star with arbitrarily many missing edges never supports even symmetry
once d≥2. The only 2-primary phase is the final single missing edge between
twins; completing that edge destroys its transposition. This is about
**unrooted** partial graphs; root-preserving groups were odd all along.

## 5. Reproduction and exact finite evidence

Run:

    python3 check_partial_exposure.py
    python3 -O check_partial_exposure.py

All checks use explicit exceptions, not removable asserts. Both outputs
agree. Every base class through order6 is generated by extension followed
by complete degree-filtered canonicalization. Every ternary partial state,
every exposure/erasure incidence, and every completion is inspected. The
aggregate polynomial is also calculated independently by odd partitions,
not by the generated tournament bank. Unrooted partial groups are separately
enumerated through base order4. Universal results rest on the proofs above.

Aggregate exposure layers r=0,...,n:

| n | coefficients of L_n(z)=sum_[T] L_T(z) |
|---|---|
|0|1|
|1|1,2|
|2|1,4,4|
|3|2,8,16,12|
|4|4,24,64,88,48|
|5|12,96,352,688,704,296|
|6|56,592,2816,7360,11008,8912,3040|

Numbers of partial classes with completion-forgetting collisions are
0,0,0,2,12,90,792 for n=0,...,6. The exact receipt records the check count.
The verifier checks path identity (2) at full rank; its intermediate-rank
version is proved by the same labeled-history argument above.

## Outcome and boundary

The additional useful classifications are exposure rank plus exact ternary
orbit, stabilizer action on the unknown set, retained-versus-forgotten exposure
flag, and intrinsic-versus-marked center. These distinguish three operations
that otherwise look like one vertex-growth recurrence: expose edges, forget
the exposure flag, and forget the completed root. All three quotients are
needed if a structural generator is to produce A000568 rather than a mass or
pointed count. The closed all-order count remains the inherited Burnside
response recurrence; this packet supplies exact structural states and failure
controls for attempts to compress that recurrence geometrically.
