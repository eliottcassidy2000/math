# Dense orthogonal closure outside dimensions2,4,8

Status: PROVED subject to independent audit; separately identified matrix-Lie-group strengthening of the elementary main theorem.

For even n>=6, n!=8, the topological closure of W_n on V={sum(x)=0} is O(V), so every nonzero orbit is dense on its sphere.

For adjacent positive supports the angle theta has 2cos(theta)=2−8/n nonintegral rational. The root-of-unity argument shows theta/pi is irrational. Thus the reflection product generates a dense subgroup of rotations in span(h,k), fixing its orthogonal complement. K=closure(W_n) contains every such full circle, so its Lie algebra contains h∧k for each edge of J(n,n/2). Here u∧v sends x to u<v,x>−v<u,x>.

We give the required Lie-generation argument explicitly. The Johnson graph is connected and its vectors span V. Start with an independent adjacent pair and U its two-dimensional span; the associated circle supplies so(U). Traverse a spanning tree, retaining in U the span of visited vectors. If the next vector v lies in U, no enlargement is needed. Otherwise its parent u is a nonzero previously visited vector in U. Write v=p+w with p in U and nonzero w perpendicular to U. Subtract u∧p, already in so(U), from the available u∧v to obtain u∧w. For every a in U perpendicular to u, commutation with a∧u gives, up to sign, ||u||²a∧w. Together these give all U∧w. These and so(U) span so(U+Rw). Induction yields so(V).

The closed subgroup theorem makes K a compact matrix Lie subgroup of O(V). Since its Lie algebra contains so(V), its identity component contains SO(V). K also contains a reflection and therefore equals O(V). Transitivity on each sphere yields the orbit conclusion. W_n itself is countable; equality is asserted only for its closure.
