# Order-compatible Hamilton paths merge across an overlap

Status: PROVED candidate, pending independent audit. No eventual-coverage theorem is asserted.

Inheritance: Session197/lowcycle/root-amalgamation/ROOT_AMALGAMATION.md already gives the exact ordered shared-root law, and Session197/path-polynomials/ROOTED_ORDERED_AMALGAMATION.md gives its positive retained-state representation. Those formulas are not claimed as new here. The present argument is instead an order-preserving injection with arbitrary cross-edge orientations and arbitrary shared vertex set, plus the low forward-sector rigidity in Section5. No literature-priority claim is made.

Let T be a finite tournament. Let X,Y be vertex subsets with X union Y=V(T), and let S=X intersection Y. For each linear order pi of S, let h_X(pi) count Hamilton paths of T[X] whose restriction to S has order pi. The order pi need not itself be a directed path of T[S]. Define h_Y(pi) in the same way. The empty tournament has one empty path.

## 1. Overlap-compatible product injection

For arbitrary orientations between X minus S and Y minus S,

    H(T) >= sum_pi h_X(pi) h_Y(pi).                 (1)

Here is a direct injection proving the claim. First, two vertex-disjoint directed paths can be merged while preserving both internal vertex orders. Repeatedly compare their current first vertices and emit the one pointing to the other. If one path empties, emit the rest of the other. After a vertex is emitted, it points both to the next vertex on its own path, when present, and to the unchanged first vertex of the other path. Therefore each consecutive pair in the output is an arc. Both input paths are recovered by restriction to their vertex sets, so this is an injection from ordered input-path pairs to merged Hamilton paths.

Now fix a compatible pair P_X,P_Y inducing the same order pi=(s_1,...,s_k) on S. Split both paths at these shared vertices. In each of the k+1 open gaps, their two private subpaths are vertex-disjoint. Merge those subpaths by the preceding rule, and reinsert s_1,...,s_k between the resulting gaps. A nonempty gap begins with one of the two original first gap vertices, which is beaten by the preceding shared vertex when one exists. It ends with one of the two original last gap vertices, which beats the next shared vertex when one exists. If both input gaps are empty, the shared-vertex adjacency was already an edge of both input paths. Thus the full output is a Hamilton path in T.

Restriction to X and Y recovers the original pair. This also distinguishes inputs from different pi. The construction proves (1) without any assumption on cross edges. It includes S empty and X=Y as boundary cases.

## 2. One shared vertex is still multiplicative

If |S| is zero or one, there is only one order pi. Consequently

    H(T) >= H(T[X]) H(T[Y]).                       (2)

This holds for arbitrary cross-edge orientations, not merely an ordered join of the two private parts. Identifying one vertex therefore does not make a generic additive construction possible: a proposed output bounded by H(X)+H(Y)+C, for a fixed C and arbitrarily large independent input counts, contradicts (2).

The product can be sharp even for a strong output. Take vertices v,a,b,c,d with v->a->b->v and v->c->d->v, and orient {a,b}->{c,d}. Each three-vertex input has H=3, while their union has H=9. (Its rows in order v,a,b,c,d are [10,28,25,16,1].)

## 3. The first order disagreement is a real escape

For S={x,y}, (1) becomes

    H(T) >= h_X(xy)h_Y(xy)+h_X(yx)h_Y(yx).          (3)

The constituent total counts do not determine this expression. In particular the majority order of Hamilton paths need not agree with the common arc x->y. An exact six-vertex example has arc mask 91, with unordered pairs listed lexicographically and bit one meaning i->j. Its rows are [54,8,2,5,14,30]. For x=0,y=1, its sectors are h(01)=7 and h(10)=8.

Therefore a replacement of (3) by a claimed fixed-edge-majority product bound is not justified. A potentially additive construction must exploit at least two shared vertices and genuinely control these order sectors, together with the additional Hamilton paths not in the merging injection. Small overlap alone supplies no exact-count formula.

## 4. Scope and remaining problem

The map is an explicit path-merging injection. It preserves both input paths as subsequences and the order of every shared vertex. It generally misses other output Hamilton paths. Thus (1) is a lower bound, never an asserted exact graph compiler. It rules out unrestricted additive growth from zero- or one-vertex overlaps and exposes an order-profile sidecar for larger overlaps. It proves neither all large odd counts occur nor that a suitable two-overlap additive operation exists.

## 5. The first two forward sectors are rigid

Write f(T;u,v) for the number of Hamilton paths placing u before v when u->v is an arc. Then:

- f=1 if and only if T is transitive; in that case the reverse sector is zero.
- f=2 if and only if T consists of one directed triangle containing u,v, together with a transitive source prefix and transitive sink suffix. In that case H(T)=3 and the reverse sector is one.

To prove this, Hamilton paths of an induced subtournament extend injectively to the full tournament while preserving their vertex order. Insert each omitted vertex into the first allowable gap: put it first if it beats the first path vertex; otherwise choose the first later vertex it beats and insert immediately before that vertex, or append it if no such vertex exists. Restriction recovers the original path.

If a directed triangle C does not contain both u,v, merge its three Hamilton paths with the one path u,v using Section1. The overlap has at most one vertex, and the outputs all put u before v. Extend them to T. Thus f>=3. If a triangle does contain u,v, two of its three Hamilton paths place u before v, so f>=2. Every nontransitive tournament has a directed triangle. This proves the f=1 characterization.

For f=2, every triangle must contain both u,v. A nontrivial strongly connected component contains a triangle, so there is exactly one such component C and it contains u,v. Every vertex w of a nontrivial strong tournament lies in a directed triangle: if no edge goes from an out-neighbor of w to an in-neighbor of w, the out-neighbors cannot reach w. Therefore each w in C minus {u,v} lies in the triangle u->v->w->u. Put W=C minus {u,v}. We have v->W->u, and W is transitive because any triangle within W would miss u,v. All other strong components are singletons, hence neutral transitive padding.

If |W|=m>=1, exactly one Hamilton path places v before u, namely v,W,u in the transitive order of W. Every other Hamilton path has u immediately followed by v; each vertex of W can independently be assigned to its increasing block before u or after v. Thus the forward sector is 2^m and the reverse sector is one. The condition f=2 forces m=1. Conversely a triangle with the stated padding has sectors (2,1), proving the second characterization.

This prevents an arbitrary large reverse sector with forward sector one or two. It does not assert a bound for any larger fixed forward sector.
