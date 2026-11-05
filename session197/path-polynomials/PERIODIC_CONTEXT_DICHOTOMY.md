# A geometric-or-superexponential dichotomy for periodic whole-module contexts

Session 197. **PROVED for every fixed finite word of actual rooted tournament contexts.** This is a limitation on a physical construction language, not an obstruction to the tournament Hamilton-path spectrum.

## 1. Scope and statement

For a rooted tournament (Q,r), write C_Q[T]=Q[r:T]: replace r by the entire current nonempty tournament T, with every exterior incidence uniform on T. Fix a finite list of these contexts and repeat precisely that list, in the same order, m times. The initial tournament T_0 is arbitrary and nonempty. Put h=H(T_0) and h_m=H(T_m), so h_0=h.

Call a context active if its root belongs to a directed triangle, and passive otherwise. Then exactly one of the following descriptions applies:

1. Every context in the word is passive. There is a positive integer c such that h_m=h c^m for every m>=0.
2. Some context is active. For every m>=0,

       h_m >= h (1+floor(m/4))^m.                         (1)

   In particular h_m^(1/m) tends to infinity.

The empty word is passive, with c=1. Either parity of the total number of vertices added per period is allowed. The input is an intact module at each application; this theorem does not cover nonuniform births at exposed ports, changing the period, choosing a new context adaptively, or rebuilding a graph.

## 2. A positive triangular sub-operator in every active context

Use the full bin-profile and central-factorial conventions of `UNIVERSAL_CONTEXTS.md`. The nonempty odd and even basis vectors are

    q_k=x product_(j=1)^k(x^2-j^2),
    r_k=x^2 product_(j=1)^k(x^2-j^2).

Every physical nonempty B_T has nonnegative integer coordinates, and its coefficient at index0 is H(T). Evaluation at1 picks precisely that coefficient. The reviewed universal context theorem proves entrywise nonnegative integer matrices for every context, including on nonphysical basis columns.

Let N=|Q|-1. Its exact skeleton expansion always contains the packing in which all exterior vertices are unused. That contribution is X^N, where X multiplies by x.

If r,o,i,r is a directed triangle, retain additionally the skeleton consisting of that single crossing triangle and all other exterior vertices unused. It has N-2 unused vertices, one crossing cycle, one marked run and exterior length2. In the notation of the skeleton theorem its parity is epsilon=1, and L_(1,1)=delta, where delta f(x)=f(x+1)-f(x-1). Its contribution is exactly X^(N-1) delta, with coefficient1. All other skeleton contributions are nonnegative matrices in the extended central-factorial basis. Thus

    A_Q >= X^N+X^(N-1)delta = X^(N-2) D,                 (2)
    D=X^2+X delta.

The comparison is entrywise in the central-factorial basis, not coefficientwise in ordinary monomials. It is valid on every basis column. The X powers act after the shift; commuting them would be incorrect. An active context has N>=2, so the exponent N-2 is nonnegative.

For any context, active or passive, A_Q>=X^N by the same unused-exterior skeleton. Therefore a fixed context word containing an active context dominates a physical D/X word containing D: replace a selected active context by D followed by N-2 single-vertex paddings, and replace every other context by N paddings. Nonnegative matrix multiplication preserves this entrywise comparison, including across parity changes. Repeating the comparisons preserves it at every period.

## 3. A uniform path lower bound, including changing parity

Single-vertex padding has the exact cross-parity actions

    X q_k=r_k,
    X r_k=q_(k+1)+(k+1)^2 q_k.

The first matrix is the identity and the second is lower bidiagonal with diagonal (k+1)^2 and subdiagonal1. The triangle wrapper D has, in either parity, positive transitions from k to k+1, from k to k, and from k to k-1 when k>=1. Their weights are respectively

    1,
    k^2+6k+3 (odd) or k^2+6k+5 (even),
    2k^2(2k+1) (odd) or 2k(k+1)(2k+1) (even).

In particular its diagonal weight is at least (k+1)^2. For a period containing D, keep every operation at the current index except one chosen D. This proves that the complete one-period transfer, from either starting parity, has

- a k-to-k+1 entry at least1;
- a k-to-k-1 entry at least1 for k>=1;
- a k-to-k entry at least (k+1)^2.

The remaining diagonal factors are positive integers. These statements hold separately for both possible starting parities, so no parity-returning hypothesis is required. When the period changes parity, use the appropriate alternating matrices rather than a single fixed matrix.

The initial profile has at least h times its parity vacuum vector at index0. Fix 0<=k<=floor(m/2). In the m-period matrix product, retain the path which climbs one level per period for k periods, stays at level k for m-2k periods, then descends one level per period for k periods. All transition factors are at least1, and every stay has weight at least (k+1)^2. Consequently

    h_m >= h (k+1)^(2m-4k).                              (3)

Set k=floor(m/4). Since 2m-4k>=m, equation(1) follows. This is a positive matrix-path lower bound, not a claim that the intermediate basis states are physically realizable graphs.

## 4. Exact passive factorization

For a root r put O=N^+(r) and I=N^-(r). The root belongs to a directed triangle exactly when some cross edge is directed O->I. Therefore passivity is equivalent to every cross edge being I->O. For every substituted input T the graph is literally the ordered join

    Q[I] -> T -> Q[O].

A Hamilton path of an ordered join consists of a Hamilton path in each nonempty block in that order. Empty blocks contribute1. Hence

    H(C_Q[T])=H(Q[I]) H(Q[O]) H(T).

Multiplying these fixed positive integer factors over the period gives c and proves the passive alternative. This also shows why an outside cyclic component may produce a nontrivial geometric sequence: passivity does not mean the entire exterior is transitive.

## 5. Consequences for inverse construction

No fixed periodic whole-module word containing an active context can eventually have a count bounded by C b^m for fixed finite positive C,b. In particular it cannot eventually equal a finite sum of fixed-base exponentials, a polynomial times such exponentials, or a nonconstant affine clock

    h_(m+1)=a h_m+b

with fixed real a,b and finite initial value: every such recurrence has at most exponential growth (including the linear-growth case a=1).

A passive period may have an affine recurrence, but its exact formula is already h c^m. If this nonconstant sequence also satisfies h_(m+1)=a h_m+b eventually, comparing two consecutive distinct h_m forces a=c and b=0. Thus no nonconstant returning clock with a nonzero additive term is realized by any fixed whole-module period, active or passive.

This does not rule out a finite desired value or a graph assembled with a nonperiodic/adaptive word. It also does not rule out the known returning ternary port clocks: those operations do not retain the old whole graph as a uniformly substituted module. Distinguishing these languages is essential.

## 6. Reproduction and named dependencies

Run `python3 verify_periodic_context_dichotomy.py` and its `-O` form. The checker tests the exact sub-operator inequality on virtual odd/even columns of every rooted labelled context through order4, verifies passive full-profile factorizations, checks both parity versions of D/X period transitions, and directly builds and counts repeated physical context words, including a parity-switching active context and an exterior-cycle passive context. Finite controls supplement the all-order proof.

Frozen dependencies:

- `UNIVERSAL_CONTEXTS.md`, SHA256 `8615845e7bf06a6cdb1925fdbd8198a77b1a8864aa9c22704b0459678ff2ed51`.
- `PROOF.md`, SHA256 `b88c86ac0afa64d3a1bf399c56fa7f8aad046fb93424c70a9ace8d164708183b`.
- `verify_extension.py`, SHA256 `3b537a9cbfea66aabe8b1d4cf59c84de45cc41d2ea6f80d053df8c083b1cad76`.
- `verify.py`, SHA256 `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`.

This proof uses positivity and the literal triangular skeleton only. It does not rely on Hankel total nonnegativity, which may fail for an arbitrary initial profile.
