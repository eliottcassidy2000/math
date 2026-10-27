# Exception debt is a paid native clock, not a bound on grammar size

Status: PROVED elementary potential bounds; FINITE-EXACT controls; no accessibility theorem.

Let f(T) be the minimum number of backward arcs of a finite tournament over all total vertex orders. For an order r write D_r(T) for those arcs, viewed as an undirected exception graph. Retain named vertices, active legal Hamilton paths, every shared port, source weights and the supplied integer clock separately.

## 1. The potential

For any directed Hamilton path P on n vertices, let tau_P reverse every arc not on P. For every rank r the exact frame transport is

D_reverse(r)(tau_P(T)) = D_r(T) symmetric-difference E(P).

Consequently f(tau_P(T)) <= f(T)+n-1. The same P remains legal after tau_P and tau_P is an involution, so

|f(tau_P(T))-f(T)| <= n-1.

The shared-port lane sharpens this tariff to h(n)=n-1-min(2,floor(n/2)). Indeed choose an optimal source rank with k=f(T) exceptions and s of those on P. If s>=1, the transported frame has at most k+n-3 exceptions. If s=0, P is increasing in the source rank and therefore is exactly that rank; in the reversed rank all n-1 path edges are backward. Swapping floor(n/2) disjoint consecutive rank pairs removes exactly that many backward edges without affecting other pairs. Thus f(tau_P T)<=k+h(n); apply the involution for the reverse inequality. This includes h(1)=h(2)=0 and h(3)=1; h(n)=n-3 for n>=4. The weaker n-1 tariff remains valid.

Full converse preserves f. A re-marking changes no arcs and preserves f, though it need not have zero physical cost in the user's chosen move language. For disjoint ordered join J(A,B), f(J)=f(A)+f(B): restriction proves the lower bound and block-ordered optimal ranks attain it. For the star join with backward cross arcs except the one forward hinge, ranks with B before A show

f(star(A,B)) <= f(A)+f(B)+1.

The exception graph's cycle rank beta(D)=|D|-n+c(D), where isolated vertices count as components, satisfies for n>=1

f(T)-n+1 <= min_r beta(D_r(T)) <= f(T).

Negative lower bounds are replaced by zero. Thus quadratic f forces quadratic minimum cycle debt, but no conclusion about all other representations follows.

## 2. A lower bound for every monotone compiler expression

Consider a finite expression built from n named singleton leaves, disjoint binary J/star joins, lawful Hamilton re-markings, full converses, and legal tau complements. There is no free duplication, induced-vertex deletion, hidden relabelling that merges ports, or extra primitive prime-substitution rule. Every internal complement i acts on its actual subobject with n_i vertices. Let s be the number of star nodes. Structural induction gives

f(T_out) <= s + sum_i h(n_i) <= s + sum_i(n_i-1).

This includes complements below later joins: the induction sums both children's bills. There are n-1 binary joins, hence s<=n-1, and each n_i<=n. The inequality is a necessary construction budget, not a sufficient reachability condition. Physical elementary relocation costs may only increase the actual clock; they cannot be subtracted from this bound. Source weights are carried unchanged because no step here merges source records.

For a fixed-size process beginning instead with an arbitrary forest frame on n vertices, the same bound holds with initial allowance n-1 and no joins. More generally use the actual initial f or certified exception count.

## 3. Ternary cyclic towers separate two meanings of economical

Define T_0 as a singleton and T_k=C3[T_(k-1),T_(k-1),T_(k-1)], n=3^k. The shared-port lane supplies the same construction independently; the following proof is self-contained.

Inductively decompose every internal block's arcs into directed triangles. For equal block size m, label blocks A_x,B_y,C_z and add triangles (A_x,B_y,C_(x+y mod m)) for all x,y in Z/m. Each cross-block arc occurs exactly once. All unordered pairs are thereby partitioned into n(n-1)/6 cyclic triangles. Every rank has at least one backward arc per triangle. Conversely block ranks A,B,C at every ternary node yield recurrence f<=m^2+3 f(T_(k-1))=n(n-1)/6. Therefore

f(T_k)=n(n-1)/6.

Every monotone J/star/tau expression that does realize this target must satisfy

sum_i h(n_i) >= (n-1)(n-6)/6,

and for n>=4 the number of complement operations must be at least max(0,ceil((n-1)(n-6)/(6(n-3)))). The weaker max(0,ceil((n-6)/6)) bound also follows for n>1; a negative right side is vacuous. The same lower bound applies when the starting n-vertex object has a forest frame. This does NOT prove that any such expression or fixed-size path exists. The grammar enlarged by an arbitrary C3 substitution primitive instead stores the tower in O(n) nodes (O(n log n) bits with arbitrary named labels), while its minimum rank-exception debt grows quadratically. Cheap description in the enlarged language is not a cheap native-complement clock in the smaller language.

Concrete sizes: n=9 has f=12 and cycle-debt lower bound4; n=27 has f=117 and lower bound91; n=81 has f=1080 and lower bound1000. The shared-port lane sharpens these particular lower bounds to exact minima for k>=2: split an optimal rank of block A into nonempty L,R across a backward edge, and order L,B,C,R. The cross exceptions are exactly L×C and R×B, totaling m². Each is connected and the chosen internal A exception joins them; all vertices therefore belong to one defect component. The feedback count is still optimal and beta=f−n+1. The checker here certifies the lower bounds and the FAS optimum, not this connected reframe; that explicit construction and its checks belong to the shared-port lane. Natural block rank need not minimize beta.

## 4. Verification and inherited scope

check_compiler_budget.py uses independent subset dynamic programming and permutation enumeration on all1099 labelled tournaments of orders1 through5; checks7887 legal path complements; verifies triangle packings and explicit rank upper certificates for k=0,...,4 (n<=81), and independently computes tower minimum f through n9. Normal and optimized outputs agree exactly. Explicit exceptions protect checks under -O.

The frame XOR rule and marked J/star operations are inherited from weighted_renewal_moving_modulus_20261026, T1/T2, pinned at92e7f12f3b95840cff92c5e6cf40b2ae8952cae4. No claim of a new general feedback-arc invariant or classical potential-method principle is made. The new application is the explicit paid clock lower bound and the representation/language distinction for these exact tournament operations.
