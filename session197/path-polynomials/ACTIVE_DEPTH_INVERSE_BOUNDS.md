# Certified active-depth and internal-padding bounds for whole-module inverses

Session 197. **PROVED for arbitrary, including adaptively chosen, sequences of actual rooted contexts.** These are necessary conditions and exact pruning rules. They do not assert that every target is reachable or that arbitrary context templates have a finite size bound.

## 1. A weighted local inequality

For a rooted tournament (Q,r), put D=Q-r. Let its root triangles be the ordered pairs (o,i) with r->o->i->r, and define

    d_Q=H(D),
    w_Q=sum_(root triangles (o,i)) H(Q-{r,o,i}).

Empty tournaments have H=1. For every nonempty input T of order n, with h=H(T) and b=B_T(2), the physical whole-module substitution obeys

    H(Q[r:T]) >= d_Q h + w_Q b.                         (1)

A context is active exactly when w_Q>0. Every term in the defining sum is a positive integer, so w_Q is at least the number of root triangles. An active context consequently satisfies

    H(Q[r:T]) >= h+b.                                    (2)

### Exact operator justification

The reviewed universal skeleton formula is an identity on the entire polynomial space, with every skeleton contribution nonnegative in the extended central-factorial basis. Keep all skeletons with no crossing cycle: together they give multiplication by B_D(x). For each chosen root triangle retain the single crossing triangle and an arbitrary packing of exterior odd cycles disjoint from it. The triangle has one marked run and two exterior vertices, so its factor is x delta, with delta f(x)=f(x+1)-f(x-1). Its remaining exterior packings give B_(Q-{r,o,i})(x). Hence

    A_Q >= B_D(X) + X sum_(root triangles) B_(Q-{r,o,i})(X) delta. (3)

Here B_D(X) means multiplication by B_D(x), and the multiplication factors act after delta. The retained skeleton classes are disjoint, and the omitted contributions remain entrywise nonnegative in the central-factorial basis. This is not an ordinary-monomial inequality.

Physical nonempty input profiles have nonnegative central-factorial coordinates; evaluation at1 reads coordinate0 in the appropriate parity. In particular it is a positive functional on this cone. Since B_T(0)=0, evaluating (3) at1 yields exactly (1). This proves the inequality for arbitrary internal structure of T, not just transitive clones.

## 2. Two input-profile bounds

The two elementary bounds used below are

    B_T(2) >= (n+1)H(T),
    B_T(2) >= 2^n.                                      (4)

For the first, take each Hamilton word of T and each cut after0,1,...,n vertices. Put the prefix in labelled bin1 and the suffix in labelled bin2, allowing either to be empty. This is an injection into the two-bin configurations: concatenate the two labelled paths to recover the original word and cut. Thus there are at least (n+1)H(T) configurations. More generally the same fixed-label argument gives B_T(d)>=binom(n+d-1,d-1)H(T), but only d=2 is needed.

For the second, retain the empty odd-cycle collection in the bin identity. Its contribution at2 is 2^n and every other term is nonnegative. Combining (1)-(4), an active step has both bounds

    H(Q[r:T]) >= (n+2)H(T),
    H(Q[r:T]) >= H(T)+2^n.                               (5)

The stronger context-specific multiplier d_Q+w_Q(n+1) and the exact weighted budget in (1) may be retained in an implementation.

## 3. Adaptive active-depth cutoff

Start with an arbitrary nonempty tournament T_0 of order n_0 and count h_0. Apply any finite sequence of actual whole-module contexts; they may depend on all previous choices. Let s be the number of active steps, with any number of passive steps interleaved.

Every context is nondecreasing in Hamilton-path count: (1) gives H(Q[r:T])>=H(T), since d_Q>=1 and w_Q>=0. An active context adds at least two vertices. Before the j-th active step, the current order is therefore at least n_0+2(j-1). Passive steps may increase this order further and cannot reduce it. Applying (5) at the active steps proves

    H_final >= h_0 product_(j=1)^s (n_0+2j).              (6)

Thus for a prescribed positive target M the exact necessary active-depth cutoff is the largest integer s>=0 for which the product on the right is at most M. If M<h_0 there is no construction in this language. The cutoff is obtained by integer multiplication until the next factor would exceed M; no asymptotic estimate or numerical logarithm is required.

This applies to adaptive context choices. It neither removes nor transports internal passive padding. Such padding only strengthens the order bound at the next active step.

## 4. Internal-order and context-budget pruning

Suppose a partially built graph T is to feed a later active step in a construction whose final count is at most M. All intermediate counts are at most M by monotonicity. The second inequality in (5) gives

    2^|T| <= M-H(T) <= M-1.                              (7)

Therefore every actual input to an active node has order at most floor(log_2(M-1)). This floor is computed by the integer bit length, and the stronger bound with M-H(T) can be used when the current count is known. If M<=1, no active node is possible.

For a proposed context Q and a known input pair (h,b), equation(1) also gives the exact necessary weighted-triangle budget

    w_Q <= floor((M-d_Q h)/b),                            (8)

provided M-d_Q h>=0; otherwise the context is immediately excluded. Since w_Q is at least the number of root triangles, this also bounds the number of crossing edges O->I at the substituted root. It does not bound the number of exterior vertices with no such incidence.

The internal-order bound is not a bound on the final graph order. Arbitrarily many transitive vertices may be joined after the last active step without changing H. An unrestricted family of context templates is likewise not made finite merely by (7)-(8). A complete inverse algorithm must still specify and justify its template/terminal-shell enumeration.

## 5. Physical checks on padding and periodic scope

Transitive input Tr_n has H=1 and B(2)=2^n. The directed-triangle wrapper therefore has exactly

    H(D[Tr_n])=1+2^n.

For example Tr4 has count1 before wrapping and count17 afterward; deleting three transitive padding vertices beforehand leaves K1, whose wrapped count is3. Thus count-preserving internal padding cannot simply be stripped before a later active step. Inequality(7) deals with this correctly by retaining the actual order.

The first triangle iterations from K1 have counts1,3,15,123, while (6) gives successive lower bounds1,3,15,105. The bound is exact for the first two active steps and remains a proved lower bound thereafter.

If a fixed context period contains a>=1 active steps, after m periods formula(6) applies with s=am. It immediately gives superexponential growth. This strengthens the separate periodic-context dichotomy quantitatively and does not need a parity assumption or any Hankel-positivity assertion. A period with no active steps still has the exact ordered-join geometric formula from the passive decomposition.

## 6. Reproduction and dependencies

Run `python3 verify_active_depth_inverse_bounds.py` and its `-O` form. It reconstructs the exact weighted sub-operator for every rooted labelled context through order4, checks virtual columns and physical inequalities, independently enumerates literal fixed-label cuts for all labelled inputs through order5, checks integer target cutoffs and the internal-padding examples, and tests adaptive words of actual contexts.

Named frozen dependencies:

- `UNIVERSAL_CONTEXTS.md`, SHA256 `8615845e7bf06a6cdb1925fdbd8198a77b1a8864aa9c22704b0459678ff2ed51`: full-profile skeleton identity and positivity.
- `PROOF.md`, SHA256 `b88c86ac0afa64d3a1bf399c56fa7f8aad046fb93424c70a9ace8d164708183b`: physical bin profiles and central-factorial readout.
- `verify_extension.py`, SHA256 `3b537a9cbfea66aabe8b1d4cf59c84de45cc41d2ea6f80d053df8c083b1cad76`.
- `verify.py`, SHA256 `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`.

The cut injection is proved directly above. The separate periodic dichotomy is a consequence/companion, not a dependency of this proof or verifier.
