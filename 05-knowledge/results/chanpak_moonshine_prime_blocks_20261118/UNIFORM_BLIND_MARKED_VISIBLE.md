# Uniform contexts cannot see this switch; a marked contact programs its effect

**PROVED all-order context law and literal five-vertex control.** This combines the repository's slot-polynomial/three-probe machinery with the retained-root binary interface. The general inadequacy of unmarked counts and slot readings is inherited; the explicit low-order combination and its 2W response are recorded here without a minimality or historical novelty claim.

## Actual adjacency matrices

Rows are source vertices, columns destinations, labels0..4. Let

    L = 0 0 1 0 1       R = 0 0 0 0 1
        1 0 0 0 0           1 0 0 0 0
        0 1 0 0 0           1 1 0 0 0
        1 1 1 0 0           1 1 1 0 0
        0 1 1 1 0           0 1 1 1 0.

Their outgoing masks are(20,1,2,7,14) and(16,1,3,7,14). They differ only by reversing arc0→2. Exhaustive120-permutation comparison shows they are neither isomorphic nor converse-isomorphic.

Both have the same FULL optional path-slot polynomial

    f(q)=q^5+6q^3+2q.

For L the directed triangles are021,041,043; for R they are041,042,043. Both have exactly one directed5-cycle,04321. No two nontrivial odd cycles fit disjointly on five vertices. The inherited odd-cycle packing formula therefore gives the displayed polynomial: the empty packing contributesq^5, each triangle contributes2q^3, and the5-cycle contributes2q. Independent ordered path-cover convolution verifies its values separately.

## Exact all-order composition law for uniform contexts

Let C_k(T) count ordered covers by k NONEMPTY directed paths. Then

    f_T(q)=sum_k binom(q,k) C_k(T),
    C_k(T)=sum_(j=0)^k (-1)^(k-j) binom(k,j) f_T(j).

Now substitute T for one vertex of any fixed tournament host, giving every vertex of T identical contacts with each outside vertex. Compress each maximal T-run in a Hamiltonian path to one slot. For each outside skeleton with k slots, completions are exactly C_k(T). Thus

    H(host[T])=sum_k a_k C_k(T),  a_k nonnegative integers.

The same argument works for q ordered external path slots: order the internal runs first by external slot and then by position. Therefore the entire output slot polynomial depends only on f_T, not on the internal tournament representative.

Consequently L and R remain indistinguishable under ANY finite chain of such uniform one-module substitutions. In particular they have identical outputs under every P/Q/R probe word from `fifteen_probe_descent/README.md`, no matter how long. This is an exact contextual equality, not an asymptotic statement.

For reference the three inherited physical probe operators are

    P f(q)=q^3 f(q)+q^2[f(q+1)-f(q-1)],
    Q f(q)=q^3 f(q)+q[(q+1)f(q+1)+(1-q)f(q-1)],
    R f(q)=q^3 f(q)+2q^2[f(q+1)-f(q-1)]
              +q[f(q+1)+f(q-1)-2f(q)].

Their equality-preserving composition is immediate by induction. First-probe counts are93,177,243 for BOTH inputs. The skeleton argument explains why this is not peculiar to those three operators.

## One marked nonuniform contact exposes a programmable difference

Mark vertex3. Literal rooted profiles(h,e,t,d) are

    L:(9,1,23,5), R:(9,1,25,5).

Append a binary word using the retained-root operations:0 adds an ordinary sink;1 adds a vertex beating only the marked root. These contacts are deliberately NOT uniform across the substituted module. For the chronological bits read least significant first as W,

    H_L(W)=9+23W-5K(W)-(W mod2),
    H_R(W)=9+25W-5K(W)-(W mod2),
    H_R(W)-H_L(W)=2W.

The sign is R minus L. At W=1 the exact counts are31 and33; at W=2 they are55 and59; at W=3 they are67 and73. The original single arc0→2 is still the only arc differing between the two constructed graphs. Hence one arc switch initially invisible to EVERY unmarked uniform context can be made to contribute any prescribed nonnegative even difference2W by a marked context.

This does not select an arbitrary absolute center. Once W is chosen from the required gap, the lower count is still forced. It identifies a concrete resource for the selector: two physical sources with the same global polynomial can carry different marked response data. Uniform quotienting erases precisely that resource; adding more uniform probe words cannot recover it.

## Provenance and scope

The slot polynomial, odd-cycle formula, physical P/Q/R operators and their increasing read depth are inherited from `fifteen_probe_descent/README.md` and `third_cover_binary/README.md`. The binary transfer and K checksum are inherited from `ans_inverse/ANS_INVERSE.md`. The literal pair was found in this packet's small source search. The distinction between full slot equality and richer marked/deletion information is already a repository theme; no claim is made that this is its first occurrence. The present result isolates and independently checks this particular nonconverse five-vertex pair, states its uniform-context factorization, and combines it with an exact programmable2W marked response.

Replay: `python3 finite_edit_recontact/uniform_blind_replay.py`. This is standalone: it enumerates cycles, isomorphisms, rooted profiles and optional covers without importing earlier verification code; literally constructs the three uniform probes and marked-word outputs; and checks64 arbitrary three-outside-vertex uniform hosts.
