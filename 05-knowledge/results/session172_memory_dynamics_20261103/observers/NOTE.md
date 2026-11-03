# The future-observer metric: exact memory without arithmetic control

**Status: PROVED.** General contextual quotient and affine full-branch observer theorems, with finite-exact verification. These are elementary constructions; the Collatz parity conjugacy is classical, not a novelty claim. The proposed research connection is to treat response memory, precision consumption, and arithmetic height as separate resources. S171 supplies contextual response quotients; this packet makes their infinite dynamical analogue explicit. It does not prove Collatz, and makes no inference from 2-adic generic behavior to positive-integer behavior.

## 1. Nonlinear safe forgetting and the observer metric

Let X be a set, T:X→X a deterministic map, and o:X→A an observation into a finite alphabet. Define x~y when o(T^j x)=o(T^j y) for every j≥0. Define x~_k y when this equality holds for 0≤j<k, with ~_0 universal. Set R(x)=(o(T^j x))_{j≥0}.

**Theorem 1 (minimal autonomous response record).**

1. The relation ~ is the largest T-forward-invariant equivalence relation contained in equality of current observations. The quotient X/~ has a well-defined dynamics and observation.
2. If q:X→S is surjective and there exist B:S→S and d:S→A with qT=Bq and o=dq, then q(x)=q(y) implies x~y. Hence S surjects uniquely onto X/~, compatibly with the dynamics and observation. A finite exact autonomous observer exists iff X/~ is finite.
3. An arbitrary record preserving the next k outputs needs at least |R_k(X)| states, where R_k is the length-k prefix of R. Equality is attained by storing that prefix. This finite-horizon optimum need not admit a same-horizon autonomous update.
4. Let h(x,y) be the length of the common prefix of R(x),R(y), possibly infinity, and let d_R(x,y)=b^(-h(x,y)) for any fixed b>1, taking b^(-infinity)=0. Then d_R is a pseudoultrametric, a metric precisely when R is injective. If h(x,y)≥1 and is finite, d_R(Tx,Ty)=b d_R(x,y). Always d_R(Tx,Ty)≤b d_R(x,y).

**Proof.** Equality of infinite output strings persists after removing the first symbol. Conversely any forward-invariant relation within observation equality preserves every later observation by induction. This proves (1). An intertwining compression preserves every future output, giving (2). Distinct prefixes must have distinct records, proving (3). Two strings agreeing respectively with a third for h_1 and h_2 positions agree with each other for at least min(h_1,h_2) positions. This is the strong triangle inequality. Shifting a common prefix of positive finite length removes exactly one symbol; at prefix length zero the asserted bound is automatic. ∎

This is the deterministic, single-operation version of the S171 linear response-radical framework. In finite X, partition refinement eventually stabilizes. Infinite X admits no finite-dimensional stabilization argument. A finite output alphabet does not imply finite response memory.

## 2. A whole class of arithmetical systems has exactly the same observer geometry

Fix an integer b≥2. The ring Z_b is the inverse limit of Z/b^k Z, with metric d_b(x,y)=b^(-v_b(x-y)), where v_b(z) is the largest k with z∈b^k Z_b. This definition works also for composite b; no field assertion is needed.

For each residue a∈{0,...,b−1}, choose integers u_a,v_a satisfying gcd(u_a,b)=1 and u_a a+v_a≡0 mod b. Define

    T(x)=(u_a x+v_a)/b when x≡a mod b,
    o(x)=x mod b,
    g_a(y)=u_a^(-1)(b y−v_a).

Division by b on bZ_b is well-defined and unique. Each g_a maps Z_b bijectively onto the a-th residue cylinder and is an exact contraction:

    d_b(g_a(y),g_a(z))=b^(-1)d_b(y,z).

**Theorem 2 (full-branch observer conjugacy).** For every such family:

(a) Every word w=a_0...a_(k−1) is the exact length-k observation of one residue class modulo b^k. Its cylinder is

    g_(a_0)∘...∘g_(a_(k−1))(Z_b).

(b) Two points have the same first k outputs iff they are congruent modulo b^k. Thus d_R=d_b, with the metric base b, and there are exactly b^k length-k response classes.

(c) Every infinite word occurs for exactly one x∈Z_b. The response map R is an isometric bijection from Z_b to the one-sided b-symbol shift space with its prefix metric, and R T=shift R.

(d) No finite-state autonomous observer preserves all future residue outputs. Any exact length-k observer requires at least b^k states, or at least ceil(k log_2 b) fixed-length binary bits.

(e) The residue class modulo b^k determines T(x) modulo b^(k−1), for k≥1, but does not determine T(x) modulo b^k. Hence these optimal finite observations form precision-losing transitions, not same-level finite dynamical quotients.

**Proof.** The unit u_a has an inverse in Z_b. The congruence imposed on v_a gives g_a(y)≡a mod b; direct substitution shows Tg_a=id and g_aT=id on that branch. Multiplication by a b-adic unit preserves v_b, while multiplication by b raises it by one. A composition of k inverse branches therefore maps Z_b bijectively onto a single residue class modulo b^k. The b^k words have disjoint cylinders because their outputs differ, and together cover Z_b by taking the first k actual outputs. This proves (a) and (b).

The cylinders attached to successive prefixes of an infinite word are nested residue classes. Their intersection is a unique inverse-limit point. This proves surjectivity and injectivity in (c); the shift identity follows from the definition. Apply Theorem 1 to prove (d). For (e), if x≡y mod b^k then their first branches agree and T(x)−T(y)=u_a(x−y)/b, which is divisible by b^(k−1). But y=x+b^k is in the same input class and gives difference u_a b^(k−1), not divisible by b^k. ∎

One can express (c) by the numerical digit map Q(x)=sum_(j≥0)o(T^j x)b^j. Then Q is a b-adic isometry and Q T=S Q, where S deletes the least significant base-b digit. It is generally not a ring homomorphism or an order-preserving real map.

The inverse limit of finite response partitions is therefore a complete state record even though none of its positive-depth finite levels is an autonomous exact quotient. This distinction is the central memory lesson: progressively finer snapshots and a finite closed state machine are different promises.

## 3. Explicit finite-word arithmetic and a non-Collatz example

For a word w=a_0...a_(k−1), put U_0=1, C_0=0 and recursively

    U_(j+1)=u_(a_j) U_j,
    C_(j+1)=u_(a_j) C_j+b^j v_(a_j).

On that cylinder,

    T^k(x)=(U_k x+C_k)/b^k,
    x≡−C_k U_k^(-1) mod b^k.

The congruence is the word's unique cylinder residue: it follows alternatively by composing inverse branches, not by assuming that mere final divisibility independently verifies all intermediate branches. The inverse-branch proof establishes that equivalence.

For b=3 take (u_0,v_0)=(1,0), (u_1,v_1)=(2,1), (u_2,v_2)=(4,1). Then every ternary output word of length k occurs on exactly one class mod3^k, and no finite-state observer can reproduce all future residues. This is a literal ternary arithmetic dynamics, not a claim that all theorems are special to parity.

The coefficients affect archimedean growth and the arithmetic embedding of integers, but not the observer geometry in Theorem 2. More generally the inverse-branch proof applies whenever a compact complete ultrametric state space is partitioned into b branch images, and its inverse branches are exact contractions of ratio 1/b with disjoint first-level images at mutual distance 1. The affine hypotheses are a concrete, fully checkable realization.

## 4. Collatz, easy descent, and a positive counter-cycle in indistinguishable geometry

Consider three maps on Z_2, with even branch x/2 and odd branch (m x+1)/2, where m=1,3,5. All satisfy Theorem 2 and all are isometrically conjugate to the same binary shift by their own parity observers. Nevertheless their positive-integer arithmetic is different:

- For m=1, T(x)=ceil(x/2) for positive integers. Every x>1 strictly decreases, so every positive orbit reaches the fixed point 1.
- For m=3, this is the usual shortcut Collatz map, whose positive-integer convergence problem is not resolved here.
- For m=5, there is the exact positive cycle

      13 → 33 → 83 → 208 → 104 → 52 → 26 → 13.

  Thus convergence of every positive orbit to the cycle containing 1 is false for m=5.

The m=5 cycle has repeating parity word 1110000. In the m=3 system, the same periodic parity sequence belongs to the rational 2-adic point 19/101: its seven-step affine map is (27x+19)/128. This point is not a positive integer. In the m=5 system the seven-step map is (125x+39)/128, with fixed point 13.

These computations demonstrate, without relying on any unresolved assertion, that shift conjugacy does not preserve the distinguished positive-integer subset. The arithmetic content has moved into Q(Z_{>0}), which differs between maps. A universal full-shift theorem alone cannot establish the Collatz descent claim. Observer completeness is not a Lyapunov function.

## 5. Why finite parity cannot certify arbitrary tail properties on Z_2

Every finite binary prefix has extensions that are eventually alternating, eventually constant, and non-eventually-periodic. By Theorem 2 all occur as actual 2-adic points in its cylinder. Consequently neither eventual alternating behavior nor eventual periodicity can be decided from a fixed finite parity prefix uniformly over Z_2. The nonperiodic claim follows, for example, by appending bits with ones at factorial indices: their unbounded zero gaps and infinitely many ones preclude eventual periodicity.

For m=3 the alternating tails are exactly eventual entry to the 1↔2 cycle: Q(1) has bits 1010..., Q(2) has bits 0101..., and injectivity identifies every point whose shifted response equals one of them. Hence every 2-adic cylinder contains points eventually reaching that cycle and points not eventually reaching it.

**Crucial scope restriction.** This does not assert that every residue class contains both kinds of positive integer. The non-eventually-periodic or eventually alternating extension may invert to a nonintegral 2-adic value. A finite certificate consisting of an actually computed positive-integer orbit hitting 1 is entirely valid. The theorem rules out only a uniform fixed-depth decision of the tail property over the full 2-adic domain. It does not exclude variable-depth certificates, induction, or a new number-theoretic descent theorem.

Likewise the infinite-state obstruction concerns a finite autonomous record with exact parity readout at every future time. It does not prohibit an algorithm using unbounded integers, a finite program with unbounded working memory, an observer supplied with fresh external input, lossy or probabilistic prediction, or a finite-state model restricted to a finite invariant set.

## 6. Connection to tournament memory and a useful research split

S171's full-profile-blind storage and fixed contextual Hamilton read head can store and retrieve finite output words. The session’s Collatz construction develops that literal construction. The general theorem here supplies its missing contract: an exact k-symbol output tape stores one of b^k distinguishable finite responses; deleting its head is a correct update only at precision k−1. Keeping a full k-symbol forecast requires new information about the next symbol. A complete static tournament encoding is therefore not automatically a finite autonomous simulator of the original arithmetic.

Three resources should be measured separately in future constructions:

1. **Response memory:** number of distinguishable future output prefixes, exactly b^k here.
2. **Precision consumption:** one b-adic digit per dynamical step at finite resolution.
3. **Arithmetic control:** an independently proved bound or descent on positive-integer height.

The first two are completely solved for the affine full-branch family, including examples with incompatible positive-integer behavior. A proposed Collatz advance must contribute the third, or prove a new restriction on the arithmetic image Q(Z_{>0}); renaming parity memory or establishing symbolic completeness does not do so.

## 7. Sources and reproducibility

Sources: [S171 index](../../../../00-navigation/S171.md), [response-interface framework](../../session171_knots_information_20261103/framework/NOTE.md), [tournament information](../../session171_knots_information_20261103/information/NOTE.md), and [contextual recoverability](../../session171_knots_information_20261103/synthesis/NOTE.md). Classical precedent for the m=3 conjugacy: D. J. Bernstein and J. C. Lagarias, *The 3x+1 Conjugacy Map*, Canadian Journal of Mathematics 48 (1996), 1154–1169, [author-hosted original paper](https://websites.umich.edu/~lagarias/doc/bernstein.pdf). The present proofs are self-contained. No originality claim is made for the parity conjugacy or symbolic coding by contracting inverse branches.

Run `python3 check.py > observed.json`, then `cmp observed.json results.json`; also run with `python3 -O`. The checker independently verifies cylinder bijections, explicit residues, lower-precision updates, same-precision failure, prefix metrics on finite pairs, and the 5x+1 cycle. Finite tests are transcription controls for the all-order proofs, not evidence for Collatz convergence.
