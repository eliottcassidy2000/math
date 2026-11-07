# Session223: two root changes, exact 3-adic factors and cutoff boundaries

**COMPLETE; independently audited results. Chan–Pak, the two-root-change cutoff, and the height-window multiplicity factor remain OPEN.** No new coverage census or literature-priority claim is made.

The session starts from main `59ff4cd82c66e89cc7abf47483956ae0a7e5fc51` and incorporates the recent root-change audit. Its main positive result is an exact all-precision local factor in a precisely stated ensemble, not a fit to the earlier height-window means.

## 1. Exact 3-adic multiplicity factor at fixed length

Sample a source word of length L uniformly among words containing exactly two root changes J. For each fixed k>=1 and nonzero residue r modulo q=3^k, put v=v_3(r). Then

    Pr(F=r mod q)=f(v)/q+O_k(1/L),
    f(0)=47/54,
    f(v)=157/108−3^(−v)(v/6+85/144), v>=1.

The zero-residue factor is

    Pr(q divides F)=C_k/q+O_k(1/L),
    C_k=157/108−3^(−k)(k/8+11/24).

The valuation factors strictly increase and approach 157/108. Their saturation ratio to the unit factor is **157/94**, about 1.6702. Modulo 3 the probabilities are **34/81, 47/162, 47/162**. The same leading law holds for the actual counted non-O terminal words of a fixed length with at most two Js, since the zero/one-J layers have only O(1/L) of their total mass.

**Mechanism:** classify every O/Z communicating class modulo 3^k by the retained d coordinate, g=gcd(d,U,3^k), and D modulo g. Root changes move mass between these classes. Two successive projections yield an explicit 3-adic probability density, whose geometric sums give the formulas. Periodic components are retained; the initial aperiodic block isolates the leading generating-function pole. The error is for fixed k, not uniform growing precision.

[All-k proof](threeadic/ALL_K.md) · [Short mod-3 derivation](threeadic/PROOF.md) · [Independent all-k derivation](audit/HIGHER_AUDIT.md).

Exact graph controls through modulus 81 agree, including **1,056,321 normalized states** at modulus 81. A second implementation independently verifies all residue probabilities modulo 3, 9 and 27. All normal and optimized replays agree.

## 2. Why this does not yet give the observed height-window factor

The earlier M_2(N) means were computed by numerical height, whereas the theorem samples by word length. The normalized real state determines which words enter a height window and is absent from a residue-only quotient.

A new exact hostile rules out a blanket independence shortcut. When the middle O/Z block B is nonempty, the real ratio just before the second J D/U is at most 1/2 precisely when the previous O/Z letter is O, and is greater than 1/2 when that letter is Z. Conditioning on these two cases gives different limiting mod-3 laws:

    O: (71/162, 91/324, 91/324),
    Z: (65/162, 97/324, 97/324).

Thus the final residue is correlated with real-state history. This does not decide independence from the final scalar height alone. A joint real/residue transfer theorem, integrating the actual height-window indicator, is still required. The fixed-length ratio 157/94 neither proves nor refutes the earlier finite-height ratio near three.

[Proof, independent calculation and exact scope](audit/THREEADIC_AUDIT.md).

## 3. The cutoff search: an exact equation and a real method boundary

Every exactly-two-J source has the form A J B J C, with A,B,C O/Z blocks. An exact middle-block formula gives a necessary-and-sufficient digital equation for membership, including all zero-length boundaries. It retains the middle binary integer, its reversal, a weighted digit-pair sum and an adjacency correction. These coordinates are coupled; treating their ranges as an independently fillable interval loses the target condition.

A separate all-order theorem shows that every mixed repeated O/Z block has a positive linear drift alongside its exponential growth, and that this drift survives every fixed physical suffix. Consequently any finite-state transducer that selects literal bounded-J witness words together with their exact binary output values selects only polynomially many in log X distinct values through X. It cannot certify positive-density coverage in this encoding. This extends the earlier fixed-bank obstruction to literal bounded-root-change words; asynchrony was already covered by Session220.

This does **not** refute a cutoff. Unbounded arithmetic carries, different witness encodings and other all-target arguments remain available. No constructive descent closing the exact digital equation, no new largest hole and no infinite family of holes was found.

[Digital equation, positive drift and finite-state proof](cutoff/REPORT.md) · [Independent review](audit/SCOPE.md#digital-cutoff-reformulation-and-selector-obstruction-pass).

The inherited independently reproduced census has **140 odd holes through 10^9**, largest **3,780,967**, with two source root changes. That remains finite evidence. A supplied source witness with at most two Js, closed by one more J, gives a strong one-cycle tournament of order at most floor(log2(N−2))+5. An effective tail theorem and small-target repairs are still needed for Chan–Pak.

## 4. Non-preserving histories, Fermat factors and square-sum cancellation

The full two-run rectangle satisfies

    F(Z^a O^b)=(2^a+1)(2^b+1), a,b>=1.

Swapping the runs preserves this scalar, but appending O separates the next values by **2^b−2^a**. The known chamber square is the equal-factor diagonal. Equal scalar values are not interchangeable rooted states.

The infinite square-sum graph on all factors **2^a+1** has exactly one edge, **3–33**. Its prime subfamily has none. Prime factors require power-of-two exponents, but the converse fails. This is a precise boundary to the Fermat/square-sum analogy, not a transfer of either prime conjecture.

The inherited full-state relation JZ=OJ does preserve every suffix. Its histories redistribute positive growth rather than subtracting births. In the square-sum graph, reflections x↦s²−x have an exact alternating-square cancellation law; closure still requires separate range, distinctness and connectivity checks to become a Hamiltonian construction.

[Proofs and controls](bridges/REPORT.md) · [Independent bridge review](audit/BRIDGE_AUDIT.md).

## 5. Reproduction, audits and continuity

Run `python3 verify_all.py`. It checks every published SHA256 and replays all nine verification programs in normal and optimized mode, including the modulus-81 graph. Every check remains active under optimization. The independent mathematical audits are recorded separately from replays of the authors' scripts.

The audit repaired an intermediate assumption of pointwise mixing: one mod-3 component has period two. The published resolvent proof handles it. It also separated square-sum switch legality from degree balance and credited Session220's prior asynchronous extension. These were draft corrections before promotion, not evidence against the final scoped theorems.

[Checkpoint 1: product rectangle and square-sum boundary](https://github.com/eliottcassidy2000/math/commit/8865bb666b3e6cebba43606d47b3c1b4557d522b). Final checkpoint adds the all-k formula, the real-state dependence theorem, the digital cutoff condition, independent audits and this synthesis. Publication uses non-forced expected-head updates and exact committed-file readback.

**Concurrent same-prompt work preserved:** commit `0cda2446321f7df69f6f2551580f44d496747cfe` added the [S15 preregistration](../root_change_local_factors_20261107/PREREGISTRATION.md), residue-automaton experiments, and planned height-window tests and a census to 2^32. It independently conjectured the mod-3 limits proved here. Its retained-letter budget n has the same dominant two-J ensemble at L=n+2; lower-layer weighting differs at finite length. The present all-k theorem supplies exact asymptotic constants, while its pending height tests remain separate. No unreported census outcome is assumed. The finite-automaton representation itself is inherited, not a new method claim.

Research posture: anchor = constructive cutoff; niche = actual modular multiplicity mechanism; wildcard = non-preserving histories and square-sum reflections. The live comparisons are full state versus scalar, words versus graph classes, length versus height, cancellation versus occupancy, and finite evidence versus all-order existence. The result is a proved local law and clearer constructive obligations, not a Chan–Pak proof.
