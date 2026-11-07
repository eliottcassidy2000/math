# A density-zero barrier for the universal binary cycle comb

**Status: PROVED all-order statements, using the inherited odd-cycle formula; FINITE-EXACT checks separately below. No literature-priority claim.** The Chan-Pak unrestricted Hamiltonian-path spectrum conjecture remains OPEN.

## 1. Inheritance and scope

The source family is Session212, `session212/cycle-bridge/PROOF.md`, sections 2--5, whose audited construction realizes every positive Hamiltonian-cycle count. We keep that definition and notation literally. The odd-cycle formula is inherited from `01-canon/theorems/THM-4876-deletion-and-flip-increments-of-the-hamiltonian-path-count-an-increment-of-two-forces-at-most-three-paths-the-increment-ladder-thm-013-for-all-orders-and-what-extension-spectra-forbid-for-chan-pak.md`; its Setting and item (4) state the formula and its older repository provenance. The source files and SHA-256 digests are pinned in `SOURCES.json`.

The closest proved mechanism is binary coding of a marked start count, followed by a sole-exit closure. The hostile is treating the total Hamiltonian-path count as that marked start count. The corrected near miss is the attempt to use Session212's complete cycle spectrum as a complete path spectrum. The least-used sidecar here is a pair of disjoint odd cycles, one at each apex.

The concept board is: binary code; total path count; odd-cycle collections; trailing zero bits; scalar growth; numerical density. Anchor: the all-odd path-spectrum question. Niche: quantify the gap between the two counting problems on the successful cycle family. Wildcard: a repeated total-path value whose binary-growth children separate.

**Concurrent independent proof.** Session213's context-transform packet independently obtains a stronger bound by two explicit classes of Hamiltonian paths: H >= (1+tau)2^m + 2 sum_i epsilon_i 3^i 2^(m-1-i). It implies the same density exponent with a smaller explicit constant. The density-zero conclusion is therefore a shared session result with two independent mechanisms. This packet retains the disjoint-odd-cycle certificate and the scalar-growth counterexample; it does not claim the stronger path-class bound as its own.

The result concerns only this one binary-comb family, including both values of its auxiliary bit and both orientations of the apex arc. It does not rule out all two-apex tournaments, nor other growth gadgets. Arbitrarily many additions of universal sources or sinks, or isomorphisms and tournament converse, do not enlarge this family's set of total path counts, so the same conclusion applies after those count-preserving operations. Arbitrary compositions or other attachments are outside scope.

## 2. Definitions

For a tournament T, let H(T) count directed Hamiltonian vertex words. Let HC(T) count directed Hamiltonian cycles modulo cyclic rotation. These are different observables.

For m >= 1 and q = sum_{i=0}^{m-1} epsilon_i 2^i, let B^s(m,q,tau) have a transitive chain x_0,...,x_{m+1} and two apexes a,b. Every earlier chain vertex beats every later one. The other arcs are:

- a -> x_1 and x_i -> a for i != 1;
- b -> x_0;
- x_1 -> b iff tau = 1;
- x_{i+2} -> b iff epsilon_i = 1;
- a -> b for s = +, and b -> a for s = -.

Here tau is 0 or 1 and 0 <= q < 2^m. The order is m+4.

Session212 proves

    HC(B^-) = q + tau,     HC(B^+) = q + tau + 1.

It also proves that inserting a new bit vertex immediately after x_1, with choice epsilon, sends (m,q,tau) to (m+1,2q+epsilon,tau) while retaining every old arc. In the reversed family with tau=0, these are the literal cycle maps C -> 2C and C -> 2C+1.

We use the known odd-cycle formula

    H(T) = sum_F 2^{|F|},

where F ranges over collections of pairwise vertex-disjoint directed odd cycles, each of length at least three, counted modulo cyclic rotation. The empty collection contributes 1. We do not re-prove or claim this inherited formula.

## 3. Two explicit lower bounds

**Lemma 1.** For every allowed m,q,tau,s,

    H(B^s(m,q,tau)) >= 1 + 2^m.                         (1)

**Proof.** For every odd-cardinality subset A of the m bit vertices x_2,...,x_{m+1}, there is the directed odd cycle

    a, x_1, A in increasing order, a.

It has length 2+|A| >= 3. There are 2^{m-1} such cycles. Their singleton collections and the empty collection contribute 2 * 2^{m-1} + 1 in the odd-cycle formula. All their arcs are present independently of q,tau,s. QED.

Suppose q > 0. Write j = floor(log_2 q), so epsilon_j = 1, and put t = m-j-1 >= 0. The terminal bit vertex is z = x_{j+2}. There are j bit vertices strictly before z and t bit vertices strictly after z.

**Lemma 2.** For every allowed q > 0,

    H(B^s(m,q,tau)) >= 1 + 2^m + W(j,t),               (2)

where

    W(j,t) = 2^t(3^j+1) - (1+(-1)^j) [t=0].          (3)

In particular,

    H(B^s(m,q,tau)) >= 2^t 3^j.                       (4)

**Proof.** In addition to the singleton collections from Lemma 1, take collections of two cycles

    (a, x_1, A in increasing order, a),
    (b, x_0, B in increasing order, z, b).

Here B is an even-cardinality subset of the j earlier bit vertices. A is an odd-cardinality subset of all bit vertices other than z and those in B. Thus A and B are disjoint. The first cycle is odd and has length at least three. The second cycle has length 3+|B|, also odd. Its closing arc z -> b is present because epsilon_j=1. None of the other bit orientations, tau, or the apex arc occurs in these cycles.

For each earlier bit vertex, choose A, B, or unused. For each later bit vertex, choose A or unused. Project onto A odd and B even. If 0^0 is read as 1 for this calculation, four parity evaluations give the exact number of choices

    P(j,t) = [3^j 2^t + 1^j 2^t - 1^j 0^t - (-1)^j 0^t]/4
           = W(j,t)/4.

Each choice determines a distinct unordered cycle collection, because its two cycles have different apexes and recover A,B. They are disjoint from the singleton and empty collections used before. Every two-cycle collection contributes 4, proving (2).

If t>=1, W(j,t)=2^t(3^j+1) >= 2^t3^j. If t=0, W(j,0)=3^j-(-1)^j >=3^j-1, and the baseline 1+2^m pays the possible deficit. Thus (4) holds also at j=t=0, where P=0. QED.

Both bounds are lower bounds, not formulas for the entire path count. Most odd-cycle collections have deliberately not been counted. The parity projection counts only the displayed certificates, not all cycle pairs.

## 4. Density-zero theorem

Let S_comb be the set of all integers H(B^s(m,q,tau)), over every parameter allowed above. Put alpha = log(2)/log(3), approximately 0.63093, and let

    N_comb(X) = #{h in S_comb : h <= X}.

**Theorem 3.** For every real X >= 3,

    N_comb(X) <= 4 floor(log_2(X-1))
                 + [8/(1-2^{-alpha})] X^alpha.        (5)

Consequently

    N_comb(X) = O(X^alpha),
    N_comb(X)/X -> 0,
    N_comb(X)/#{positive odd integers <= X} -> 0.      (6)

In particular the whole Session212 binary-comb family does not realize every sufficiently large odd integer as a total Hamiltonian-path count.

**Proof.** Count parameter presentations, which can only overcount distinct numerical values.

For q=0, Lemma 1 forces 1 <= m <= floor(log_2(X-1)). There are four choices of (tau,s), giving the first term of (5).

For q>0, partition presentations by the nonnegative t above the highest set bit. For highest bit j there are exactly 2^j binary codes q. By Lemma 2 any code whose total path count is <=X must satisfy 2^t3^j <=X. At fixed t, the number of possible q is therefore at most

    sum_{0 <= j <= floor(log_3(X/2^t))} 2^j
      <= 2 (X/2^t)^alpha.

An empty sum is zero. There are four choices of tau,s, and m is determined by j,t as m=j+t+1. Summing the upper bound over all t>=0 gives

    8 X^alpha sum_{t>=0}2^{-alpha t}
      = 8 X^alpha/(1-2^{-alpha}).

This proves (5). Since log X=O(X^alpha) and alpha<1, (6) follows. The same conclusion relative to odd integers holds because their counting function is asymptotic to X/2. QED.

The exponent is an upper bound, not an asserted optimal exponent or asymptotic formula. Distinct parameter choices can produce the same total path count.

## 5. A finite witness against scalar path-count growth

Even on this family, the binary growth operation cannot be encoded by a function of the total path count alone, retaining only m,tau,s and the new bit. In the reversed, tau=0 family,

    H(B^-(3,6,0)) = H(B^-(3,7,0)) = 81,

but inserting the zero bit gives

    H(B^-(4,12,0)) = 225,
    H(B^-(4,14,0)) = 245.

Thus no single-valued function f with H(new)=f(m,tau,s,epsilon,H(old)) can describe every such marked growth step. This is a finite-exact counterexample, not a claimed minimal collision or a no-go for richer states. Its four tournaments, adjacency rows, dynamic-program values and independent odd-cycle-formula values are in the receipt. The older two-apex transfer work already proves that richer state is generally needed; this witness specializes that distinction to the cycle-universal family itself.

## 6. Verification and consequences

Run from this directory:

    python3 verify.py > RESULT.replay.json
    cmp RESULT.json RESULT.replay.json
    python3 -O verify.py > RESULT.optimized.json
    cmp RESULT.json RESULT.optimized.json

The verifier uses explicit exceptions, so checks remain enabled under optimized Python. It checks every family member for m=1,...,8 by a subset dynamic program for total H and rooted HC; all 2,040 parameter combinations are included. Complete path permutations independently check all family members through order seven. A separate odd-cycle enumeration and collection recurrence independently checks all family members through order eight and the four growth-collision witnesses. It enumerates every proposed pair certificate for all j,t with j+t<=9 and validates its parity count and actual arcs. It also checks the two inequalities, the inherited HC formulas, the tournament encoding, old-arc preservation under binary insertion, and the collision.

The computation is finite implementation verification. The density-zero conclusion is proved by the explicit all-order cycle collections and geometric series above.

The result explains a concrete obstruction: the bits independently control the spanning-cycle count, but they also allow many smaller odd cycles to coexist. Their positive contributions force the unrestricted path count to grow on a ternary scale, too quickly for the binary code to cover a positive-density set. Any proof of path-spectrum completeness inspired by the cycle compiler must change the family or the gadget, rather than merely reinterpret its output count.
