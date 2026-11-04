# Primitive-order domains for odd perfect supports

Session 179. **Proved, using classical Bang–Zsigmondy.** This adds a sparse primitive-factor layer to Sessions 177–178, not a new nonexistence proof or a new prime-count record. Historical novelty is not established. The powerful input is classical; the contribution here is its exception-safe integration, exact algorithm, and test certificates.

## 1. Sparse-domain theorem

Let S be a set of s≥2 odd primes and q∈S, q≡1 (mod 4). A hypothetical odd perfect number with this support has exponents a_q≡1 (mod 4), a_p positive even otherwise. Set n_p=a_p+1 and

D_p = {ord_r(p): r∈S, r≠p}.

Define C_p to contain the n∈D_p with n>1, the appropriate parity (n≡2 mod4 at q; n odd otherwise), and

    every divisor d>1 of n belongs to D_p.                 (PC)

**Theorem.** Every hypothetical exponent length n_p belongs to C_p. In particular:

1. There are at most s−1 candidate lengths per source.
2. a_p≤max(S)−2.
3. τ(n_p)≤s, and 2^ω(n_p)≤s.
4. Every forced divisor L|n_p must satisfy (PC). More strongly, after any exponent or congruence restrictions, retain only candidates n∈C_p satisfying all those restrictions. An empty domain rejects the branch.

**Classical input.** For odd p and d>1, p^d−1 has a prime divisor absent from p^k−1 for all 1≤k<d, except possibly when d=2 and p+1 is a power of two. The other positive-base Bang–Zsigmondy exception is (p,d)=(2,6). We use the standard theorem, not an experimental generalization. See Győry and Smyth, *The divisibility of a^n−b^n by powers of n*, Integers 10 (2010), Theorem 13: https://math.colgate.edu/~integers/k27/k27.pdf .

**Proof.** Fix d>1 dividing n_p. If p≠q then n_p is odd, so d≠2. If p=q, q≡1 mod4 implies q+1≡2 mod4 and q+1≥6, hence q+1 is not a power of two. The base-two exception is unavailable. Thus a primitive prime r exists. It is odd because 2 already divides p−1. Also r≠p, and ord_r(p)=d. Since d|n_p, r divides p^{n_p}−1; since d>1, r does not divide p−1, so r divides σ(p^{a_p}). Perfectness implies r∈S. Hence d∈D_p. Taking d=n_p gives n_p∈D_p; taking all divisors gives (PC). The size bound follows from |D_p|≤s−1. Further n_p|r−1 for its top-level witness r, so n_p≤max(S)−1. Different d require different witness primes because a fixed r has only one order. Therefore τ(n_p)−1≤s−1. Finally τ(n)≥2^ω(n), and every divisor of a forced L also divides n_p. QED.

The bound is not uniform in S. Neither max(S) nor s is bounded globally here. The candidate-count bound is not a bit-complexity bound: constructing p^n and computing orders can be expensive.

## 2. What the earlier graph and LTE layers do not retain

Session 177 retains order parity and exact LCM implications among known target primes. Session 178 retains incoming valuation capacities and a finite numerical exponent interval. Neither alone explicitly demands a primitive target at **every cyclotomic divisor level**, including levels without any known in-support target.

The classical factorization

    σ(p^{n−1}) = product_{d|n, d>1} Φ_d(p)

explains the refinement. A primitive prime gives a witness labelled by its unique level d. Order-one denominator edges can exist, but cannot provide any of these witnesses. Thus an odd-order outgoing edge is weaker than a primitive outgoing edge.

For S={3,5,7,11,13}, the only possible Euler primes are 5 and 13. For either choice the S177 maximal graph is strongly connected. Yet the source p=7 has order labels

    ord_3(7)=1, ord_5(7)=4, ord_11(7)=10, ord_13(7)=12.

Its only allowed odd-order graph edge is 7→3, of order one. There is no odd order greater than one, so C_7 is empty. Both Euler branches are rejected immediately, without an LTE exponent box. S178 already rejects this five-prime support by valuation propagation; the new certificate identifies a simpler distinct missing primitive-level reason. This small support is already excluded by general known prime-count results, so it is only a mechanism comparison.

## 3. Two boundaries with exact arithmetic controls

**Top level alone is insufficient.** For source 5 and putative length n=9, prime 19 has ord_19(5)=9. But if S={5,13,19}, there is no target of order 3. Therefore n=9 is excluded by (PC), despite having a top-level witness. Directly,

    σ(5^8)=31·19·829,
    ord_31(5)=3, ord_19(5)=ord_829(5)=9.

The missing lower cyclotomic level is essential.

**All levels present is still insufficient.** For source 5, n=5 and S={3,5,7,11,13}, every nontrivial divisor level of n is represented by ord_11(5)=5. Nevertheless σ(5^4)=11·71 has an outside factor. Thus C_p is a sound domain, not a characterization of closed sigma rows.

**Role exceptions cannot be dropped.** σ(7)=8 supplies no odd primitive target of order 2; this is why the Euler-source residue is used. The familiar σ(2^5)=63=3²·7 has no prime of order 6 modulo which 2 is the base. Exact order modulo a composite (for instance ord_9(2)=6) does not replace a primitive prime witness.

## 4. Divisor ideals and a precise topological reformulation

Adjoin 1 and define I_p={m≥1: every d>1 dividing m lies in D_p}. This is a **finite divisor ideal**: m∈I_p and e|m imply e∈I_p. Its nontrivial elements are contained in D_p. Thus I_p is a finite multicomplex under prime-exponent coordinates. Applying the source role gives exactly C_p=I_p∩{role-compatible lengths}.

For squarefree lengths, define K_p on odd prime labels ℓ by saying a finite set T is a face when product(T)∈I_p. Divisor heredity makes K_p an abstract simplicial complex. Minimal nonfaces give irredundant conjunction obstructions: if every proper subproduct is permitted but the full product has a missing order level, that collection cannot all divide n_p. For Euler lengths, use the slice 2·product(T)∈I_q, provided 2∈I_q; this is again downward closed.

This is a combinatorial encoding, not a claimed homology or Euler-characteristic obstruction. It is useful for caching minimal forbidden conjunctions. The all-power multicomplex, rather than just its squarefree skeleton, is required to retain constraints such as the missing level 3 below n=9.

A tempting Hall matching construction gives **no extra strength** here: once d|n_p, every in-support target of order d is forced, not just a freely chosen witness. If all corresponding incoming lower loads meet capacities, one can select one witness from each nonempty level and automatically meet those capacities. The genuine extra content is nonemptiness of each required level, not matching terminology.

## 5. Complete exact decision without scanning an exponent interval

Compute C_p for every source. For each n∈C_p, construct σ(p^{n−1}) and repeatedly divide by 2 and each prime in S. Retain the row only if the remaining cofactor is 1 and the 2-adic valuation is 1 at q, 0 elsewhere. Store the remaining odd-prime valuations w_pr.

Select one retained row per source. Accept if and only if, for every r∈S,

    sum_{p∈S} w_pr = a_r.

There are at most (s−1)^s row combinations before further pruning. Necessity follows from the theorem and exact divisor-sum balance. Conversely the row closure and column balances give σ(N)=2N exactly, with no unidentified factors. This is a complete finite decision procedure for the prescribed support and Euler prime, just as S178's dense-box procedure was complete. The improvement is the provably sparse domain representation. The code deliberately strips known primes instead of requiring a general factorization routine.

Intersecting these domains with S178's LTE bounds and S177's branch requirements is sound and may be stronger than either representation separately. No runtime benchmark improvement is asserted beyond the explicit candidate-count bound.

## 6. Targeted repository recovery and provenance

The search deliberately avoided the prohibited mistakes ledger. Selected older threads, not bulk ingestion, supplied these tools and boundaries:

- **mac-mini-2026-10-14-wilson, THM-4640**, §5(a), `01-canon/theorems/THM-4640-the-repair-algebra-at-minus-one-collisions-are-increasing-two-three-representations-ladders-and-gersonides-carries.md`: uses Bang–Zsigmondy to rule out a purported new power-of-three cyclotomic quotient, retaining the exceptional base-two exponent six. This is the most direct recovered mathematical tool.
- **codex-planar-jc-dynamics-20260825, THM-4139/THM-4146**, especially THM-4146 §order-six finite-ring discussion, and `07-reflections/three-kinds-of-six-and-the-last-max-eight-wall-codex-20260825.md`: separates prime primitive divisors from lifted composite-modulus period. This supplies the hostile control ord_9(2)=6, with no primitive order-six prime.
- **mac-mini S75d (2026-06-28)**, `07-reflections/the-recursion-modes-are-cyclotomic-factors-depth-times-phi-d.md`: organizes x^n−1 by its divisor-indexed cyclotomic factors. Only that exact factorization is imported; its speculative LRC synthesis is not used as a theorem.
- **codex S645**, `07-reflections/perfect-aliquot-fixed-point-carrier-s645.md`: recommends retaining prime-power and residue data rather than inferring fixed-point rigidity from scalar near-perfectness. The row-domain theorem implements that direction.
- **codex S644**, `07-reflections/vieta-perfect-aliquot-carriers-s644.md`: distinguishes exact aliquot fixed-point data from triangular numerical coincidences. No triangular or tournament heuristic is used in the proof.
- **codex S540**, `07-reflections/lrc-support-flow-cut-flow-duality-s540.md`: separates existence of a flow from the fully assembled cancellation constraint. This motivated testing, and rejecting as redundant, the tempting primitive-packet Hall layer.
- **opus S704**, `07-reflections/the-cyclotomic-witness-tower-and-the-depth-wall-s704.md`: warns that per-configuration finiteness does not imply a uniform global theorem. Our support-specific sparse search has exactly that boundary.
- **S177/S178**: retain their Euler-aware graph, exact order-one denominator exception, LCM closure, LTE valuations, and exact row/column balance. The primitive existence theorem complements rather than replaces these.

These are conceptual provenance links, not evidence that unrelated conjectures transfer. Historical novelty of the OPN corollaries has not been surveyed exhaustively.

## 7. Reproduction

Run `python check_primitive_domains.py` and `python -O check_primitive_domains.py`. The standard-library test uses explicit exceptions rather than disableable assertions. It checks 189 support/Euler cases; independently preserves 124 directly verified closed rows; checks 129 primitive-factor levels by repeatedly stripping gcds with all earlier powers; confirms both strong-graph/empty-domain examples and both strict-boundary examples. `RESULTS.json` is deterministic. Finite tests corroborate the implementation; Bang–Zsigmondy and the proof establish the unbounded result.
