# Exact Collatz certificates beyond context-free slope languages

Status: elementary theorem and exact certificate construction; reproducible finite checks. Collatz remains open. No claim of literature priority. Sources were inspected in repository snapshot `1b91742ce55a202b73c80d56561ba56d00b503a6`. This is an algorithmic/resource lane, not an asserted major-problem proof.

## Useful outcome first

**Theorem.** Let D be the language of nonempty binary words w with `3^a < 2^L`, where a counts ones and L is word length. Every context-free sublanguage K contained in D has some constant `delta>0` such that

`3^a / 2^L <= exp(-delta L)` for every w in K.

Consequently D is not context-free, and its count set is not Presburger-definable. In particular, a finite family of context-free templates cannot represent the complete sharp slope-descent criterion. This is stronger than merely ruling out a fixed finite-state monitor, but much narrower than a lower bound for proving Collatz.

The constructive replacement is a compositional exact arithmetic certificate retaining `(A,B,C)` for the full branch map `(A n+C)/B`, its integer congruence, and the strict positive-height gate `(B-A)n>C`. All three coordinates are necessary to distinguish a descent from a positive cycle at threshold equality. A tiny certificate checker is included.

## 1. Proof of the grammar theorem

Write a Parikh vector as `v=(a,z)` with length `ell(v)=a+z`, and define

`h(v)=a log(3)-(a+z)log(2)`.

By Parikh's theorem, the count set of any context-free K is a finite union of linear sets `b + N p_1 + ... + N p_s` with nonnegative integer vectors. The cited theorem is used only for this semilinear conclusion.

Suppose every represented word belongs to D. Each base b belongs to that count set, so `h(b)<0` and b is nonzero. Each nonzero period p must have `h(p)<=0`, because `b+t p` belongs for every integer t>=0. In fact `h(p)<0`: equality would give `3^a=2^(a+z)`, impossible for nonzero nonnegative `(a,z)` by unique prime factorization.

There are finitely many bases and nonzero periods. Take delta to be the minimum of the positive numbers `-h(v)/ell(v)` over these vectors. Linearity gives `h(b+sum t_i p_i)<=-delta ell(b+sum t_i p_i)`. This proves the displayed uniform bound, even with prefactor 1. Zero periods can be discarded; the empty language is vacuous.

For each a>=1, take `L=bit_length(3^a)` and `w_a=1^a 0^(L-a)`. Unique factorization ensures

`3^a < 2^L < 2*3^a`, hence `-log(2)<h(w_a)<0`.

The lengths grow without bound, contradicting any uniform negative rate. Hence D is not context-free. The same argument directly rules out semilinearity of its count set, and thus Presburger definability.

**External source:** Rohit J. Parikh, *On Context-Free Languages*, Journal of the ACM 13(4), 570–581 (1966), https://doi.org/10.1145/321356.321364. The standard semilinear formulation was also verified against Esparza, Ganty, Kiefer, Luttenberger, *Parikh's Theorem: A simple and direct automaton construction*, https://arxiv.org/abs/1006.3825. The uniform-gap deduction here is proved above, not attributed as a quotation.

## 2. Exact sign-aware certificate calculus

Use shortcut Collatz `T(n)=n/2` for even n and `(3n+1)/2` for odd n. A chronological word w has triple

`A_w=3^(number of ones), B_w=2^(length), C_w>=0`.

The empty word has `(1,1,0)`, symbol 0 has `(1,2,0)`, and symbol 1 has `(3,2,1)`. For chronological concatenation uv,

`(A,B,C)_(uv)=(A_v A_u, B_v B_u, A_v C_u+B_u C_v)`.

Therefore a straight-line expression with concatenations and repeated subwords compiles using exact integer arithmetic; repetitions can use binary powering. This is a compression of branch composition, not a polynomial-bit-complexity claim. B alone has L+1 binary digits, so expanding a word whose length is exponential in its description may still require exponential space/time.

For a word of length L>=1, its unique compatible integer residue is

`r_w = -C_w * inverse(A_w) mod B_w`.

**Legality proof without simulation.** Final divisibility by B implies prefix divisibility: reduce the final numerator modulo the preceding power of two and use invertibility of the final branch multiplier 1 or 3. Induct backward. At the final step, integrality of y/2 forces y even, and integrality of (3y+1)/2 forces y odd. Induction therefore recovers exactly the claimed itinerary. Conversely a legal itinerary obviously makes every composition integral. Thus this one residue condition certifies every parity bit of the word.

For a positive integer n, the certificate is exactly

1. `n mod B_w = r_w`, and
2. `(B_w-A_w)n > C_w`.

For a triple authenticated by compiling the supplied word or straight-line expression, these conditions are equivalent to following w and ending at a positive integer strictly smaller than n. The low-level `certified(t,n)` routine assumes that compiler precondition; arbitrary externally supplied triples are not authenticated by its two arithmetic tests. An untrusted certificate must include its word/expression and have its triple recomputed. Positivity of intermediate states follows from the actual Collatz branches on positive integers. In particular, only A<B is possible. Conversely every word with A<B certifies all sufficiently large positive integers in its compatible residue class.

Let d=B-A>0. The smallest positive certified source is exactly

`n_min = r + B max(0, floor((C-d r)/(d B))+1)`.

The nonnegative threshold C/d automatically makes n_min positive, including r=0. Every larger source in the same progression is also certified; the preceding positive source is not. This formula makes the Archimedean side of a dyadic class explicit rather than silently dropping it.

**Hostiles.** For w=10, `(A,B,C)=(3,4,1)`, so n=1 gives equality and is correctly rejected; the first certified source is 5. Under the distinct `5n+1` shortcut, the actual 13-cycle has word 1110000 and triple `(125,128,39)`, again threshold equality at 13. Any test replacing strict descent by slope contraction alone wrongly certifies this nontrivial positive cycle. No assertion about the 5n+1 conjecture is required.

## 3. What the theorem does NOT prohibit

It does not rule out an intelligently chosen context-free language of sufficient certificates. For example,

`K={nonempty w : 2*(number of ones) <= 3*(number of zeros)}`

is context-free and has multiplier at most `(27/32)^(L/5)`. Indeed its odd density is at most 3/5, and 27<32. Its context-freeness follows from a pushdown counter for the signed weighted difference `2*ones-3*zeros`, accepting nonpositive difference. If every positive n>1 has a legal prefix in K satisfying the exact strict threshold, strong induction proves Collatz. Conversely Collatz implies that property: after reaching 1 or 2, extend along the 1↔2 cycle until `2*ones<=3*zeros`, and end at 1<n. Thus conservative grammar coverage can be equivalent to Collatz; the present theorem only excludes completeness for all sharp slope-contracting words.

It does not prove that long positive trajectories have high language complexity. It does not establish a new runtime lower bound for orbit evaluation. It does not identify a finitely checkable global coverage invariant. Arbitrarily long all-odd positive shadows remain a compulsory control.

## 4. Repository inspiration and disposition

All paths below are relative to the inspected math repository. These are precise source locations, not blanket endorsements of every assertion in a document.

1. `01-canon/theorems/THM-3480-rule30-staircase-transducer-entropy-and-nonrectangular-macroblock-compiler.md`, §§2–6: exact nonrectangular boundary state and charged preprocessing. Transfer: compile exact compositional summaries, account for bit growth. Do not transfer Rule-30 entropy to Collatz.
2. `01-canon/theorems/THM-3512-rule30-van-der-put-haar-cocycle-and-profinite-automaton-boundary.md`, status and §1: finite modular shadows need not give a finite exact synchronous machine. Transfer: distinguish precision-limited output from exact arithmetic state.
3. `01-canon/theorems/THM-2211-carry-regime-root-transducer-and-infinite-autonomous-index.md`, status and §1: finite controller with external carry stream is not an autonomous finite machine. Transfer: count unbounded integer sidecars as resources.
4. `05-knowledge/results/crossroads_20260926_automata.md`, §1: exact Collatz section `(k,u)` and native-prefix obstruction. Transfer: retain slope as well as carry. The pair `(2,2)/(1,2)` is a concrete carry-only collision.
5. `05-knowledge/results/nextforest_20260926_automata.md`, §§2–3: reachable residual states form ternary unit levels; a finite ternary fraction keeps denominator depth. Transfer: arithmetic compact coordinates are not finite-state quotients.
6. `05-knowledge/results/crossroads223_20260926_automata.md`, §§3–5: complete finite residue lifts preserve maximum cycle means, vector rotation sets, and exponential weighted path growth; word repetition preserves rational fixed-point quotient. Transfer: do not try to eliminate expansion by merely adding complete residue fibers.
7. `01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md`, mechanism and scope: finite polynomial-valuation corrections cannot repair a global logarithmic Lyapunov function. Transfer: selected returns and exact thresholds are outside the stated obstruction; arbitrary finite prime banks are not a new proof.
8. `01-canon/theorems/THM-4594-the-maximal-class-decided-collatz-sieve-and-the-sign-barrier.md`, §§4 and remaining conjecture: positivity and strict threshold matter; a positive cycle can sit exactly on a descent threshold. Transfer: include the strict Archimedean gate explicitly.
9. `01-canon/theorems/THM-4072-mahler-safe-terminal-fibre-product-and-finite-state-obstruction.md`, status and §1: safe symbolic language and eventually-zero native tape are distinct predicates. Transfer: grammar membership never substitutes for positive integer source legality.
10. `01-canon/theorems/THM-4654-safe-mahler-codes-contract-collatz-logarithmic-memory-and-complexity-approximation.md`, statements 4–7: real contraction supplies bounded entry, 2-adic separation supplies memory. Transfer: separate metrics and quantifiers; contraction alone is not integer coverage.
11. `01-canon/theorems/THM-3485-periodic-polynomial-fourier-jordan-recurrence-classification.md`, statement and §1: exact recurrence order needs all surviving Fourier/Jordan channels. Transfer: a fitted small recurrence is not a sufficient certificate; arithmetic composition identities are exact.
12. `01-canon/theorems/THM-1670-toral-recurrence-order-is-D-and-the-descent-is-order-2-only.md`, status and §(a): distinguish tested generic-order claims from uniform theorem wording. No universal minimal-order claim from that file is used here.
13. `01-canon/theorems/THM-4250-odd-core-suffix-lift-automaton-and-density-one-compiler.md`, status and §1: a sufficient finite carry compiler can be useful without being necessary or complete. Transfer: conservative grammar languages are legitimate even though D is not context-free.
14. `07-reflections/computational-irreducibility.md` and `07-reflections/computational-irreducibility-and-tournaments.md`: speculative analogy board only. Claims such as 'must simulate' or full isomorphism requiring all permutations are not accepted as complexity theorems. The useful prompt is to identify an exact compressible observable, which is `(A,B,C)` here.

## 5. Exact checks and open bottleneck

Run `python check_certificates.py`; repeat with `python -O check_certificates.py`. Both outputs byte-match the supplied JSON.

- All 16,382 binary words of lengths 1–13, including exact run-power compilation.
- 49,146 independently iterated source lifts, checking parity legality, affine endpoints and certificate equivalence.
- All 13,522 contracting words in that universe, with exact first certified positive source and minimality.
- First-descent certificates for n=2,...,10,000, including n=27 at shortcut depth 59 and n=703 at depth 81.
- Conservative `2*ones<=3*zeros` strict descent certificates for all n=2,...,10,000, with the exact fifth-power contraction envelope checked.
- 500 near-critical words generated without floating arithmetic, demonstrating the family used in the proof.
- 4,096 instances of one semilinear conservative family, checking a uniform exact envelope.
- Threshold-equality controls at the 3n+1 cycle 1 and 5n+1 cycle 13.

**Open bottleneck:** produce a finite, checkable *coverage argument* that every n>1 admits one of these strict descent certificates. Composition and verification are solved; global coverage is not. A proof-producing search over macro words, retaining their affine carries, may organize exploration, but enumerating more certificates or finding a compact grammar alone does not settle that bottleneck. The value here is a precise arithmetic interface and a proved boundary for attempts to make the fully sharp slope language context-free.
