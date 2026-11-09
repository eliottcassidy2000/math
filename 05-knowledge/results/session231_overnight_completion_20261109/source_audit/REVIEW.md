# Independent audit of Session231 exact-bank source cylinders

Verdict: PASS for the arbitrary-integer inverse theorem and the supplied-bank capped executable. This is an exact membership algorithm, not a proof that the bank covers all odd targets. Physical realizability remains inherited from the rooted-bank theorem and audit. No finite residue calculation is promoted to integer coverage.

## Independent derivation

Use the forward bit state (q,h,p) initially (0,1,0). A bit 1 sends it to (2q+1,h,1); a bit 0 sends it to (2q,h+q+p,0). This also handles leading zeros. In a suffix of k bits, the zero at position j from the left contributes 2^(j-1) times the original u. Thus its total linear coefficient is C=2^k-1-rev_k(r). Only a zero in the first position sees the original parity, so its coefficient is E=[r<2^(k-1)]. Everything else is independent of u. Substituting u=0 gives constant g(r)-1. This proves the producer's identity, including r=0 and u=0.

For clarity, the elementary bound g(u)<=u+1 can be sharpened. Define f(0)=0, f(2s)=f(s), f(2s+1)=f(s)+s+(s mod 2). Induction gives g(u)+f(u)=u+1-(u mod 2), while f is nonnegative. Hence 1<=g(u)<=u+1. This avoids an overly naive even-step induction which would leave an unwanted parity term.

With a physical label (a,b,c), c>=a>=b>=0 and a positive, the constant Q=a(g(r)-1)-b e(r)+c(r+1) is nonnegative, and L=cM+aC>=cM>0. The exact equation becomes x=L u+a g(u)+aE(u mod 2). Its lower envelope is L u+a and its upper envelope is (L+a)u+2a. Therefore precisely the necessary interval printed by the producer follows. If x<a the interval has no solution. Otherwise its unclipped length is

 a(x+L-a)/(L(L+a)) <= aN/L^2+a/L <= 1/4+1/4.

The last inequalities follow from c^2 M^2>=4aN, cM>=4a and L>=cM. An interval of length at most one half has at most one integer. Clipping cannot increase that number. Every v>=0 has exactly one Euclidean decomposition v=Mu+r, so testing that candidate for every r is a complete membership decision. Equality, not congruence or nonempty interval alone, is required.

The least dyadic M satisfying both thresholds is O(1+a/c+sqrt(aN)/c). Direct bit evaluation costs O(log N) arithmetic steps per branch; reversal and other branch calculations cost at most that much. A fixed finite bank gives O_Y(sqrt N log N) arithmetic operations and O_Y(log N) working bits when streamed. This is numerical complexity with O_Y(log N)-bit operands, not polynomial bit complexity in the input length.

## Independent checker and implementation audit

`independent.py` evaluates g by the forward append state, independently of the producer's backward accumulator. It compares exact witness SETS for all 2,732 labels and every target 1..96 (262,272 comparisons), checks 4,000 randomly chosen wide suffix identities with source arguments up to about 88 bits, and compiles the original C++ with undefined-behavior sanitization. It verifies the executable output's cardinality, target ordering, bank membership and exact witness equality on every target 1..4096 plus two large easy-witness targets. The only missing odd targets in the checked small interval are 7 and 21. All printed CASES witnesses and both printed NONE cases were checked independently. The receipt records PASS.

At N=10^15 the checker bounds all 2,732 labels separately using M at that cap. Monotonicity bounds all smaller inputs. Conservative per-label envelopes are L<=(a+c)M-a, Q<=a(M-1)+cM, candidate v<=floor(N/c)+M, and evaluation terms <=(a+c)(v+1)+b. Maxima over labels appear in RECEIPT.json; even the evaluation envelope is below 2.01e15, far below signed 64-bit range. The producer's lo and hi calculations, ceiling-division addition, g accumulator and branch counter consequently fit as well. Threshold products fit signed 128 bits. Reversal increment is the standard reflected carry for sequential r and reaches zero after the last suffix. There is no unbounded source-length loop.

This safety conclusion requires the supplied validated bank and valid decimal input 1<=N<=10^15. The executable does not validate arbitrary external banks, missing arguments or malformed input, and must not be advertised as hardened against those inputs. None is part of its stated mathematical implementation contract.

The initial producer verifier omitted an explicit count/order check for subprocess output; the producer subsequently fixed it. The independent checker also checks count and order explicitly. CASES are witness evidence, not timings or interval coverage beyond 4096. Both original CHECK files agree and all producer SHA256SUMS entries match.

## Provenance and inherited scope

The bank is byte-identical to the Git blob `04-computation/experiments/rooted_triple_bank_order7_20261104.txt` at commit 9c48b9ac1e2a9c22b0814686c76e3358ad56051d. Its SHA-256 is 19f29d80b09f513d5ad36cc9537fb62a0037edd28bc304165955527dc8b5fd70. Cardinality, uniqueness and all physical cone inequalities were independently checked. This audit does not reconstruct rooted tournaments or independently count their Hamilton paths.

THM-4875(D) explicitly supplies the F2 carrier and e(0)=1. Its S83 addendum supplies the asserted 2^35 F2-only finite coverage but explicitly notes that cutoff was one run, not independently recomputed there. The present packet does not re-audit that cutoff.

Session217's `inverse/SUFFIX_INVERSE.md` already gives a square-root, magnitude-preserving inverse for the different fixed-core family, with a source-length loop. Its general strategy, limitations and numerical complexity are inherited. Session231's simpler closed suffix coefficient and elimination of that length loop are a specialization/adaptation, not a new global-coverage mechanism.

Session220's `selectors/PROOF.md` Sections 1–7 really concern this F2 family. The positive mixed-word drift rules out exact ordinary-binary rational witness relations, including asynchronous relations, unlabeled finite-bank choices and finite-state decodable recodings. It does not rule out output-support recognizers, arbitrary arithmetic registers, or non-rational encodings. The producer's discussion is correctly qualified by exact witness relations; read the sentence about residue automata in that context, not as a ban on every conceivable use of residues in a proof. The fixed-suffix dyadic-closure obstruction is also inherited there.

No all-target occupancy implication follows from these checks. The remaining statement is exactly that, for every admissible odd target, at least one physically valid label/cylinder candidate satisfies the integer equation. Neither necessary intervals nor modular hits discharge it.

## Portable replay

Run `python3 -I -B source_audit/independent.py` from any working directory. The checker resolves the producer packet as sibling `../source` after integration (with the original scratch layout as a fallback). Set MATH_REPO to the repository root if it cannot be discovered among the script’s parent directories. Provenance is frozen to commit 9c48b9ac1e2a9c22b0814686c76e3358ad56051d, not the moving HEAD. Normal and optimized receipts agree. The generated UBSan binary is a local replay artifact, not a publication input.
