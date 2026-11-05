# Independent audit of arbitrary-gadget and complete binary-ready return barriers

Session 195, 2026-11-05. **PASS.** Both all-order proofs are mathematically valid, all replay stages completed, and no unresolved correction remains. These results exclude a marked return architecture, not the full scalar tournament spectrum.

## Frozen reviewed inputs

`SOURCE_SHA256SUMS` binds both source packets:

- One-gadget arbitrary-size barrier proof: `a4f2e2517d410a29a4b07dca70ab57bfe698a97e5927689ff8bcdb9380150389`.
- Complete all-word return proof: `540f4c7993103eae6d32149c0221d10d347647e4287b73dc563a7325d116e3e4`.

The full proof's stray equation reference was corrected to “by the K transition” before this snapshot. No mathematical claim or verifier changed in that correction.

## One-gadget theorem and all padding cases

The SCC factorization establishing the endpoint-complete dichotomy is valid: e=d puts the root in the last component; if its deletion there were strong, an out-neighbor could be a Hamilton endpoint and would prevent endpoint completeness. The S194 weak-deletion ratio therefore supplies h>=3d outside the universal-sink case.

For every physical gadget, the missing-first-bin identity

    B-L = sum_(b->u) L_u(G-b)

is an exact complement of the cover-deletion injection. It yields V>=d_G for a non-sink root, while the sink identities handle every remaining case, including the singleton. The stated K-offset and unified DK-offset decompositions are correct. Their nonnegative coefficients exclude either sign ±1 without any gadget-order or target cutoff.

The suffix reduction covers all finite K/D words, not just short suffixes. The first K follows D^m. Values m=0,1 use the explicit K and DK formulas; m>=2 gives offset at least the output deletion count. After that first K the root is endpoint-complete, and the padding lemma excludes all delayed repairs. Empty suffixes and words ending in D are treated separately. Universal sinks are handled by the divisible-coordinate argument, so that branch has not been lost.

The one-gadget sparsity corollary retains both necessary finite choices: a finite ready seed bank with h>=3 and a fixed finite gadget bank. It does not assert a polylogarithmic bound for arbitrary target-dependent gadgets in that first-reset setting.

## Complete all-word theorem

The additional endpoint-gap lemma is valid at every order: h-e is 0 or at least 3 when h>=9, by the terminal-SCC argument and positivity at every endpoint of a strong component.

The class with h>=9,d>=3 and each of e,c=R-h,f=R+e-h-d in `{0} union [3,infinity)` is invariant under K,D and every physical root switch. The zero-c D case correctly invokes endpoint completeness; the singleton root-switch case correctly invokes the endpoint-gap lemma. All other gadget cases use the missing-bin bound or the exact sink formulas.

Crucially, this gap class itself can contain ready states. The proof explicitly notes that fact and does not infer exclusion from invariance alone. Exclusion instead uses the **last** root switch:

- On gap inputs its GK and GDK offsets are nonnegative integer sums of gap coordinates, hence zero or at least three, except for the sink GK value at most -9.
- Longer initial D-blocks and all delayed padding repairs are excluded as in the one-gadget theorem.
- If the last switch is the first non-K operation, its input is still ready and the one-gadget theorem applies.
- If an earlier non-K operation occurred, the last switch starts in the invariant gap class.
- Words with no root switch but at least one D are covered by the padding lemma.

These cases exhaust every finite word, including adaptively chosen arbitrarily large gadgets. Therefore, from a physical binary-ready input with d>=3, the exact ready-to-ready language is K*, with unchanged sign and d and count `2^k(h+sigma)-sigma`.

The hypotheses and scope are essential. The physical C3 A,D,K reset at d=1 is retained as a positive control. The theorem does not cover arbitrary rerooting at an existing vertex, old-arc changes, reconstruction, or scalar targets reached at an unready root. The O(log X) ray consequence from finitely many mature ready seeds is valid even with arbitrary gadget choices, because the theorem has already excluded every non-K return.

## Completed executable review

Both frozen producer verifiers were rerun under normal and optimized Python. Each pair's stdout agrees byte for byte. The full verifier reproduces all 33,867 labelled tournaments through six, 202,013 rooted cases, 273 distinct control states, 165 gap inputs, and eight mature ready seeds. Its bounded marked-history saturation reproduces 16,519 retained histories and 1,139,811 checked transitions, including neutral transitions. This is a control for the proof, not a scalar coverage census.

`verify_independent.py` additionally uses the separately written combined-growth audit helper, not producer implementation code. It checks all 67 gadget profiles and the missing-bin identities, independently enumerates every rooted labelled tournament through five, and adds the actual H17 ready state. Its controls include:

- 1,599 endpoint-complete structural cases and 2,720 endpoint-gap cases;
- 414 gap invariance checks and 136 first-entry checks;
- both the genuine boundary state `(13,3,15,3)` with c=f=2 and `(17,3,21,3)`;
- 68,474 one-gadget/KD-padding controls through suffix length eight;
- 520,149 marked whole-word controls through length three, distinguishing pure K histories from histories containing a disruptive operation.

Independent normal/optimized receipts are identical. None of these bounded controls is used as the proof of an unbounded quantifier.

## Reproduction and release

Run `verify_independent.py --source /path/to/binary-return-barrier` under normal and optimized Python, comparing stdout. It imports the independent helper from `../combined-growth-audit/verify_independent.py` in the published layout. The producer replays use their own frozen standard-library scripts.

Publish `AUDIT.md`, `verify_independent.py`, `receipt.json`, `SOURCE_SHA256SUMS`, `replay-one-normal.json`, `replay-full-normal.json`, and `SHA256SUMS`. Frozen source copies, optimized duplicate receipts, and stdout support files need not be published.
