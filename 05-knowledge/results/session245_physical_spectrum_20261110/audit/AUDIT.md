# Independent Session245 mathematical audit

**Verdict: PASS for the scoped proved statements. Chan–Pak remains OPEN.**
No substantive mathematical correction is required in the frozen supplied texts.
Audit covers interfaces/PROOF.md (including Section 8), forest/PROOF.md,
forest/RUN_BUDGET.md, descent/PROOF.md, and their proof-facing verifiers.
SOURCE_SHA256.txt records the final proof-facing sources using paths relative
to this packet. Publication normalized provenance wording and relative paths
without changing mathematical claims or executable algorithms.

## Interface mechanisms

The converse map exchanges A/G and S/E while preserving h,q,T; it retains
rather than complements the selected cut. The ordinal law and protected
closure use different operations. Appending a sink to the deletion graph
makes its endpoint cut coordinate zero and doubles A; only then adjoining
the selected-cut root gives (h+2A−T,h,h+2A,h). The empty separator is not an
identity for the endpoint formulas. The reversal theorem correctly requires
t>=1, and its cancellation equation and mixed E1 E2 term are correct.

Both order 19 reversal witnesses pass literal Hamilton-path and endpoint
subset DP in the supplied self-contained executable. The order 22 witness
is physically constructed, strongness checked, and composed from recounted
interfaces. Its full order 22 count has NOT been independently subset-recounted;
the proof and documentation correctly rely on composition there.

The independent.py audit uses permutation enumeration of paths on every
vertex subset, rather than the supplied endpoint subset recurrence. It
recounts every coordinate of both h13 and h37 core interfaces and then
recounts the closed graph's full count, root-ending count, root-start
ordered-two-path response, and deletion count. Results:

- h13 cut 3: (71,13,77,13), strong.
- h13 cut 6: (69,13,77,13), strong.
- h37 cut 17: (217,37,229,37), strong.
- h37 cut 33: (197,37,229,37), strong.

The h13 pair is not a pure full-interface T switch: G and S change.
Sink protection removes the dangerous E*S interaction. Constant A and h
then make all binary chain choices affine. The repeated h13 divisor law
and mixed h13/h37 gcd 2 proof are correct for all stated lengths. In the
mixed proof A2=47552≡7 mod 37 is the essential boundary calculation.

Section 8 is also valid: b_(j+1)=704 b_j+14492 h_j gives successive positive
coefficient ratio >704/37>19. The initial coefficients −74S,124S have
minimum subset-sum spacing 50S, S=37^(k−2); all subsequent differences are
larger by the proved domination bound. Exactly2^k distinct counts and
minimum spacing 50*37^(k−2) follow for every k>=2. Primitive gcd does not
imply binary subset-sum interval coverage.

## Fixed-apex letter spectrum

Each single-gap assignment acts on its own coordinate, with its eligibility
condition also depending only on that coordinate. It is an idempotent even
on unreachable states and with the zero-vector rejection convention.
Different gaps commute for the SAME letter. Simultaneous kernel/image
splitting proves the square-free annihilator with roots 0,...,k+1 over Q.
This is valid for arbitrary apex orientation and arbitrary neighborhood
letters, not merely the executable's star alphabet. Different letters need
not commute. The m=0 zero-eigenvalue correction is correctly retained.

Independent.py tests the gap-assignment bijection against literal path
permutations on all apex tournaments through k=3 and all length-two words
over the full 2^k-signature alphabet (549 cases, including k=0). It also tests
the full annihilator in nonempty prefix/suffix contexts for every letter
and every such apex tournament. Signature bits in this audit mean the leaf
beats the apex; Session221 uses the opposite bit convention. This is an
explicit encoding choice, not an identification of their raw bit words.

The top coefficient k! proof is sound: the product of gap projections fills
every gap with the bad nonfinal status, and each suffix vertex is then
forced to the unique gap preceding its matching root. This leaves exactly
one accepted assignment per root permutation. Finite formula checks do not
replace this all-order proof.

Direct source inspection confirms that the k=2 two-star formula is the
p=q=0 crossing specialization of THM-4748,
01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md.
Session221 05-knowledge/results/session221_computation_to_construction_20261107/automata/PROOF.md explicitly attributes the positive tape transfer
to Sessions215/216. The present texts accurately preserve that inheritance.
No exhaustive literature novelty search was performed in this audit.

## Run-budget deduction

Session221's finished-word SCC deletion preserves the apex core and H and
produces a subsequence of length<=L(N). A subsequence cannot increase run
count or introduce letters. Counting j-run words of length at most L gives
q(q−1)^(j−1) binom(L,j), so the stated distinct-count upper bound has the
correct inequality direction. For fixed q and r=o(log N), the logarithm
of this bound is o(log N). The empty representative contributes at most
one count, even when the fixed apex graph itself has count>N. A finite
number of fixed cores preserves the conclusion. This does not cover growing
alphabets, growing apex sets, unrestricted run budgets, or all tournaments.

## Synchronized-block descent

After O^k the physical state is (1,0,2*2^k,1), including k=0. The formal
coefficient extraction b+x*a is legitimate because the operations are
linear, without claiming either formal vector physical. The positive final
O block resets E and gives exactly the stated quadratic in x. For binary
middles, A=R is equivalent to palindrome by fixed-length binary uniqueness;
the pair-distance correction C includes the extra adjacent-pair term.
The empty middle and k=0 boundaries are handled correctly. The three
root-changing-middle examples have the displayed coefficient triples.
The inherited start-one invariant ensures that terminal J gives a single
Hamilton cycle and hence strongness. The executable independently confirms
one Hamilton cycle on all sixteen literal examples.

This is an exact family compiler, not a freely adjustable anchor or a
coverage theorem. The Mersenne family is correctly credited to THM-4748.

## Executable receipts

All supplied interfaces/verify.py, interfaces/verify_chain.py,
forest/verify.py, and descent/verify.py completed with PASS. Ordinary and
optimized receipts agree byte-for-byte. Both supplied --fail sentinels
(interfaces and descent) exited nonzero. One initial verify_chain comparison
used a pre-Section 8 source and therefore differed; it was rerun after the
source was frozen and matches. This was source drift, not an optimization
or mathematical failure.

Audit-owned independent.py and independent.json provide the additional
permutation and arbitrary-signature controls. Per-lane receipts and AGGREGATE.json retain the verified results; verify_all.py
repeats both modes. These tests have explicit finite universes and are not
claimed to establish the symbolic universal results.

## Final-source amendment: half-word norms and complete family inverse

A later descent update added Section 2's half-word formulas and target-first
inverse after the initial audit snapshot. The current source was reread,
the normal/optimized descent receipts replaced with fresh matching outputs,
and SOURCE_SHA256.txt regenerated. This amendment supersedes the earlier
snapshot for descent only; the other proofs were unchanged.

The added claims PASS. With a half-word of length t and chronological bits
v_i, the cross-half exponent is 2t−2−i−j, so its weighted sum is R²;
an inserted middle vertex doubles that cross-half term. The even central
adjacency contributes e. A central Z contributes 2R and two adjacency
corrections 2e. Internal-half contributions are 2C because reversal preserves
distances. These establish all three norm identities, including empty V.

The all-Z word maximizes C at fixed l. Its increments are 2^l, giving
C<=2^l−2 for l>=1. Combining l<=m−2 with the exact count formula gives
h−1<=2^m<=4(h−3)/3 for h>=9. For l=0 the exceptional m=2 output is 5;
remaining empty-middle outputs also obey the interval. No other output
below 9 exists. The interval's endpoint ratio is below 2, so it contains at
most one power of two. The code uses (h−2).bit_length() to find precisely
the least m with 2^m>=h−1. Its parity, length and defect tests exhaust all
palindromic middles capable of giving h. Summing their half-word counts is
geometric and O(2^(m/2))=O(sqrt(h)); per-candidate exact arithmetic is only
polynomial in log(h). This is a complete inverse for this synchronized
palindrome family, not the general O/Z/J middle family or all tournaments.

The updated verifier checks 511 mirror half-words and the inverse against
the forward image for all odd h<=5001. The forward controls' length<=10
and k<=5 exhaust this finite target range: inequality(7) forces m<=12,
then l<=10 and k<=5. All accepted recipes are reapplied physically.

The added forest/verify_run_budget.py was read and replayed in normal and
optimized modes with identical outputs: 24 exact word-count cases and 2730
literal fixed-apex deletion/run controls. Its signature convention is apex
beats leaf, and its SCC reachability deletion implements the stated kernel.
These are finite checks; the all-order bound remains the proved deduction.
