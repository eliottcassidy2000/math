# Gilbreath depth: record staircases, a fixed-gap obstruction, and corrected near misses

**Status:** PROVED deterministic statements and FINITE-EXACT prime-prefix computation. The actual prime-depth asymptotic remains OPEN. This note repairs specific logical and empirical claims in HYP-9336, `05-knowledge/hypotheses/HYP-9336-gilbreath-depth-tracks-the-maximal-prime-gap.md`.

## 1. The exact prefix staircase

Let P be a fixed prime list of length Nmax and a_r its absolute-difference row. After checking its leading entries are 1, let

    B_r = r+1+min{i≥1 : a_r(i)>2},

using zero-based columns and setting B_r=Nmax+1 if there is no defect. Then a prefix of N primes is unsafe at row r exactly when N≥B_r. Thus its sea-entry depth is the least r with N<B_r. The thresholds B_r are nondecreasing: a prefix already in the binary sea stays there. This turns one full-width row scan into an exact staircase for every prefix, not just selected powers of ten.

`prime_depth_records.py --limit 200000000` sieved all 11,078,937 primes below 200,000,000, evolved the integer rows until depth239, checked each computed leading value, and retained every threshold and its local source cone. The resulting files are `prime_depth_200000000.json` and `.out`. All arithmetic is integer; NumPy is used only for speed.

The familiar endpoint is reproduced: N=11,078,937, p_N=199,999,991, k=239, G=248. But the narrow ratio band inferred from six sample sizes is not a uniform description even inside that range:

- At N=31,569, p_N=370,609, k=69 and G=112, giving k/G=69/112≈0.616071.
- For N≥1,000 in the scanned range, the maximum ratio is35/34≈1.029412 and the minimum is69/112.
- For N≥1,000,000, the minimum is162/220≈0.736364. The depth162 plateau runs from N=1,094,426 to5,418,138 while the gap record grows from180 to220.
- The depth175 plateau runs from N=5,418,152 to10,655,462, while G grows from220 to234. At N=10,655,463, the gap248 triggers a depth jump to224, followed by later steps reaching239 at N=10,655,500.

The historical sampled band `[0.8,1.05]` is therefore sampling-sensitive and not a valid all-prefix bound. The broader proposed falsification band `[0.5,1.5]` survives for N≥1,000 in this computation; that is finite evidence, not a theorem. There is no claim that these statistics were previously unknown: Odlyzko's [1993 paper](https://www-users.cse.umn.edu/~odlyzko/doc/arch/gilbreath.conj.pdf) already investigates prime-difference depths.

A scalar standard-library implementation independently reproduces the key N=31,569 witness, as well as the N=1,000 and N=2,814 controls.

## 2. Maximum gap does not control general sea-entry depth, even with a perfectly safe edge

**Theorem.** For every N≥5, take the first N terms of

    2,3,5,9,15,17,19,21,23,... .

This is a strictly increasing integer sequence with first term2, every later term odd, maximal consecutive gap6, and left edge1 at every positive row. Nevertheless its sea-entry depth is exactly

    k(N)=N−3.

**Proof.** The gap row is

    1,2,4,6,2,2,2,...,

and its next row is

    1,2,2,4,0,0,0,... .

The two initial2s protect the leading1 by THM-4511, so the unit edge persists. The entry at zero-based column3 is4 and has an all-zero suffix to its right, so it remains4 for as long as that column exists. Hence no sea-entry occurs while the row length is at least4. Once only three columns remain, the protected columns0,1,2 contain1 followed by0/2 only. This is row N−3. The small N=5 case is the displayed length-three second row directly. ∎

This defeats every bound `k≤f(G)` based only on maximal gap G for the larger class of increasing prime-parity seeds, even if all unit-edge conclusions have already been verified. It does not refute a prime-specific relationship: the displayed seed contains composites, and prime arithmetic must supply any missing control.

**Connection contract.** Integrating a positive gap word preserves order, parity, maximal gap, and the entire difference triangle after its first row. It does not preserve primality or consecutiveness of primes. Those are precisely the missing sidecars. The equality family shows why one cannot omit them.

## 3. An isolated gap does not die globally in three rows

The inherited argument observed that a large gap g between neighboring gaps a,b produces neighboring differences g−a and g−b, whose central difference is |a−b|. That identifies one cancellation, not the lifetime of all descendants.

For a direct hostile, take an infinite gap row equal to2 except for one6 at column j≥1. Its first difference has4s at j−1 andj and zeros elsewhere. Thereafter it is the XOR automaton in units of4. In particular the entry at column j remains4 forever, because its right neighbor is always0. The central descendant may disappear while the outer descendants survive.

Thus the claim “an isolated large gap dies in three rows” and the inference “long persistence must come from clusters of unusual gaps” are not valid general CA statements. The actual prime environment might destroy isolated defects quickly, but that requires an environmental argument.

## 4. A later safe row does not certify the earlier edge

The assertion that `k(N)<N` alone proves every earlier leading entry is1 needs an additional prefix check. The nonconsecutive prime seed

    2,3,7,13

has left edge `(2,1,3,1)` and enters the safe singleton row at depth3<4. It has already failed at depth2. A correct finite certificate combines direct checking of all leading entries through the entry row with the absorbing-sea argument afterward. The prime computation in Section1 performs that earlier check.

No counterexample to the consecutive-prime conjecture is claimed. The example isolates a logical gap in the proposed certification implication.

## 5. Reproduction and next targets

Run from `04-computation/experiments/gilbreath_alphabet_three_20261019/`.

    python prime_depth_records.py --limit 200000000 > prime_depth_200000000.out
    python verify_depth_controls.py > depth_controls.json
    python -O verify_depth_controls.py > depth_controls_optimized.json
    cmp depth_controls.json depth_controls_optimized.json

The full-width scan uses exact integer arrays. Independent controls use a different standard-library sieve and scalar list differences. The fixed-gap equality family is checked at every length5 through1000, but its proof holds at all lengths.

The useful repaired question is prime-specific: which source arithmetic prevents the stationary XOR tails or other long-lived low-alphabet structures? The new alphabet-three shield supplies a finite certificate for bounded rows; it does not supply that prime arithmetic. No bound k≤G+O(1), asymptotic ratio, or all-N result is proved here.
