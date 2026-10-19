# Independent review: guard-refined first-crossing repair

Date: 2026-10-19. Verdict: **PASS**. Three minor wording corrections raised
during review have been applied and verified.

The main 25-episode / N+183 theorem and the separate 1000-episode / N+273
corollary are supported by the stated symbolic argument and independently
reproduced exact arithmetic. No gap in all-lift coverage or first-crossing
scope was found. This is a restricted bank-language result, not a general
Collatz or pointwise bank-entry theorem.

## Checks performed

- Both supplied normal and Python `-O` runs reproduce the stored receipt
  byte-for-byte.
- `independent_check.py` imports no supplied verifier. It derives affine
  constants by a closed sum over odd-bit positions, computes source residues
  by reversing the branches over the rationals, and replays every stored
  source guard as an integer parity word.
- All 1,066,572 source/offset pairs for the two base cylinders and offsets
  1 through 533,286 were recomputed. The 33,990 first-contracting parity
  classes, 15 exceptional offsets, complete offset digest, and exact threshold
  maxima agree with the published receipt.
- Main maximum repair lengths are 135 on A and 183 on B. The exact main
  threshold is 810692053755261919/1738366812781. The full base maximum is
  19411997879661734889626095/420491770248316829.
- All 838 stored forest nodes were independently reconstructed: 202 branching
  nodes, 275 offset-repair leaves, and 361 source-descent leaves. Every parent
  link, both roots per forest, all four children at each branch, and terminal
  status were checked. The count identity is 838 = 30 + 4*202.
- Every nonterminal and every recorded offset-repair leaf has multiplier at
  least one. Every source-descent leaf is an explicit F/H endpoint with
  contracting multiplier and strict all-lift descent. Thus offset guards in
  these actual forests occur strictly before the declared first crossing.
- Every offset leaf satisfies the authenticated guard, source-modulus
  divisibility, and strict original-source inequality. The largest repair
  length is 273 and largest source guard is 394 bits.
- Exact carry budgets, the expansion center of V, every nonempty K/V prefix
  multiplier, and finite controls for the negative rational K fixed point
  were recomputed.

The independent run performs 1,398,025 explicit checks. Its separate receipt
is `independent-receipt.json`.

## Symbolic coverage review

The partition is exhaustive. Starting in K, five K blocks or a V after one
through four K blocks supplies K^5. A reset after one through three K blocks,
or H after four, is one of the seven descending exits. The remaining K^4F
is precisely V and is expanding. Starting at this V, a V before three
additional K blocks or three additional K blocks supplies VK^3; the other
six possible early exits are VK^iF/H for i=0,1,2. All thirteen exits have
strict descent even at their least positive lawful residue. No lift or
short first-crossing word is omitted. The cylinders A and B differ in their
guard after the first 120 bits.

For each continuation of length j, source lifts differ by a multiple of
2^n with j<=n. Therefore their translated points x+d have the same j-bit
guard. The inequality (2^j-a)r > a*d+b at the least lawful residue then
extends monotonically to every positive lift. The numerical census is over
displacements, not an assumed exhaustive census of integer sources.

The finite-forest criterion is valid: a first-crossing word cannot finish at
a nonterminal node, because its slope is contracting, and cannot continue
past a source-descent leaf, because that would already be its first
crossing. Completeness and finiteness force a descent leaf or an offset
certificate within the word. In the general criterion an offset leaf may
also occur exactly at N; that remains valid. In all stored forests the
offset leaves are noncontracting and hence strictly before N.

The predecessor's carry argument remains applicable under the exact K^4F=V
rewrite. The rewrite preserves N, since all nonempty V prefixes expand and
the absorbed reset endpoint could not have been the first contraction.
The inherited zero-displacement minimum/least-period proof covers interior
steps as well as block boundaries and supports the distinct smaller
companion-cycle alternative.

## Rational K-infinity obstruction

The K fixed point is X=-4912087879/88519643. Canonical residues of X modulo
2^(30L) are nested, nondecreasing integers. Boundedness would make them
eventually constant, forcing X to be a nonnegative integer, a contradiction;
hence they are unbounded. Their translated residues agree with X+d to any
fixed parity depth for sufficiently large L.

A finite offset-repair certificate on the pure K path therefore gives a
contracting coefficient prefix of X+d. Conversely, any such contracting
prefix has a finite threshold, and sufficiently long K prefixes have both
enough guard bits and residues larger than that threshold. This proves the
claimed equivalence. Since 55<c<56, d>=56 gives the indicated positive
rational subproblem. The necessary-condition statement for universal finite
forests is correctly limited; it does not claim that this rational condition
is necessary for every possible proof about positive integer sources.

## Resolved wording corrections

1. Section 3 now correctly says the checker also runs the larger Section 5
   interval, rather than describing it as optional.
2. Section 5 now distinguishes the stored offset-leaf residue and threshold
   from the source-descent word and length, whose arithmetic is reconstructed
   by the verifier.
3. Section 4 now correctly allows the general criterion's offset guard to
   occur no later than N, including exactly at N. The actual stored forests
   have all offset guards strictly before N.

These prose corrections do not affect the certificates or theorem.

## Reproduction

After copying this directory into the result's `independent-review/`:

    python independent-review/independent_check.py
    python -O independent-review/independent_check.py

Or supply the result directory explicitly from any location:

    python independent-review/independent_check.py /path/to/collatz_guard_refined_repair_20261019

The reviewed supplied verifier SHA-256 is
`afde91bff09638a59fad89840783326c300aaa16f43487c9152c1313d8d0d9bb`.
The reviewed supplied receipt SHA-256 is
`b650a1cfa5803cb6235b2fba64aa458b0fd92400d25cf529d64c542fbe7f5998`.
The final reviewed proof SHA-256 is
`82552e278c5cd7cfac59942cd86beeb0f5d04695718e2e24076c2162b2464bfc`.
The complete base row digest is
`5708046e711978079dc77c04756291a4b460ee7b5c962a0a540755b79bfa37e0`.

No checkout files or Git state were modified by this review.
