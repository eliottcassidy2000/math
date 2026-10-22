# Session 26 source-bound archival repairs (proposal only)

No repository files or Git state were changed. `proposed.patch` contains four source-bound edits; `manifest.json` gives original SHA-256 and Git blob hashes and proposed SHA-256 hashes. Full original/proposed files accompany it. The publisher must verify hashes against its current target and reconcile any intervening changes. No mistake ID has been reserved; `mistake-entry-proposed.md` is one proposed entry for atomic publisher assignment.

## Sources and correction search

- Read all four source files in full from the local checkout (HEAD observed: 5887baff8ed); sampled files 19 and 23 independently identify the two reflections. Sampling manifest names commit a58ad91ac9050feb86c686a117759b90b1b37333.
- Searched `01-canon/MISTAKES.md` for reciprocal-odd, q_ct, hyperbola-triangle, q-gauge, Coxeter, and h+1/primality terminology. No targeted existing correction for A or B was found in that ledger. This is a targeted result, not a corpus-wide absence claim.
- Read `05-knowledge/results/hyperbola_triangle_q_gauge_s395.out` in full: the table labels the untruncated function as q_ct (its value at 2 is 1.5), and the global sinh wording repeats that conflation. Preserve this historical output, with explicit warning in THM-368. The cited `04-computation/hyperbola_triangle_q_gauge_s395.py` was absent from the local checkout; no regeneration or execution of it is claimed.
- Read MISTAKE-097, 098, 096, and neighboring corrections. 097 explicitly retracts both the deep-well census-invisibility claim and a subsequent compressed >=7-far uniform bounded-q claim. Only the former has a positive replacement at t=3/40. 098 rejects scale-based narrowing from weak chains. Do not repeat approximate q~log M empirical scalings as a universal proved upper bound.
- C was also fetched read-only through GitHub default-branch connector: https://github.com/eliottcassidy2000/math/blob/main/07-reflections/the-endgame-is-lrc13-plus-the-compressed-case.md, blob cbc4d06cb312f87d763e1f8d36a9774330aed279. It matches the local blob and has NO correction banner. The root handoff initially supplied this path for the runner's different banner sample; that attribution is now resolved. C inserts a separate warning, while D repairs the actual banner in the 13-comb reflection.

## Independent checks and preserved scope

Run `python check_math.py`; saved output `checks.json`. Exact rational checks verify q_ct's counterexample and untruncated parity, and the deep well's 3/40 witness. Cutoff identity checked numerically in addition to its elementary antiderivative proof. A fully explicit compressed primitive lcm construction proves absence of a uniform denominator bound, as detailed in the proposed mistake entry; the Q=50 instantiation checks every rational a/q for q<=50. No inference of non-loneliness is drawn.

B changes no historical prose or numeric table. It adds an epistemic warning only. Primality of 7,13,19,31 is checked independently; no extra niche Lie statement is introduced or relied upon. No external Lie reference is needed to observe that the source contains no bridge proof; root-system embedding/divisibility claims are marked unsupported here rather than replaced by a new unverified classification assertion.

A preserves all true area/cutoff/gauge statements, geometry, and Mellin formulas. It distinguishes value transformation from density transformation and avoids assuming inversion preserves (0,2]. C preserves the conditional reduction and leaves its Lean proof status untouched: no Lean files were read or kernel verification claimed in this repair packet. Its historical universal-computation language and narrowing are explicitly scoped by the inserted warning.

## Resolved C/D source attribution

The root corrected its initial path attribution: the runner's existing banner belongs to
`07-reflections/the-13-comb-lever-is-the-eisenstein-resonance-types-of-numbers-and-tournaments-opus-20260703.md`, not the endgame reflection. Read that full file and verified its banner explicitly repeats the retracted bounded-q statement. D is the narrow replacement of its final two banner lines with the second-retraction warning, preserving the original body and surviving deep-well correction. `banner-only.patch` isolates that requested repair. The full `proposed.patch` now contains A, B, C, D (four files); C is separately proposed historical warning for the endgame reflection, not the requested banner replacement. Treat C as an optional additional source-bound correction and D as the precise requested banner repair. The manifest now includes all four sources. No root-system classification, Cayley trichotomy, antipode terminology, equidistribution, or other claims in D's historical body were independently certified by this narrow edit.

## Strengthened exact hostile check

For every integer Q>=14, put L=lcm(1,...,Q) and S={1,...,6,L,2L,...,7L}.
Not only does every rational census with denominators <=Q fail: this family
has exact loneliness maximum M(S)=1/8>1/14. At t=(L+1)/(8L), since 8 divides L,
the seven fast phases are j/8 modulo one. The slow phases are
`i(1+1/L)/8`, i=1,...,6, all in [1/8,7/8] because L>=6. Thus the minimum
distance is exactly 1/8. Conversely, among the eight points
0,Lt,...,7Lt modulo one, some circular gap has size <=1/8. The corresponding
pair has index difference j in {1,...,7}, so ||jLt||<=1/8. This bounds the
fast subfamily, hence the whole family's minimum distance, by 1/8 for every t.
If points coincide the same conclusion holds with distance zero.

The displayed witness is already reduced: gcd(L+1,L)=1 and L+1 is odd,
so gcd(L+1,8L)=1 and its denominator is exactly 8L. This denominator is a
witness denominator, not a claim that the minimal lonely denominator is 8L.
The exact argument establishes arbitrarily long finite censuses can miss a
uniformly comfortable lonely family. It neither closes general LRC(14) nor
claims novelty. `check_math.py` verifies the witness and reduced denominator
for Q=14,20,50,100, retaining exhaustive blockage checks through Q=50.

The four proposed source patches are unchanged by this strengthening; the
independent check and unnumbered mistake entry carry the added explanation.
`FREEZE-SHA256.txt` records finalized artifact hashes for publisher review.
