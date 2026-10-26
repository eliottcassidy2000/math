# Independent adversarial audit: source-fixed shape switching

**Verdict: PASS after two scoped repairs.** The producer's elementary
proofs are valid in their stated domains. The final packet continues to
leave all-prefix root-27 survival and an arbitrary-shape global rank OPEN.
No repository edits or publication were performed by this auditor.

## Reviewed object and integrity

Target: `session66-root27/{REPORT.md,check.py,normal.json,optimized.json,
SOURCE_PINS.json,SHA256SUMS,AUDIT_REPAIR.json}` and the ten pinned source
snapshots. `MANIFEST.json` in this audit records exact reviewed hashes.

All ten source snapshots match their SHA-256, computed Git blob SHA-1,
recorded tree blob, and `git show COMMIT:PATH` bytes in `session66-repo`.
The root-27 anchor and ancestry file are unchanged between `63f52936` and
`2807ee83`. The repaired THM-4771 is correctly pinned to `2807ee83`; it
retains the returned-state valuation sum beginning at index 1.

The original producer normal/optimized receipts replayed exactly. The root
then applied this auditor's proposed repairs, preserving original bytes in
`original-producer/` and recording lineage in producer `AUDIT_REPAIR.json`.
The final producer's normal and optimized runs again replay byte-for-byte,
and every entry of its repaired `SHA256SUMS` verifies.

## Corrections and their scope

1. **Native versus formal blocked-leaf replay.** The original helper
   `literal()` required only `q>0`; it therefore permitted arithmetic
   evaluation of the blocked `q=1,c=0` leaf. Its surrounding comments and
   claimed native macro counts correctly excluded that formal last update,
   so this was an executable guard-separation weakness, not a refutation of
   the refill theorem. The repaired helper rejects the blocked pair before
   mutation. Each family control now executes only `word[:-1]`, checks its
   endpoint has `q=1` and even parity, and separately requires attempted
   blocked execution to raise. There are nine rejected final leaves and
   exactly 837 legally executed family macros. No extra native macro is
   counted or evaluated.

2. **No invented arithmetic Frobenius map.** The last paragraph now says
   the sibling lane has a signed-port permutation and orientation-defect
   specialization, not an identified arithmetic Frobenius action. This
   preserves the useful general permutation-defect analogy without asserting
   a nonexistent arithmetic identification with this affine cocycle.

Both repairs are present in the final reviewed bytes. Original producer
receipt hash: `99db08ce61a5906e63466bfadeead0fa80db6b628b3270d76b6268a45018389f`.
Final producer receipt hash:
`7c405b49acdcbbf895d4db2f8a6b8bd1ff1d96a03a89bcbb8027db3f80ad3c92`.

## Mathematical audit

- **Signed cocycle and cancellation.** Direct coefficient comparison gives
  `B_P D_P N_Q(F_P(z))=A_P D_Q N_P(z)+B_P Delta_(P,Q)` with exactly the
  stated sign. The commutator is `-Delta/(B_P B_Q)`. Both summands have
  valuations `h_P(z)` and `L_P+v_2(Delta)`; division by the left-side factor
  removes `L_P`. Equal valuation gives at least one extra bit, with no
  universal bound on cancellation surplus. Zero defect is its own branch.

- **Carry injectivity and primitive words.** At fixed length, two words
  first differing at position `i` have different parity cylinders modulo
  `2^(i+1)`. Equal carry therefore cannot occur for different words.
  Affine commutation gives equal carries of `PQ` and `QP` at the same
  length, hence literal word commutation. Prefix cancellation/Euclidean
  descent gives exactly the common primitive-word classification.
  No infinite-orbit periodicity theorem is smuggled into this step.

- **Original source and clock.** The source pullback is
  `rho_(w,P)=(2^n r_P-C_w)/3^n`, so
  `z_n-r_P=(3^n/2^n)(s-rho_(w,P))`. This proves (5), with the plus sign
  on `h_P(z_n)+n`. The difference of two pullbacks is
  `-2^n Delta/(3^n D_P D_Q)`, proving (6). Retaining the same center
  through its own block preserves `rho`; switching can refill precision.
  Equal periodic centers do not identify durations or affine maps.
  Equality of clocked affine interfaces on every source forces equal
  lengths and then equal carries, hence identical words.

- **Actual root mapping and hostile.** Independent native iteration gives
  `(7,1),(11,3),(17,5),(26,7),(39,6)` in `(z,q)` coordinates; at cut 3
  this is `(x,a)=(27,6)`. The odd-z and odd-x section clocks are not
  conflated. At macro 86, the independently reconstructed state is
  `(10952205219944621,47)`, followed by `1001010001011101000`.
  After `P=100`, the state is `36963692617313097` with credit 47 and
  fuels `(h_P before,h_P after,h_Q after)=(5,2,9)`. The two cocycle
  summands have valuations `(5,5)` and their sum valuation 12, producing
  the claimed +7 switch gain. The nineteenth old odd-x return hostile
  and the repaired twenty returned-state valuation identities also pass,
  reconstructed by native macro cuts rather than by the producer's
  closed return formula.

- **All-m family and rank quantifiers.** `P Q^m 0` is a full single-tree
  word. Its nonempty residue class yields positive physical-sector
  integer sources. At the state after P, the cylinder for `Q^m0` forces
  `h_Q=6m`. Since this is above defect valuation 2, the cocycle forces
  `h_P=5`. Every proper prefix is legal and the final symbol is blocked.
  Therefore fixed finite shape offsets cannot make `h_P+w(P)`
  nonincreasing universally over positive integer sources. The family
  varies its source; it does not refute a source-27-specific richer rank,
  and the actual macro-86 witness alone does not refute every weighting.

- **Guard and unresolved target.** The root all-prefix condition remains
  the original source-7 tree-cylinder exclusion. The exact switch ledger
  is a carrier, not a well-founded rank. Any global argument still needs
  trajectory coverage, native guards, and control of cancellation gains.

## Expanded independent exact checker

`independent_check.py` imports no producer or inherited research code. It
uses exact integer/Fraction arithmetic, literal parity reconstruction from
all source residues for the small-word census, direct carry sums, brute-force
balanced-word filtering rather than recursive tree grammar, a native guard
that refuses the blocked update, and bit-by-bit lifting for refill seeds. The default execution is offline: it
checks sealed source SHA-256 and computed Git blob hashes without requiring
a checkout. Historical tree verification is preserved separately in
`HISTORICAL_GIT_PIN_VERIFICATION.json` and can be rerun explicitly.

Normal and optimized outputs are byte-identical PASS. Closed universes:

| Control | Exact universe / count |
|---|---:|
| All residue itineraries through word length 12 | 8,191 |
| Balanced words through five internal vertices | 344 |
| Ordered balanced block pairs | 118,336 |
| Common-center ordered pairs | 368 |
| Valuation branches: equal / unequal / zero defect | 17,420 / 100,548 / 368 |
| Pair pullback separation checks at three prefixes | 353,904 |
| Marked source's 111 prefixes against 344 centers | 38,184 |
| Source pins checked by SHA-256 and computed Git blob hash | 10 |
| Unbounded-family finite m controls | 1,2,3,4,5,8,16,32,64 |
| Legally executed family macros / rejected leaves | 837 / 9 |
| Rational balanced centers, three literal periods each | 344 |
| Legal-prefix balanced completions, seeds 1..256, depth <=80 | 7,740 |
| Seeds whose lifetime was resolved by that cap | 163 |
| Seeds not halted by that cap (not nonhalting claims) | 93 |

The all-m statements and the following two lemmas are proved symbolically;
none is inferred from these finite caps.

Reproduction, from this audit directory:

```bash
python3 -B independent_check.py --offline --producer ../session66-root27 > independent-normal.json
python3 -O -B independent_check.py --offline --producer ../session66-root27 > independent-optimized.json
cmp independent-normal.json independent-optimized.json
```

The producer path is configurable so the audit can be relocated alongside a
published packet with no Git repository available. Optional historical
verification is `python3 -B independent_check.py --verify-git-pins
--pin-check-only --producer PATH --repo CHECKOUT`; its separate receipt
records the original committed-byte verification and unchanged anchor/ancestry.

## Separately audited root-proposed lemmas

[`CENTER_BANK_APPENDIX.md`](CENTER_BANK_APPENDIX.md) supplies complete
proofs and exact domains for both additional statements:

- **A:** A finite nonempty bank closed under a nonempty family's inverse
  balanced-block maps must be a singleton common fixed center. Injectivity
  makes each finite inclusion an equality; strict real contraction then
  forces diameter zero. Thus closure under *all* balanced blocks is
  impossible, since `100` and `101000` have different centers.
- **B:** At initial credit `q=1`, the supremum of `v_2(z-r_P)` over all
  nonempty balanced blocks is exactly the number of legal macros before
  the blocked symbol, including lifetime 0 and infinity. Legal prefixes
  have balanced zero-completions, while periodic balanced centers cannot
  share the actual blocked symbol. Finite lifetime is attained by its
  whole legal-prefix block, or by any center in the lifetime-zero case.

The credit restriction in B is necessary: `(z,q)=(4,2)` has lifetime 1,
while all balanced centers have depth 0 at `z=4`. An infinite center bank
therefore recodes the q=1 halting question; it is not a proof of
uncomputability, undecidability, nonhalting, or independently bounded rank.

## Delivery and blocker ledger

Recommended publication allowlist is in `MANIFEST.json`; local original
snapshots, draft patches, and duplicate producer replay receipts need not be
published. Final producer bytes are pinned there, so a later producer edit
requires re-audit rather than inheriting this verdict silently.

**Completed:** proof audit, both repairs verified, final producer replays,
independent expanded normal/-O replays, ten Git-source pins, two appendix
proofs, exact hostile controls. **Denied actions:** none.
**Pending/uncertain access:** none. **Open mathematics:** root-27
all-prefix survival and arbitrary-shape global source-specific rank.
