# Exact synchronization language for four physical carrier pairs

**PROVED, independent audit PASS.** For each of the four
Session203 paired carriers at bases 677, 759, 1547 and 821, the common
O/Z/J words preserving equal adjusted bases are exactly the O/Z words.
Every first J produces a permanent signed separation; no later common
suffix repairs it. This classifies the whole word language, not just a
one-step failure or words through a finite cutoff.

The result does not exclude asynchronous words, source replacement, new
carriers or different constructions. It does not prove eventual odd-spectrum
coverage. Its value is to identify exactly why adding common J macros to
these existing constructive families cannot supply a missing induction.

The independent mathematical review required no corrections. Its separate
direct-word implementation recounted all 104 literal profiles and checked
1,062,880 common-word cases through length eleven, with exactly 4,095 O/Z
words synchronizing per pair. The proof's earlier candidate-status line is
retained to preserve the exact audited bytes; this README records the
completed review. The release's `PAIRED_LANGUAGE_AUDIT.md` supplies the report.

- `PROOF.md`: all-order signed-cone proof and precise scope.
- `verify.py`: standalone Python verifier with explicit error gates.
- `RESULT.json`: deterministic normal/optimized receipt.
- `dependencies/`: prior audited physical constructions and transfer laws,
  verified unchanged on main at `4c63fd58b566f3e006eb87eaedcf48824741cc84`.

Run from this directory:

    python3 verify.py > /tmp/paired-language-normal.json
    python3 -O verify.py > /tmp/paired-language-optimized.json
    cmp RESULT.json /tmp/paired-language-normal.json
    cmp RESULT.json /tmp/paired-language-optimized.json

The verifier checks 54 exact linear-basis identities, 104 literal graph
profiles on the eight sources and their children through order ten, and
354,292 paired-word cases through length ten. The word cases are finite
controls; the unbounded theorem is proved in the text. The small census
experiments outside this packet are exploratory and are not dependencies.
