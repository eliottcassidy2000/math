# A sharp two-port decoder for the leaf atom grammar

The Sessions 14/15 singleton/C3/S4 leaf grammar and Session 27's marked response
interface fit together in an exact binary experiment. In a known dominance
join with s independently reversed C3/C5 factors, all one-port count tables
agree, but exactly s two-port Hamilton-count queries modulo 4 recover all
s labelled orientation bits. Over EVERY labelled paired-black leaf, the full
one-port fibre is exactly the independent triangle reversals; its SCC
partition and S4 center are already determined. The matching lower bound applies to this
binary oracle, not to exact-integer queries or arbitrary reconstruction.

- [Proof and freshness ledger](PROOF.md): general module hiding/visibility
  is inherited from later repository work; the scoped decoder is proved here.
- [Standalone exact checker](verify.py), [results](RESULT.json), and
  [normal output](verify.out).
- [Pinned primary sources](SOURCES.json).

Run `python3 verify.py` or `python3 -O verify.py` with standard-library Python.
Each run writes RESULT.json beside the script. The normal and optimized
outputs agree. The independent audit, if included by publication, has its
own explicit scope.

The 2^s variants are labelled presentations of one unmarked tournament,
not 2^s metagraph nodes. The decoder does not reveal complement wiring,
construct new Hamilton-count values, or solve the eventual odd threshold.
