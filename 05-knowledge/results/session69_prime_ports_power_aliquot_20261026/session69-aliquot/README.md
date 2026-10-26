# Session69 aliquot/prime/ABC packet

Read `REPORT.md` first. `check.py` independently checks a copied, pinned Session66 inverse algorithm against an ordinary divisor sieve, then tests the new elementary branch identities and hostile boundaries.

```
python3 check.py > replay.json
python3 -O check.py > replay-optimized.json
cmp replay.json normal.json
cmp replay-optimized.json optimized.json
python3 check.py --negative-control   # expected exit 1: missing source 25
```

Finite scope: complete fibers for targets 2–1000, prime-target witnesses through 1,000,000, guarded core/prime branches for cores 2–200 and primes ≤199, twin-prime collision witnesses for p≤997, and all labelled tournaments on four vertices. Normal and optimized receipts match. The source checkout was read-only; there are no repository edits or publication actions in this packet.

Load-bearing mathematical boundaries: prime-only untouchable uniqueness is not proved here; distinct Goldbach is sufficient but not asserted equivalent; fixed-core radical constants are not uniform in the core; raw H and orbit-normalized tc are distinct.
