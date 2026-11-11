# Reproduce Session274 structural controls

Python3 standard library only. From this packet directory:

```
python3 constructions/check.py > /tmp/session274-constructions.json
python3 geometry/check.py > /tmp/session274-geometry.json
python3 moments/check.py > /tmp/session274-moments.json
python3 review/independent.py > /tmp/session274-review.json
```

Compare the JSON objects with each lane's committed results (key ordering is
not mathematical). The exact declared universes and proof-versus-test
boundaries are in the lane proofs and independent review. No large coverage
range is being tested. No external mathematical package is needed.

The moments script additionally compares its explicit run-sum with the
three-state recurrence at every n=1..15. Geometry distinguishes framed tiles,
cut sectors and complement lines from ordinary tournament isomorphism.

`review/REVIEWED_SHA256SUMS` pins reviewed proofs. Source code hashes in
`SHA256SUMS` identify the final delivered packet. Published results depend on
the proofs for all orders, not on extrapolating these small tests.
