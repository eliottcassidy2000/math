# Receiver inversion and renewable nine-bank

[Proof and scope](PROOF.md) · [Exact verifier](verify.py) · [Results](RESULTS.json)

From this directory:

```
python3 query.py 6 --universe h3
python3 query.py 1906 --universe h5
python3 query.py 139878 --universe nine-bank-cores
python3 verify.py
python3 -O verify.py
```

The first two query universes are global fixed(h,V) scalar fibers by the
proved/certified low-count core classification. The last is exactly the
stated finite-core universe with unbounded transitive padding, not every
h15 deletion. Each query returns explicit incoming cuts and padding words.

Main outcomes: exact all-order receiver inverse; global h3 no-three theorem;
global h5 inverse; literal renewable NINE-bank and eighteen consecutive odd
counts at infinitely many scales. Chan–Pak completeness remains open.
