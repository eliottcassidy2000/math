# Reproduce the endpoint invoice and scoped searches

From repository root:

```
python3 05-knowledge/results/session290_core_descent_20261113/repairs/verify.py --fullscan
python3 -O 05-knowledge/results/session290_core_descent_20261113/repairs/verify.py --fullscan
```

The input is the inherited complete ordinary order-nine census
`05-knowledge/results/blue_quotient_20261020b/n9_classes.tsv.gz`.
If a sparse checkout omits it, materialize that tracked file or use
`--census /path/to/the/exact/file`; the receipt records its SHA-256.
No new isomorphism census is claimed. Python3 and g++ with C++17 are required.
The script checks the literal endpoint data and every declared slope chart.
Omitting `--fullscan` checks the endpoint proof controls and saved near-match
records only; it does NOT replay finite-search exhaustion.

The separate initial two-port parity probe is reproducible by compiling
`two_port_scan.cpp` with g++ -O3 -std=c++17 and feeding its standard input
the census's `code` column, one decimal integer per line. Stdout followed
by stderr gives TWO_PORT_SCAN.txt. The expanded even-branch scan uses
scan.cpp and is invoked automatically with --fullscan. These are discovery
filters for exactly stated construction universes, not all-order coverage
proofs or claims that finite samples certify exponential identities.
