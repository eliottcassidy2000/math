# Twenty: exact recovered overlap, typed small atlas, and an arity extension

**PROVED, with inherited attribution:** the five numbers
`20,104,464,650,1952` are all solutions to `sigma(n)=2n+2` with
`tau(n)<=17` (at most sixteen proper divisors). The inherited packet closed
tau<=12; the present proof excludes the remaining exponent shapes.
`130304=2^8*509`, with tau18, is the boundary hostile.

**Recovered correction:** `s(20)=22`, whereas `20+1=21=s(20)-1` excludes
the unit divisor. The user's lists mix prime-power shapes, literal divisor
arity, and a factor-pair compression;16 occurs in both positive and negative
lists. No invented common classification removes that conflict.

**Exact source-preserving links:** `sigma^{-1}(42)={20,26,41}` has an
exhaustive elementary proof. `s^{-1}(16)={12,26}`, `s^{-1}(15)={16,33}`,
and `s^{-1}(22)={20,38}` connect sources with different exponent patterns.
The `1+2+4`/F8/Fano connection is inherited with its marking and rank2 guard,
not presented as a normal-basis or global Collatz theorem.

Read `PROOFS.md`. Portable source snapshots are indexed by original paths,
Git blob IDs, commit, and SHA-256 in `sourcepins/index.json`. They are limited
to relevant consulted sources, not a broad repository copy.

From this directory:

```sh
python3 -B verify.py --replay
```

This checks the closed manifest, independent Git-blob hashes, and two exact
replays (normal and optimized). `check.py` compares independently computed
sigma/tau values at every n=1..200000, direct proper-divisor/compression
checks through5000, exact bounded inverse fibers, controls, tails and marked
F8 operations. All checks survive `python -O`. This finite evidence does not
substitute for the all-size case split. No global excess-two classification,
prime-family infinitude, or cross-field isomorphism is claimed.
