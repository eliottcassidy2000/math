# Digit rules, source memory and highly connected tournament storage

2026-11-03. **PROVED scoped:** DM1–DM6 in [PROOFS.md](PROOFS.md), and TM1–TM3 in
[TOURNAMENT.md](TOURNAMENT.md). **FINITE-EXACT:** three standalone checkers
and saved receipts. Proof review and separate implementations were performed
by the same agent; no independent human or subagent review is claimed.

The main new tournament result extends the incoming S172 memory carrier:
for every r>=1 and M>=2 ceil((r+1)/5), there are 2^M pairwise nonisomorphic
tournaments on 5M+r vertices, of vertex-connectivity exactly r, all having
the same complete path-cover profile. A fixed one-vertex attachment reads
any addressed bit. A transitive feedback block replaces the previous single
articulation vertex. The read survives because its two seed profiles differ
coefficientwise by 2z+14z²+12z³+2z⁴. This is an all-size proof, not an
extrapolation from the literal r<=4 checks.

The arithmetic results identify precisely when compression is safe:

* **Every base and modulus:** split divisibility into a trailing-digit part
  and a coprime periodic-weight part. DM2 gives the exact smallest future-
  suffix observer, including cases with fewer states than the modulus.
* **A sharp operational boundary:** base-p-power appends and legal division
  by p retain exactly v_p(n). Every other base, with the same division
  menu, distinguishes all positive integers. This extends, and sharply
  limits, the inherited decimal-halving decoder.
* **Collatz finite memory:** a k-step parity history together with endpoint
  modulo M has exactly sum_o binom(k,o) M/gcd(M,3^o) possible responses.
  For M=9 this is 2^k+2k+8. These classes have nonuniform source fibres.
* **Source-preserving refinement:** a fixed terminal modulus can erase an
  old phase. Given the word's odd count o, endpoint modulus 3^o M0 is
  minimal for retaining the full source-tail residue modulo M0 when
  3 divides M0; otherwise M0 suffices. These adaptive records are exactly
  the nested source residues modulo 2^k M0. The ordinary clock is k+o.

The two results fit the same precise pattern: a compressed statistic becomes
unsafe when a later permitted operation can test what it discarded. The
appropriate repair can be a named port, a changing modulus, a finite carry,
or the entire source; it is not invariably one extra bit.

The next exploration makes the port cost sharp: q arbitrary added vertices
can remove at most 2q of the storage modules' forced factors of 3, while
preserving old arcs. A single coupled reader attains that budget on two
modules and produces a literal strong tournament with prime H=30431.
On three modules one factor remains. TM2 gives an all-size sharpness family;
this is a concrete escape from the uniform-module divisor, not full spectrum
coverage or a compatible all-depth return grammar.

The same incidence argument gives a second sharp rule: a fixed q-vertex
reader's Hamilton response has Boolean interaction degree at most 2q in
the stored bits. A coupled reader attains degree two and distinguishes all
four two-bit inputs. Matching low-order marginals therefore forces equal
mean responses to small readers, without equating their pointwise answers.

The unchanged storage carriers still have H divisible by 3^M. The coupled
extension escapes that divisor at an explicitly paid boundary. The Collatz controller
retains finite source information but does not prove the actual descent
inequality (2^k-3^o)n>b_w. Universal convergence, full tournament spectrum,
LRC(14), and six-colour forcing remain outside the proved conclusions.

[CONNECTIONS.md](CONNECTIONS.md) records the inherited mechanisms, old niche
sources, source/target contracts, and three next research questions.
[AUDIT.md](AUDIT.md) lists exact universes and limitations.

## Reproduction

From this directory:

```sh
python3 check.py > certificate.json
python3 tournament_check.py > tournament_certificate.json
python3 audit_independent.py > audit_certificate.json
python3 -O check.py > certificate.optimized.json
python3 -O tournament_check.py > tournament_certificate.optimized.json
python3 -O audit_independent.py > audit_certificate.optimized.json
cmp certificate.json certificate.optimized.json
cmp tournament_certificate.json tournament_certificate.optimized.json
cmp audit_certificate.json audit_certificate.optimized.json
```

All scripts use explicit failures, exact integers/rationals, and no imported
repository implementation. No stochastic model or probabilistic sampling is
used. The recorded arithmetic source universes and graph deletion universes
are exhaustive within their stated finite ranges, not universal experiments.
