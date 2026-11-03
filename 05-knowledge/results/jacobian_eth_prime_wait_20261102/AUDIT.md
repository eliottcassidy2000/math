# Audit and reproduction ledger

2026-11-02. **Same-agent proof review, with independent computational routes.
No independent subagent review and no external Lean/kernel replay.**

## Mathematical gates

* The ambient dimension m, generic covering degree n, polynomial coordinate
  degrees and partition function p_part are distinct. W3 concerns every
  m>=3 and n>=4; JC(2) is untouched.
* W1 assumes C!=0 only in its simple-root formula. Its separate C=0 proof
  exhausts the disjoint source loci x=0 and gamma=0,x!=0. The quartic
  three-point boundary is a positive hostile against forgetting this split.
* Repeated nonzero R roots are the new lane. The designated root 1 is
  always simple, so it cannot be included in an escaping cluster. Generic
  inertia is computed away from B=0 and other residual collisions; the
  all-target fibre formulas do not have that generic qualification.
* The discriminant is computed in the inverse primitive, not in the source
  Jacobian. Its index contribution is two in the valuation formula; this
  is not a finite critical locus. Every checked source determinant is a
  nonzero constant.
* Inertia is one cyclic meridian action with several disjoint cycles. It
  is not a direct product of independently chosen branch generators.
* The invariant used for inequivalence is an unordered **multiset** of
  component data. It survives component swapping; finite-sheet counts
  distinguish identities that total permutation cycle types alone miss.
* W3 counts explicitly constructed inequivalent maps. Exhausting integer
  partitions does not exhaust arbitrary Keller maps or prove finiteness
  of the classes in one degree. Global S_n and geometric degree remain
  inadequate complete invariants.
* A prime's integer size, original prime rank, coprime-prime rank and rank
  after m-2 are separate clocks. The 290 progression certifies a lower
  bound, never primality of its first uncovered candidate. CRT proves
  unbounded p_*/m, not an unbounded normalized eligible-rank ratio.
* The Stern comparison has a real counterexample at 557: preserving ten
  failures does not preserve failure over a growing quadratic range.
* The MRB filter keeps parity explicitly. Its mean is exact and supplies
  the logarithmic-square coefficient. The orientation repair is a
  convergent new functional, not a new formula for the MRB constant.

## Exact checks

```sh
python3 -B 05-knowledge/results/jacobian_eth_prime_wait_20261102/prime_wait.py > /tmp/prime-wait-normal.json
python3 -O -B 05-knowledge/results/jacobian_eth_prime_wait_20261102/prime_wait.py > /tmp/prime-wait-opt.json
cmp /tmp/prime-wait-normal.json /tmp/prime-wait-opt.json
python3 -B 05-knowledge/results/jacobian_eth_prime_wait_20261102/weighted_fibres.py > /tmp/weighted-fibres-normal.json
python3 -O -B 05-knowledge/results/jacobian_eth_prime_wait_20261102/weighted_fibres.py > /tmp/weighted-fibres-opt.json
cmp /tmp/weighted-fibres-normal.json /tmp/weighted-fibres-opt.json
git diff --check
python3 agents/check_docs.py
```

The stored JSON files are the corresponding normal outputs. Prime checks
use sieve counts versus trial division. Fibre checks use symbolic inverse
and boundary identities versus direct finite-source enumeration. All 7,196
source triples over the three declared fields were evaluated, including
the newly admitted repeated-root seeds. Discriminants are computed with
both A and B indeterminate, avoiding an accidental special slice.

Partition enumeration is checked against an independently implemented
coin-change recurrence through degree 14. These finite checks support
the general proof; they are not its quantifier justification. MRB's h_K
shift is checked as a generating-function identity for K=1,...,8, with
the all-K argument inherited from the contour audit.

## Source and history boundaries

The numerical 290/227 search recovered older results in metagraphs,
Keller discriminants, Berggren local screens, Stern primes and filtered
prime sums. Hashes and identifier matches were not interpreted as
mathematics. The exact shared pair (290 edges,56 classes) versus
phi(290)/2=56 was tested against regularity and does not produce a Cayley
graph identification.

External ETH source statements and comparator instructions were read.
The two published time scopes are not merged. Statement-source hashes
authenticate what was read, not the validity of an unexecuted proof build.
The formalization's expected axioms are not reported as a local audit.

## Handoff

The new all-target weighted-family rule and partition bound are proved
in the companion note, with local IDs W1–W4. W4 keeps the entire positive
coordinate fixed and W1 gives a two-gcd fibre-count algorithm. Incoming
Session144's stronger MRB obstruction and corrected equivariant scope
are inherited explicitly. No scarce THM namespace was
reserved and no false universal classification was promoted. The next
stronger target is an invariant detecting moduli inside one fixed
component passport, or a polynomial-gluing criterion for a proposed
cover. For prime arrivals the next useful change is a specified sampling
measure and retained CRT state, not another interpretation of phi/2 as
an unconditional expected first rank.
