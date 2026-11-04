# Independent audit: Euler-aware support sieve and LCM neighborhoods

**Verdict: PASS within the stated scopes.** Audited 2026-11-04, Session 177.

## Mathematical audit

The Euler-aware maximal graph is a necessary supergraph of the true sigma graph. Its order tests correctly handle both v2(order)=0 at non-Euler sources and v2(order)<=1 at the Euler source. The order-one exception is handled with k_r=r, not k_r=1. The construction of simultaneous realizing exponents is valid, including the positive-exponent adjustment when the non-Euler LCM is 1. It realizes only the induced graph on the specified support, not a closed factor network or valuation balance.

The sink-component proof of strong connectivity is correct: sigma(A)|2A, combined with A<sigma(A), forces equality sigma(A)=2A, impossible for a proper full component block of perfect N. The Euler-independent first pass is also correct, since a path stopped on its first Euler visit has only non-Euler sources.

The sample support {3,5,23,31,43}, sole Euler candidate5, has the claimed proper sink {3,23,43}, all indegrees/outdegrees positive, ceiling30659/14784>2 and realizing exponents(10,461,20,1154,2). The numerical sample is explicitly pedagogical and already excluded by stronger known prime-count results.

The quadratic-reciprocity tournament statement for primes3mod4 is correct: at these targets, odd order equals quadratic residuosity, and reciprocity chooses exactly one directed edge.

The LCM closure theorem is correct in both directions. The only degenerate L=1 case is the empty non-Euler neighborhood; choosing an odd prime larger than all finite requirements realizes it. The conjunction with orders5,3,15 at11,13,4561 is exact, and4561 is prime. The implication is useful as an edge constraint but does not itself force both antecedent edges.

No claim proves nonexistence of odd perfect numbers. Novelty is explicitly unestablished. Source-based computational bounds are not imported.

## Independent computational paths

`independent_audit.py` uses modular exponentiation instead of source multiplicative-order enumeration to construct allowed edges; a modular inverse geometric-sum formula instead of binary sum composition to check realization; and BFS instead of transitive Warshall closure. It checks4758 support/Euler cases on subsets of odd primes through43 of sizes2..5,77760 ordered source-target pairs, the unique-sink first-pass equivalence, the sample's rational arithmetic, and direct geometric sums for its reported exponents.

`independent_lcm_audit.py` computes orders by testing divisors of r-1, uses Python's separate math.lcm implementation, checks4032 subsets including closure monotonicity on all one-element extensions, verifies realizing exponents by the modular-inverse geometric formula, and checks9750 actual neighborhoods with exponents through1000.

Both independent scripts use explicit exceptions, so verification remains active under `python3 -O`. Both pass normally and optimized with identical output. The author's two checkers were also inspected after assert hardening and pass normally/optimized; their output is identical in both modes.

Reproduce all four scripts normally and with `python3 -O`. No packages installed, no probabilistic primality, and no hidden computation needed. Source files under audit: `odd_order_support_theorem.md`, `lcm_neighborhood_closure.md`, `check_odd_order_support.py`, `check_lcm_closure.py`.
