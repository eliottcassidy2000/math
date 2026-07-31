---
id: HYP-9070
title: "JC(2): the leading-form Newton circuit as a search gate, and Euclidean depth as the search order"
status: >
  OPEN. Two elementary facts are PROVED + VERIFIED-EXACT here (L0/L1 tower
  identities; the Jung-van der Kulk degree lattice a,b >= 2 coprime). The
  circuit gate itself is a CRITERION/STRATIFICATION proposal, explicitly NOT
  a reduction and NOT a bridge: MISTAKE-237 retracted a previous
  NC2-to-JC(2) "bridge", and nothing here claims an arrow into or out of
  JC(2). Everything asserted is either an exact identity verified on genuine
  Jacobian pairs, or a definition plus a cheap decisive test.
source: death-star-2026-07-31-coinC2
related:
  - THM-3003  # circuit is a complete detector of reversal symmetry
  - THM-3004  # circuit sign changes <= 2K-3, K = #distinct roots
  - THM-3010  # klein: metallic recurrences attain maximal circuit alternation
  - THM-2823  # current JC(2) frontier degree 26
external:
  - "Jung (1942), van der Kulk (1953): the automorphism group of the affine plane."
  - "Lame (1844): Fibonacci pairs maximise Euclidean algorithm length."
scripts:
  - 04-computation/jc2_second_order_tower_hyp9070.py      # L0/L1 derivation + verification
  - 04-computation/jc2_leading_form_circuit_bridge_hyp9070.py  # THM-3003/3004 checks
  - 04-computation/jc2_degree_lattice_hyp9070.py          # depth order + metallic rays
outputs:
  - 05-knowledge/results/jc2_second_order_tower_hyp9070.out
  - 05-knowledge/results/jc2_leading_form_circuit_bridge_hyp9070.out
  - 05-knowledge/results/jc2_degree_lattice_hyp9070.out
---

# HYP-9070 -- a circuit gate and a depth order for the JC(2) search

## 0. Scope discipline first

MISTAKE-237 records that a set of descent programs was wrongly promoted to
an "NC2-to-JC(2) bridge" and to "equivalent formulations of JC(2)". This
file claims **no** bridge, reduction, or equivalence. It contributes (i)
two exact identities about Jacobian pairs, (ii) a lattice on which
counterexample degree pairs must lie, and (iii) two *computable invariants*
of the leading form together with the cheapest tests that would make them
useful or kill them. Anything not verified is marked CONJECTURAL.

## 1. The leading-form tower (PROVED + VERIFIED-EXACT)

Let `(P,Q)` in `C[x,y]` satisfy `Jac(P,Q) = 1`, `n = deg P`, `m = deg Q`,
and let `P_j, Q_j` be homogeneous components.

**(L0)** The degree-`(n+m-2)` part of the Jacobian is `Jac(P_n, Q_m)`, which
must vanish. In **two** variables vanishing Jacobian of two forms means
algebraic dependence, so there is a binary form `H` with

```text
P_n = c H^a,   Q_m = c' H^b,   g := deg H = gcd(n,m),
a = n/g,  b = m/g,  gcd(a,b) = 1.
```

**(L1)** The degree-`(n+m-3)` part gives
`Jac(P_n, Q_{m-1}) + Jac(P_{n-1}, Q_m) = 0`. Using
`Jac(H^a, B) = a H^{a-1} Jac(H,B)` and cancelling `H^{min(a,b)-1}` (wlog
`a >= b`):

```text
c a H^{a-b} Jac(H, Q_{m-1})  =  - c' b Jac(P_{n-1}, H).            (L1')
```

Hence `H^{a-b}` divides `Jac(P_{n-1}, H)` whenever `Jac(H, Q_{m-1}) != 0`.
The exponent pair has moved from `(a,b)` to `(a-b, b)`: **the tower runs the
subtractive Euclidean algorithm on `(a,b)`**, so its depth is the continued
fraction of `a/b`. (L0) and (L1) were verified symbolically on genuine
Jacobian pairs of several degrees, including composites of triangular maps.
*The claim that the higher orders continue to run the Euclidean algorithm
without further degeneration is CONJECTURAL and is the first thing to test.*

## 2. The degree lattice (PROVED, classical input)

By Jung-van der Kulk every polynomial automorphism of the affine plane has
`deg P | deg Q` or `deg Q | deg P`, i.e. `a = 1` or `b = 1`: the tower
terminates at the first step. Therefore

```text
a counterexample requires  a, b >= 2  and  gcd(a,b) = 1,             (D)
```

i.e. a **non-trivial continued fraction**. In every automorphism example
checked, the leading form was a pure power of a single linear form
(`K = 1` distinct root); whether `K = 1` is forced for automorphisms is
open and is a cheap experiment.

**Search order.** (D) makes Euclidean depth `depth(a,b)` (the length of the
continued fraction of `a/b`) the natural complexity of a candidate, not the
degree. By Lame's theorem, among coprime pairs with `max <= N` the depth is
maximised exactly by consecutive **Fibonacci** pairs — the golden ray — and
the constant-partial-quotient rays are the **metallic** ones (Pell for
silver, etc.). Verified: max depth at `b <= 8,13,21,34,55` is attained at
`(5,8), (8,13), (13,21), (21,34), (34,55)`. So a depth-ordered search visits
`(2,3), (2,5), (3,4), (3,5), (2,7), (3,7), (4,5), (5,8), ...` — a genuinely
different order from degree-ordered search, and one in which the golden ray
is the extremal (hardest) direction.

## 3. The circuit gate (CRITERION; the new content)

Dehomogenise the leading form: `h(t) = H(1,t)`, whose roots are the
**directions at infinity**. Then `h` carries the repo's Newton circuit
(`h_k = e_k/binom(g,k)`, `R_k = h_k^2/(h_{k-1}h_{k+1})`,
`c_k = log(R_k/R_{k-1})`), and two canonical invariants attach to any
candidate counterexample:

1. **Reciprocal stratum.** By THM-3003, the circuit is antipalindromic
   (`R_k = R_{g-k}`) **iff** the directions at infinity are reciprocal-closed
   up to scaling (`{r_i} = {mu/r_i}`). This is an `O(g)` test. Its relevance:
   the swap `sigma(x,y) = (y,x)` with `(P,Q) -> (Q(y,x), P(y,x))` sends
   Jacobian pairs to Jacobian pairs (verified: `Jac` is preserved) with the
   degrees exchanged, and acts on `H` by root reciprocation. So the
   reciprocal-closed locus is exactly the locus where the direction data is
   swap-stable.
2. **Alternation gate.** By THM-3004, the circuit's sign-change count is at
   most `2K - 3`, `K` = number of **distinct** directions at infinity, and
   the bound is attained. Contrapositive: a construction forcing alternation
   `>= A` forces `K >= (A+3)/2` distinct directions.

Verified here (corrected convention; an earlier draft of this test used
`R_k = h_k/h_{k-1}` and wrongly reported disagreement): the THM-3003
equivalence holds on reciprocal, scaled-reciprocal and generic multisets,
and the THM-3004 bound holds on its own witnesses.

**Metallic strata sit inside the reciprocal stratum.** For the metallic
number `lambda_q = (q + sqrt(q^2+4))/2` one has `lambda_q * (1/lambda_q) = 1`,
so a metallic root pair is reciprocal-closed with `mu = 1`, hence
antipalindromic. Combined with THM-3010 (metallic recurrences attain
*maximal* circuit alternation), the "maximal alternation" direction data is
a **sub-stratum of the swap-stable locus** — the two extremal notions the
owner asked to connect are nested, not independent.

## 4. Cheapest decisive tests (what would make or break this)

1. **Is `K = 1` forced for automorphisms?** Enumerate automorphisms as
   composites of affine and triangular maps, compute `H` and `K`. If yes,
   `K >= 2` is a genuine counterexample gate and the `2K-3` bound becomes a
   real constraint on the counterexample's circuit. If no, the gate is
   vacuous and this section should be retracted.
2. **Does the higher tower really run the Euclidean algorithm?** Derive the
   degree `n+m-4, n+m-5, ...` conditions and check whether `(a,b)` keeps
   reducing, or whether a degeneration (`Jac(H,Q_{m-1}) = 0`, or a repeated
   root of `H`) blocks it. A block would be more interesting than the
   reduction.
3. **Can the reciprocal stratum be excluded, or is it forced?** Either
   answer is a result; both are cheap to probe on small strata.
4. **Falsifier for the whole file:** if depth-ordering produces no exclusion
   and no prioritisation beyond degree-ordering, and if `K = 1` is not forced
   for automorphisms, then sections 2-3 carry no content beyond bookkeeping
   and should be marked SUPERSEDED.

## 5. Relation to the named repo target

`PROBLEM-LEDGER` section C already names a "golden/worst-approximable degree
corner (Lame-for-polygons)" and Zaremba's conjecture as the target for
JC(2); MISTAKE-237 downgraded the continued-fraction/Newton descent program
that motivated it. This file does **not** revive that program. What it adds
is (i) an exact derivation (L1') of where a Euclidean step actually comes
from, rather than an analogy, and (ii) two circuit invariants of the leading
form that are computable in `O(g)` and are, as far as the survey of this
repo shows, the first contact between the circuit lane (THM-3000-3004,
THM-3010) and the JC(2) lane.
