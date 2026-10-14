# Aliquot cycles with retained prime-extension guards

**Status:** PROVED elementary template criterion and order obstructions;
FINITE-EXACT bounded core census and supplied cycle receipts; CITED historical
and analytic context. Existence of an ordinary aliquot cycle of length three
remains OPEN here. No global orbit termination or complete cycle catalogue is
claimed. The elementary compiler is recorded as a reusable instrument, without
a novelty claim for affine aliquot constructions.

## 1. Inheritance and the missing information

Write \(s(n)=\sigma(n)-n\) for positive integers. A cycle in this note is a
nonempty list of **distinct positive integers** with \(s(n_i)=n_{i+1}\), including
the closing edge. Thus perfect numbers have length one; \(1\to0\) is termination,
not a positive cycle. The map is not extended at zero by the production API.

The exact dyadic factor formula is inherited from
[the aliquot/Collatz sidecar](collatz_aliquot_lehmer_five_20260927.md), subject to
[its audit](collatz_sticky_20260927_audit.md): if \(n=2^a m\), \(m\) odd, then
\(s(n)=(2^{a+1}-1)\sigma(m)-2^a m\). Its valuation cases are exact; a population
average or a finite sample of trajectories is a different assertion. The live
coordinates here are the divisor factorization, source, cyclic core word,
positive affine carry, and native prime guards.

The list of lengths \(1,2,5,28\) is not exhaustive. Moews and Moews describe the
historical five- and 28-cycles and exhibit four-cycles. Their reported absence
of three-cycles is a **bounded search result**, in particular for cycles whose
largest member is below \(10^{10}\), not a no-three theorem.
([Primary paper, §§1–2](https://sech.me/ap/S0025-5718-1991-1094955-5.pdf),
[DOI](https://doi.org/10.1090/S0025-5718-1991-1094955-5).)

The source distribution matters. Chum, Guy, Jacobson and Mosunov distinguish
uniformly sampled integers from actual iterates, whose predecessor multiplicities
and persistent divisibility conditions change their weighting. We use no
statistical drift as a source-specific conclusion.
([Primary paper, Introduction](https://arxiv.org/pdf/2110.14136).)

## 2. A complete decision for each fixed core word

Fix an ordered cyclic tuple \(c=(c_0,\ldots,c_{L-1})\), with each \(c_i\ge2\).
The intended sources are

\[
n_i=c_i p_i,\qquad p_i\text{ prime},\qquad \gcd(c_i,p_i)=1.
\tag{1}
\]

The marked prime is a unitary prime factor: it occurs to exponent one. It need
not be the largest prime, and a source may admit several markings. Set
\(S_i=\sigma(c_i)\), \(A_i=S_i-c_i>0\), and use indices modulo \(L\).
Multiplicativity gives the exact transition equation

\[
c_{i+1}p_{i+1}=A_i p_i+S_i,
\qquad f_i(x)=\frac{A_i x+S_i}{c_{i+1}}.
\tag{2}
\]

Compose in chronological order, retaining the unreduced integer carrier:

\[
f_{L-1}\cdots f_0(x)=\frac{Px+B}{Q},\qquad
P=\prod_i A_i,\quad Q=\prod_i c_i,\quad B>0.
\tag{3}
\]

It is computed from \((P,Q,B)=(1,1,0)\) by
\((P,Q,B)\mapsto(A_iP,c_{i+1}Q,A_iB+S_iQ)\).
For an independent expansion, putting \(d_i=c_{i+1}\),

\[
B=\sum_{i=0}^{L-1}S_i
       \prod_{j=i+1}^{L-1}A_j\prod_{j=0}^{i-1}d_j.
\tag{4}
\]

**Fixed-template iff.** There is a primitive aliquot cycle with precisely the
marked template (1) if and only if:

1. \(Q>P\), and \(p_0=B/(Q-P)\) is an integer;
2. successive values from (2) are integers, prime, and coprime to their cores;
3. the resulting nodes \(c_i p_i\) are pairwise distinct.

In that case it is the unique marked cycle for that ordered template.
Indeed closure implies \((Q-P)p_0=B>0\), forcing the first condition and the
unique candidate. The guarded identity (2) proves every actual edge. Conversely,
these conditions close the affine orbit and identify it with the actual aliquot
orbit. Primality can be decided by finite trial division, so this is a total
finite decision without a search bound on the unknown primes. It is not a
polynomial-time claim and does not decide all core words at once.

Two useful consequences are immediate:

- Some core must be deficient: otherwise \(A_i\ge c_i\) for all \(i\), contrary
  to \(P<Q\). This is a statement about the retained cores, not the full nodes.
- Repeating one fixed core cannot give a primitive cycle of length greater than
  one. The repeated positive-slope affine map has only fixed periodic points.
  Its only possible prime is \(\sigma(c)/(2c-\sigma(c))\), when the denominator
  is positive; the native guards still have to be checked.

The individual abundant/deficient threshold is also exact:

\[
s(cp)-cp=p(\sigma(c)-2c)+\sigma(c).
\tag{5}
\]

For deficient \(c\), the source is abundant, perfect or deficient according as
\(p\) is below, equal to or above \(\sigma(c)/(2c-\sigma(c))\). If the core is
perfect or abundant, every allowed prime extension is abundant.

## 3. The false cycle that the native guards remove

For cores \((2,2,16)\), the unique formal candidate is
\((p_0,p_1,p_2)=(23,13,1)\). It gives the formal loop

\[
46\longrightarrow26\longrightarrow16\longrightarrow46.
\]

Its last prime is not prime. In reality \(\sigma(16)=31\) and \(s(16)=15\).
The factor \(\sigma(p)=p+1\) cannot be used at \(p=1\). This is an exact example
where affine closure, positive integers, and a cycle-shaped expression all hold,
but the original map does not close. The compiler rejects it before issuing a
cycle receipt. The separate all-exponent dyadic-core exclusion is owned by
the parent packet; no such global exclusion is inferred from this example.

This is the precise transfer from the guarded Collatz work: compose the complete
affine maps only on their native arithmetic domains. Here the missing condition
is primality/coprimality rather than an exact valuation or a retained sign sheet.
There is no conjugacy between aliquot and Collatz dynamics.

## 4. Order and abundance constraints that survive the quotient

For fixed \(c\ge1,e\ge1\), and a prime \(p\nmid c\),

\[
s(cp^e)=(\sigma(c)-c)p^e+\sigma(c)\sum_{j=0}^{e-1}p^j.
\tag{6}
\]

This is strictly increasing with \(p\), except when \(c=e=1\), when it is
constant one. Thus a nontrivial cycle cannot remain in one such fixed shape.

More generally fix positive exponents \(e_1,\ldots,e_r\) and a core \(c\), with
distinct marked primes coprime to \(c\). The polynomial

\[
s\!\left(c\prod_jp_j^{e_j}\right)
=\sigma(c)\prod_j(1+p_j+\cdots+p_j^{e_j})
 -c\prod_jp_j^{e_j}
\tag{7}
\]

has nonnegative coefficients, and is strictly increasing in each coordinate
outside the same prime-only exception. If **all prime vectors of a putative
cycle form a coordinatewise chain**, distinct vectors give strictly ordered
sources and strictly ordered images. A strictly increasing permutation of a
finite ordered set is the identity, so the cycle is trivial. One comparable pair
alone is insufficient. The actual amicable pair
\(5020=4\cdot5\cdot251\), \(5564=4\cdot13\cdot107\) has crossing vectors
\((5,251)\) and \((13,107)\), so is correctly outside this obstruction.

For every authenticated positive cycle, let \(r_i=s(n_i)/n_i\) and
\(I_i=\sigma(n_i)/n_i=1+r_i\). Then

\[
\prod_i r_i=1,\qquad \sum_i\sigma(n_i)=2\sum_i n_i.
\tag{8}
\]

For a nontrivial cycle, AM–GM and \(1+r\ge2\sqrt r\) give
\(\sum_iI_i>2L\) and \(\prod_iI_i>2^L\). The minimum node is abundant, the maximum
deficient. Also \(\gcd(n_i)\) must be deficient: abundancy is nondecreasing under
divisibility, so a common nondeficient divisor would make every node nondeficient.
These identities constrain a supplied cycle; they do not prove that a trajectory
enters one. For example both 12 and 220 are abundant, but
\(12\to16\to15\to9\to4\to3\to1\to0\), whereas \(220\leftrightarrow284\).

For odd-length nontrivial cycles a stronger gcd restriction is CITED: if
\(d=\gcd_i\sigma(n_i)\), then \(d\mid\gcd_i n_i\), except that when
\(v_2(d)=1\) only \(d/2\mid\gcd_i n_i\) is required. Pollack also proves a
fixed-odd-order counting bound \(x/(\log x)^{1+o(1)}\); this is not an absence
result. We independently replay the gcd restriction on the supplied odd-length cycle.
([Pollack, Lemma 5 and Theorem 1](https://www.pollack-math.net/socodd.pdf).)

No cosmetic tournament is needed. The relevant relation is a partial order on
prime vectors; incomparable pairs are retained rather than arbitrarily oriented.

## 5. Complete core census and its boundary

The universe is all core words with \(2\le c_i\le64\), lengths one through three,
modulo cyclic rotation. Reversal is **not** identified. There is no prime cutoff.
Every core carrier is recomputed by (4) using an independent divisor-sum engine.

| Length | Templates | \(Q\le P\) | Initial nonintegrality | Native guard | Repeated nodes | Cycles |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 63 | 14 | 38 | 6 | 0 | 5 |
| 2 | 2,016 | 218 | 1,684 | 105 | 6 | 3 |
| 3 | 83,391 | 4,651 | 78,169 | 564 | 7 | 0 |

There are 85,470 templates. The eight successful marked templates represent
the six distinct cycles
\((6),(28),(496),(8128),(220,284),(5020,5564)\).
The two markings of 6 and of the 220/284 pair are retained in the JSON.

Thus a primitive three-cycle cannot have at every node a marked unitary prime
factor whose complementary core lies between 2 and 64. Equivalently, any such
cycle would have a node with no admissible marking in that range. Numbers with
no prime factor of exponent one are outside the template; arbitrary larger cores
are also outside. This is a finite-template theorem with unbounded prime
variables, not a new general size bound on three-cycles.

As a longer positive control, the supplied five-cycle
\((12496,14288,15472,14536,14264)\) is reconstructed uniquely from
\((176,304,16,184,8)\). The marked primes are computed, not guessed. Separate
factor-backed controls authenticate lengths one, two, four, five and 28. Their
actual edges are also checked by divisor-pair enumeration.

The next exact obligation is to eliminate or produce native candidates for a
specified unbounded family of changing cores. Merely increasing the finite core
bound would add a larger finite certificate; it would not supply universal
coverage or a termination argument.

## 6. Receipt API and reproduction

`analyze_core(cores)` returns an exact rejection reason or a `CoreCycle`.
`audit_core` rechecks every prime, coprimality condition, marked source, actual
edge and closing equation. `audit_cycle` accepts only distinct positive sources
with canonical prime factorizations. `formal_primes` is explicitly a rational
candidate reader, not a certificate. Exact types are checked before cached
factorization, including bool/float aliases. Production code performs no orbit
search. Known nodes are supplied data whose edges are reauthenticated.

Run from the repository root:

```text
python -B 04-computation/experiments/aliquot_guarded_cycles_20261014e.py
python -B -O 04-computation/experiments/aliquot_guarded_cycles_20261014e.py
```

The matching [JSON](aliquot_guarded_cycles_20261014e.json) stores the universe,
all successful templates, full factorizations for the known controls, rejection
histograms and the formal-cycle hostile. Its output is deterministic after LF
normalization. The independent engine checks all integers 1 through 4,096 and
every supplied/compiled cycle edge; none of these bounded controls is promoted
to an unbounded search claim.

Normal and optimized executions agree with the saved JSON: **98,731 checks**.
LF SHA256: `d9c271d0f872764e8d0d481c64e0f8183d0a6b8e40943649f0b9314802090e49`.
