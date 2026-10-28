# Session 82: the decimal triple, aliquot fibers, and square-sum interfaces

**Status:** proved elementary identities and reductions; exact finite computations with independent implementations; explicitly cited classical inputs; guarded open boundaries. No literature-novelty claim. The session separates a forward orbit, an inverse fiber, a square-sum edge, and a representation change: they are not interchangeable maps.

## 1. What the requested triple actually encodes

Put t=1000 and (A,B,C)=(t²+1,t+1,t(t+1)). Then

- A+B=C+2, not A+B=C;
- B+C=B², giving the square-sum edge1001–1001000;
- the third label1000001 is isolated in the induced three-label square-sum graph.

[DECIMAL_TRIPLE.md](DECIMAL_TRIPLE.md) classifies the extra possible edge for the whole parameter family using a negative Pell equation, including a complete two-orbit descent. It also gives a separate, invertible Euclidean construction of an AP square-sum triangle. A common translation cannot turn this three-label family into a triangle; that positive construction is genuinely a different operation.

The same note identifies the correct primitive abc maps and proves a limited decimal-family radical criterion. The anchor1000001=101·9901 is squarefree and gives no abc radical excess. An explicit congruence family k≡68 mod136 forces17² to divide10^(2k)+1 and supplies ordinary abc hits. It does not prove a fixed positive margin above quality1, refute abc, or identify the original near-additive triple with an abc triple.

## 2. The grouping by aliquot sum is an exact inverse-fiber atlas

Let F_v={n>0:s(n)=v}, where s is the proper-divisor sum. The first rows are:

| target v | complete fiber F_v |
|---:|---|
|0|{1}|
|1|all primes, infinitely many|
|2|empty|
|3|{4}|
|4|{9}|
|5|empty|
|6|{6,25}|
|7|{8}|
|8|{10,49}|
|9|{15}|
|10|{14}|
|11|{21}|
|12|{121}|
|13|{27,35}|
|14|{22,169}|
|15|{16,33}|

For v≥2 these finite cardinalities are [OEIS A048138](https://oeis.org/A048138). The prime fiber at1 is the reason that sequence starts at2. Zero cardinalities are untouchables; singleton fibers are a separate sequence. Targets2 and5 are proved empty; their conjectured uniqueness among prime targets remains open here.

[ALIQUOT_FIBERS.md](ALIQUOT_FIBERS.md) gives complete inverse algorithms based on the largest prime and its exponent, not a naive source cutoff. Two implementations and a third independent audit agree:

| requested target | number of sources | complete description |
|---:|---:|---|
|1001|29|28 distinct-semiprime sources plus11191|
|1000001|5408|5402 distinct-semiprime sources plus six listed exceptions|
|1001000|2|718720 and1000999²=1001998998001|

The semiprime branch is an exact distinct-Goldbach correspondence, with the repeated-prime diagonal expressly excluded. At fixed target it also has the square-deficit coordinate n=h²−r², where h±r must remain prime. This is not yet a square-sum condition.

## 3. A complete even-fiber census, not a truncated sieve

[EVEN_FIBER_CENSUS.md](EVEN_FIBER_CENSUS.md) recovers the classical even-source/odd-square decomposition and computes every even-target fiber through1,001,000. The source universe includes very large odd squares through their smaller square roots.

The500,500 even targets have150,390 empty fibers and642,193 sources altogether. Maximum cardinality18 occurs at exactly three listed targets. A distinct multiplicative-recurrence implementation reproduced every fiber list. This reproduces and extends a finite computation using a cited method; it does not establish the limiting untouchable density or rule out all odd untouchables beyond5.

## 4. Square-sum adjacency inside an aliquot edge has its own exceptional triple

For positive, distinct labels n and s(n), the aliquot edge is a square-sum edge precisely when sigma(n) is a square. [SQUARE_SIGMA_EDGES.md](SQUARE_SIGMA_EDGES.md) uses the explicitly cited Ljunggren square-repunit theorem to classify all prime-power sources:

    3→1,       81→40,       343→57,
    sums4,     121,         400.

These are exactly the prime-power cases, not all integers with square sigma. A self-contained parity-vector argument gives infinitely many coprime squarefree examples in the unrestricted integer class.

Every source in each of the three requested inverse fibers fails this square-sum test. For the two odd targets, the semiprime exclusion has a structural explanation: the corresponding possible Pythagorean hypotenuses501=3·167 and500001=3·166667 have no prime factor congruent1 modulo4. The remaining exceptions are checked exactly.

Touchability and this stricter square-edge touchability have different prime gaps:3,7,11,13 are ordinary touchable targets but have empty square-edge subfibers. The stronger property cannot be substituted into the conjecture about only2 and5.

### Why1001 is singled out by this mechanism

[ANCHOR_CHARACTERIZATION.md](ANCHOR_CHARACTERIZATION.md) proves a sharper statement: **1001 is the unique squarefree product of exactly three primes whose aliquot successor is a prime cube with square divisor sum**. The classical classification forces that successor to343; an elementary factor-pair argument then forces the source primes7,11,13. The full fiber F343 has19 members, so this is a restricted-class characterization, not unrestricted uniqueness. The first edge1001→343 is not square-sum; the following edge343→57 is.

## 5. A genuine translated path, with a recoverable source

[SQUARE_SUM_TRANSLATIONS.md](SQUARE_SUM_TRANSLATIONS.md) enumerates all common integer translations preserving positive labels that produce a square-sum path. There are exactly24. The least offset177612 gives

    178613 — 1178612 — 1177613,
    edge sums1165² and1535².

Both endpoints are prime, uniquely among the24 translations. They therefore share aliquot fiber1, but retaining the endpoint labels and the designated center recovers the original parameter1000 and the offset. The two summaries “both endpoints have aliquot sum1” do not retain that information.

## 6. Why fiber merging creates false paths: exact thresholds13 and37

[FIBER_PORTS.md](FIBER_PORTS.md) connects the arithmetic directly to earlier shared-port work.

The complete fibers F4={9}, F10={14}, F13={27,35} induce two disjoint square-sum edges9–27 and14–35. Merging by aliquot value creates the quotient path4–13–10, but no common middle source realizes it. Maximum label13 is the unique minimum for this precisely defined three-fiber obstruction.

At the next level,33–14–37–12 has both consecutive three-fiber pieces liftable, while the full path has no lift. Maximum label37 is the verified unique minimum for this stronger four-distinct-fiber predicate. These are exact finite minima, not claims that every occurrence of13 and37 has the same explanation.

The all-length repair is to propagate actual reachable source sets, equivalently compose Boolean incidence matrices with their middle-source labels aligned. With repeated fiber labels this decides walks; simple paths need additional used-source memory. Ambient connectivity is not being denied: the first example's endpoints connect through9–27–22–14 when another fiber is allowed.

## 7. Certified forward dynamics and the square/gcd distinction

[ORBIT_CONNECTIONS.md](ORBIT_CONNECTIONS.md) gives the exact terminating short orbits:

    1000001→10003→1437→483→285→195→141→51→21→11→1→0,
    1001→343→57→23→1→0.

The second includes the genuine square-sum edge343→57. Dividing an arbitrary additive aliquot triple by its gcd can instead create a false square signal:483+285=768 is not square, although its primitive third coordinate256 is. The removed gcd must be retained.

The much longer orbit from1001000 has **1144 certified edges**. All5806 primes used by the combined three-orbit certificate were independently checked. The unevaluated endpoint is a64-digit integer with an unresolved62-digit composite cofactor. Its final status, exact edge count, factor hints and reproduction commands are in the orbit note. No termination, cycle-absence beyond the checked prefix, or divergence theorem is claimed. Its generator creates local certificates; an independent integer-only verifier checks every prime certificate and edge. Large certificate bundles are generated, not required as repository inputs.

## 8. The organizing theory

The useful common principle is **source-preserving composition**:

- Factorization coordinates make inverse-fiber enumeration complete.
- A prime-pair gap makes the square-deficit map invertible, subject to primality.
- A labelled center and endpoints recover the translated-square construction.
- A gcd is necessary to transport square tests through primitive normalization.
- A named source at a shared fiber is necessary to compose graph edges.
- Primality certificates and exact integer clocks distinguish a verified orbit from a growth narrative.

These are precise maps with different loss ledgers. They do not assert an isomorphism among aliquot dynamics, tournament grammar, square-sum graphs, or abc. The counterexamples show exactly which purported shortcuts lose an observable coordinate.

## Reproduction and scope

Use the pinned `requirements.txt` for the code needing SymPy. The inverse producer generates deterministic local JSON; dependent checks state their command order. Other finite checkers use only the standard library. Normal and optimized runs are compared, and independent audits are distinguished from cited all-order theorems. Generated JSON, caches and operational receipts are not inputs bundled with the research packet.

The forward-orbit generator is intentionally different: factorization is resource-dependent, and a timeout never creates an uncertified edge. Its own note gives the verified from-seed replay and the exact boundary where this session stopped.

Open questions remain: the global prime-untouchable conjecture, the asymptotic untouchable density, the unresolved continuation of the large aliquot orbit, and stronger arithmetic realization/classification problems. The finite results and conditional statements above do not settle them.

## Suggested reading and exact reproduction order

Start with the decimal identity and inverse-fiber report, then compare the square-sigma and shared-port reports. The code needing SymPy uses the pinned requirements file. Run:

```sh
python3 check_aliquot_fibers.py
python3 independent_prefix_audit.py
python3 even_range_census.py
python3 check_decimal.py
python3 check_square_translations.py
python3 check_square_sigma.py
python3 check_fiber_ports.py
python3 check_anchor.py
```

The first command generates the deterministic input for the square-sigma filter and prefix audit. Each finite checker also passed under `python3 -O`. For the longer orbit replay, use the separate generator/verifier commands in ORBIT_CONNECTIONS.md; the generator requires POSIX alarm support and may take several minutes. The public packet contains nine Markdown files, ten Python programs, and one dependency file, with generated data kept local.

**October28 continuation:** [Eisenstein coordinates and overlapping aliquot properties](../grid_aliquot_overlap_20261028/README.md)
retains the owner's two descriptions of20, proves its unique excess-two
overlap among five-proper-divisor sources, and recovers the full two-prime
family. Its odd prime parameter is always5 modulo8, giving an exact
transfer to the Eisenstein fourth-power obstruction;650 limits the
factorization scope. The1001 characterization above is inherited and
rechecked, with the square-sum and aliquot edge types kept distinct.
