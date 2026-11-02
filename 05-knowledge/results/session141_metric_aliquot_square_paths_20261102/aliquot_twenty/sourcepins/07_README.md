# Aliquot fibres, power rows, and the marked 1–2–4 bridge

2026-10-26. **PROVED elementary deductions and carrier classifications;
CITED classical inputs; FINITE-EXACT controls; OPEN global spectrum and
untouchable questions; CONDITIONAL abc consequence.** No claim of historical
priority or resolution of abc, Collatz, RH, or a complete tournament taxonomy.

The strongest connection found here is that **6 is the unique perfect number
receiving a prime square under the aliquot map**. The two arrows
`25 -> 6 -> 6` combine a prime-square predecessor with a fixed point in exactly
one case. The separate observation `6²=4·9` extends to every even perfect
number, with a precise power/doubling address and an aliquot chain through its
Mersenne odd part. These are different mechanisms, and both can be proved.

## Results and their boundaries

1. **PROVED — the unique prime-square/perfect intersection.** Put
   `s(n)=sigma(n)-n` for positive integers. For prime p,
   `s(p²)=N=s(N)` holds exactly for `(p,N)=(5,6)`.
   Euclid–Euler forces `N=2^(r-1)(2^r-1)`. Apart from r=2, r is odd and
   `N-1` is a multiple of three larger than three, so cannot be p.
   The incoming condition itself rules out odd perfect targets; no assumption
   about the existence of odd perfect numbers is used.
   [Proof, cited input, and complete small fibres](aliquot/README.md).

2. **PROVED — all first rows and duplicate addresses.** In the doubling table
   `(r,j) -> 2^j r`, the first row of `n^k`, for `n=2^a u` with u odd,
   is `u^k`. It lies after n exactly when `u^(k-1)>2^a`.
   Hence 6 is the first **even** base with this property for every k>=2;
   3 is first if odd bases are allowed. The power-only and combined tables
   have different first rows. In particular 36 first appears in doubling row9,
   but at combined base3, because `36=4·3²`.
   [Complete address formulas and exact census](chains/RESULT.md).

3. **PROVED — the perfect-number diagram.** If `N=2^(r-1)M` is even perfect,
   with `M=2^r-1` prime, then

   ```text
   N² = 2^(2r-2) M²,       M² > N,
   M² --s--> 2^r --s--> M --s--> 1 --s--> 0,
   N  --s--> N.
   ```

   Thus 6 gives `36=4·9` and `9->4->3->1->0`; 28 gives `784=16·49`
   and `49->8->7->1->0`. The square row M² is not an incoming neighbour of N.
   The prime-square predecessor `(N-1)²` exists only for N=6.
   [Proof and composite-M hostile](chains/RESULT.md#4-the-even-perfect-diagram-and-its-prime-guard).

4. **PROVED — an exact shared counting language with an obstruction.**
   Proper-divisor exponents form a weighted multicomplex. Tournament paths
   are counted by the weighted independence complex of their literal
   odd-cycle conflict graph. The proper-divisor complex itself is the
   independence complex of a graph on its active prime coordinates **iff
   an integer n>=2 has at most two prime factors counted with multiplicity**.
   Beyond that boundary, repeated exponents or a higher-order forbidden face
   must be retained. Graph realizability is still weaker than realization as
   the entire odd-cycle conflict graph of a tournament.
   [Proof and finite controls](tournaments/TOURNAMENTS_124.md).

5. **PROVED — the marked 1–2–4 connection.** On the ordinary Collatz cycle
   only, Collatz acts as multiplication by4 modulo7, the inverse traversal
   of the Paley7 triangle `1->2->4->1`. Under the explicit marking
   `r -> a^r` in `F8=F2[a]/(a³+a+1)`, residue doubling is Frobenius and
   `{1,2,4}` maps to a Fano line. This orbit has rank2, not a normal basis.
   A different primitive generator gives a rank3 normal basis. The marking
   and rank are necessary data; global Collatz does not descend modulo7.
   [Exact maps and counterexamples](tournaments/TOURNAMENTS_124.md#5-a-marked-124-orbit-links-collatz-paley7-and-f8).

## The first handfuls, with complete inverse fibres

The supplied function is the inverse-fibre view of s: target m is paired with
every n satisfying s(n)=m. In particular **9 belongs to target4**, since
`1+3=4`. The target3 has source4. Zero is a legitimate output, while s(0)
is not defined in this packet; a terminal output is not a missing answer.

| Target m | Complete s-preimages | First doubling row of m² | Row later than m? |
|---:|:---|---:|:---|
|0|{1}|outside the positive-base table|—|
|1|all primes|1|no|
|2|empty|1|no|
|3|{4}|9|yes; first overall|
|4|{9}|1|no|
|5|empty|25|yes|
|6|{6,25}|9|yes; first even|
|7|{8}|49|yes|
|8|{10,49}|1|no|
|9|{15}|81|yes|
|10|{14}|25|yes|
|11|{21}|121|yes|
|12|{121}|9|no; its cube has row27|
|13|{27,35}|169|yes|

This table deliberately places the two observations on a common integer axis
without identifying their maps. Powers of two all have doubling root1; prime
squares give `s(p²)=p+1`; distinct semiprimes give `s(pq)=1+p+q`; fixed points
require the different equation `s(n)=n`. Delay alone is common: for every fixed
k>=2 the nondelayed bases through X are only `O_k(X^(1/k))`.
[The generated microcosm](MICROCOSM.json) continues the comparison through32,
including complete fibres and operation-address multiplicities.

The complete inverse bound `s(n)=m>1 => n<=(m-1)²` certifies the table, rather
than assuming a large source cutoff. The bound improves to `(m-1)²/4` outside
the prime-square branch. The recertification through target60 finds exactly
the untouchables `{2,5,52}`. The global assertion that **2 and5 are the only
prime untouchables remains OPEN**. Distinct-prime Goldbach is a sufficient
antecedent for the stronger claim that5 is the only odd untouchable; it is
not used as a theorem here.

## What survives the proposed unification

For ordinary Hamiltonian-path counts H, **7 and21 are excluded at every
order** by [THM-1370 / h-spectrum-omits-7-21-all-n](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md).
They are both aliquot outputs: `s(8)=7`, `s(18)=21`. For the quotient count
`tc=H/|Aut|`, [THM-4790 / presentation-fibres-7-and-21-absent-through-eleven-vertices](../../../01-canon/theorems/THM-4790-presentation-fibres-7-and-21-absent-through-eleven-vertices.md)
has the through-eleven scope; all-order exclusion and ordinary-H completeness
remain open. The forbidden sets express different realization problems.

The identities `7=1+2+4` and `21=1+4+16` retain their projective-plane meaning.
Grouping divisor exponents gives `s(64)=(1+2)(1+4+16)=63`, so the inner
exponent coordinate cannot be dropped. The F1 specialization sets prime
weights to1 and retains the divisor shape, giving `tau(n)-1`; it loses the
arithmetic value s(n). It is an exact polynomial specialization, not a
constructed cohomological theory or an RH implication.

The abc comparison has two exact maps. Dividing
`n+s(n)=sigma(n)` by `g=gcd(n,s(n))` gives a primitive triple, but collapses
every perfect number to `(1,1,2)`. Retaining g restores the lost scale.
Separately `1+(p-1)s(p^k)=p^k` gives an abc-conditional radical lower bound
with its unknown constant. Neither map proves an empty aliquot fibre.
[Quantifiers and losses](chains/RESULT.md#5-two-different-abc-maps-and-exactly-what-they-forget).

The [integration board](CONNECTION_BOARD.md) records sources, targets, maps,
preserved predicates, lost data, sidecars and hostile probes. It connects the
[earlier norm-seven/taxonomy packet](../norm7_taxonomy_zero_20261026/README.md)
to these arithmetic carriers without identifying its normal-basis dictionary
with the rank2 Fano-line orbit. [Corrections](CORRECTIONS.md) repair an
unconditional Goldbach sentence and a minimum-over-a-fibre/every-class error.

## Reproduction and audit

```sh
python3 -B 05-knowledge/results/aliquot_power_124_20261026/verify_packet.py
```

The verifier checks the closed file manifest, pinned inherited sources,
correction preimages, cross-audits, and normal/optimized replay receipts.
The three producer suites contain **296,946 explicit gates**: 19,801 aliquot,
268,651 power/row, and8,494 tournament/F8. A separate additive-sieve audit
adds16,005 checks; repeated replays are not counted as new evidence.
The [audit](AUDIT.md) distinguishes proof review from finite tests. No theorem
ID is reserved or promoted, and no external novelty is claimed.
