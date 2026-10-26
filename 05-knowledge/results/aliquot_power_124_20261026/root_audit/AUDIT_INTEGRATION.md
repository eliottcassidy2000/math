# Independent integration audit

**Verdict: PASS.** This is a read-only mathematical and data audit of the
integrated `aliquot_power_124_20261026` packet. No remaining substantive defect
was found in the reviewed statements. The reviewer did not modify the producer
notes or checkers. Final publication, the closed packet manifest, and its full
replay are the integrating author's responsibility.

## Reviewed boundary and frozen hashes

Paths in this table are relative to the packet root. The files were read from
`/tmp/math-wt-aliquot_power_124_20261026/05-knowledge/results/aliquot_power_124_20261026/`.

| File | SHA-256 |
|---|---|
| `README.md` | `c0dff724182e51332a2116014b05648e88a97216069e45e9be8237847b4874a0` |
| `CONNECTION_BOARD.md` | `9e9397586635848717feb9fe251ecb76244cf8febc2c4952370375dcb0612991` |
| `CORRECTIONS.md` | `51bb5b0eab0156fed88edd8958a7da8a948447bdb74dda86c2e2cd789cd6f171` |
| `CORRECTION_PREIMAGES.json` | `519d5baa7c7f88544c673316fabb1869f08dc9e6fda9b4bf70452d8c41ee3b15` |
| `build_microcosm.py` | `b5c26e709a16ce826aed17c9876b5c96f69042e6031b52e9b9c58e229e88f04d` |
| `MICROCOSM.json` | `ad7a9e03c7f9db707c4c087cfb358cfd7a8a96d01a1f929f29dcb947435a25f1` |
| `verify_packet.py` | `ef017a8f93f93f95ed70332e9db3255fbbb8ab049271a12478d88038ff63e1c2` |

The review also used the child proof notes and their existing independent
audits. It does not replace those all-parameter proof audits with a finite
calculation.

## Mathematical checks

### The triangle, circle, and box have a precise common carrier

For squarefree `n=p_1...p_r`, the proper-divisor exponent vectors are
`{0,1}^r` with the all-ones vector removed. Taking a vector's support gives
exactly the proper subsets of the `r` prime coordinates. As an abstract
simplicial complex this is the boundary of an `(r-1)`-simplex. For `r=3`
its realization is a circle, with three vertices and three edges. For `n=30`
the seven faces, including the empty face, have degree profile `(1,3,3,0)`;
their weights are the proper divisors `1,2,3,5,6,10,15`, summing to `42`.

This is the monomial-support face construction. It is not the geometric cube
boundary. The order complex of the proper-divisor poset including `1` is a
cone because `1` is a least element. For squarefree `n`, first removing both
`1` and `n` and then taking the order complex gives the barycentric subdivision
of the simplex boundary. The final board explicitly distinguishes these maps.
Prime weights and repeated-prime exponents are necessary sidecars: incidence
alone does not determine an aliquot sum, and repeated-prime examples need a
multicomplex.

The graph independence-complex criterion has the necessary domain `n>=2`.
For a prime there are no active vertices; for a prime square there is one;
for two distinct primes the complex is `Ind(K_2)`. If the total number of
prime factors counted with multiplicity is at least three, either a repeated
coordinate gives a surviving nonsquarefree monomial or the squarefree case
contains every pair while omitting a larger face. These respectively obstruct
being a simplicial complex in the stated coordinates and being flag. The
criterion does not assert tournament realization of an arbitrary graph.

### Six, twenty-eight, and operation addresses

The first unrestricted base with a delayed square is `3`, whose square is
`9`; the first even base is `6`, whose square is `36`. Under doubling alone,
`36` has first row `9`. Under the combined power/doubling carrier it has first
base `3`. The final text keeps these operations and orders distinct.

The perfect-number statements also concern distinct maps. The prime-square
incoming edge is `25 -> 6`, and `6` is fixed by the aliquot map. Squaring `6`
instead gives odd-core row `9`, with tail `9 -> 4 -> 3 -> 1 -> 0`. For `28`,
the square `784` has odd-core row `49`, with tail `49 -> 8 -> 7 -> 1 -> 0`;
there is no prime-square incoming edge to `28` because `27` is composite.
The unique-six proof uses the even-perfect classification and the mod-three
obstruction, with the prime-two case treated separately. It needs no
assumption that odd perfect numbers do not exist.

Delayed rows have density one under the stated fixed-power counting problem.
They therefore do not characterize perfect numbers. The formal composite
Mersenne candidate `N=120` supplies the necessary hostile control to the
prime-guarded aliquot tail.

### Quotients and open boundaries

The primitive normalization of `(n,s(n),sigma(n))` sends every perfect number
to `(1,1,2)` and loses the common scale. The separate prime-power identity
`1+(p-1)s(p^k)=p^k` has coprime terms and gives the stated radical inequality
conditionally on abc. Neither map supplies a missing-fibre theorem or an
effective universal cutoff. The Goldbach coverage antecedent remains
conditional; finite untouchable and fibre censuses remain finite. The claim
that `2,5` are the only prime untouchables is not promoted to a theorem.

Tournament `H` and `tc` gaps retain their distinct scopes. The explicit
`H=127`, minimum-feedback-`5` tournament invalidates identifying a minimum
over a fixed-`H` fibre with every member's minimum. The correction records
retain that mechanism. The ordinary Collatz multiplier description is confined
to the cycle `{1,2,4}`; the congruent inputs `2` and `9` disprove extension to
a global quotient map modulo seven. Evaluation of divisor variables at `1`
is an exact polynomial specialization, with no claimed RH or motivic theorem.

### Incoming marked Frobenius law

The imported portion of `rooted_lines_frobenius_sheets_20261026/ARITHMETIC.md`
is the elementary A1 orbit law at commit `563022f893`. In a chosen normal
basis, a binary word of period `d` has Frobenius orbit length `d`. Marking
one normal-basis vector forces the pair's orbit to have length `m`. More
generally the diagonal action on cycles of lengths `d,e` has `gcd(d,e)`
orbits, each of length `lcm(d,e)`.

The essential marking condition is full field degree, not normality. Each
element in the rank-two orbit `{a,a^2,a^4}` in `F_8` still has field degree
three, so individually marking it kills the Frobenius stabilizer. That orbit
cannot serve as a normal basis for the binary-word coordinate dictionary.
The final board states this distinction correctly. No global Collatz result
or reserved theorem is imported.

## Independent exact checks

A separate inline standard-library calculation, with no producer imports,
passed **351 gates**. It covered all 32 microcosm rows, independently computed
proper-divisor sums by paired divisors for sources through `961`, and checked
complete inverse fibres through target `32` against the inherited finite
bound. It checked the operation-address, first-power-base, odd-core,
prime-square-preimage and perfectness fields, including the typed `n=1`
exceptions. Combined addresses were derived from valuations; first power
bases were obtained by direct base/power iteration.

The same calculation checked the seven weighted faces for `30`, and
enumerated all `7!` vertex orders for the tournament with row masks
`(106,60,89,48,97,68,10)`. It independently obtained `127` Hamiltonian paths
and minimum feedback number `5`.

An additional **100 gates** enumerated diagonal rotation on every product of
cycles with `1<=d,e<=10`, obtaining exactly `gcd(d,e)` orbits of length
`lcm(d,e)`. These are finite controls; the orbit proof is the elementary
stabilizer calculation, not extrapolation from the census.

Both normal and optimized runs of the microcosm producer were byte-identical
to the frozen `MICROCOSM.json` listed above. Reproduction commands from the
packet directory are:

```bash
python3 -B build_microcosm.py
python3 -B -O build_microcosm.py
```

For each entry in `CORRECTION_PREIMAGES.json`, the reviewer independently
verified the current postimage hash, the uniqueness of the replacement text,
and exact reconstruction of the declared preimage hash. This establishes
the advertised correction lineage without treating earlier text as current
truth.

## Review amendments and closure

The integrating author corrected the following issues raised during review:

1. The first-occurrence comparison now says **base 3 (square 9)** precedes
   base 6; it does not compare the output `9` directly with the base `6`.
2. The graph-complex criterion explicitly excludes `n=1`.
3. The packet verifier accepts the chains receipt's exact `FINITE-EXACT`
   status and uses the correct independent tournament receipt filename.
   The frozen chains producer was not altered to satisfy the verifier.
4. The shape paragraph distinguishes monomial faces, cube boundary, and
   divisor-poset order complexes.
5. The Frobenius paragraph distinguishes full field degree, conjugate linear
   rank, and the normal-basis coordinate requirement.

The final verifier was inspected for these amendments, including its
source-lineage checks, explicit runtime gates, and normal/optimized replay
logic. This audit does not assert that the final closed manifest replay has
already run: the author is assembling that manifest after receiving audits.
The PASS verdict applies to the exact files hashed above and the stated
mathematical/data boundary.
