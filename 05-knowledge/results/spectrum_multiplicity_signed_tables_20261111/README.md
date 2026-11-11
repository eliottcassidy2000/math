# Exact returns, signed table weights, and Hamiltonian multiplicities

Date: 2026-11-11. **PROVED IN NOTES / FINITE-EXACT**, with inherited
classifications explicitly identified. **Chan–Pak remains OPEN.** The target
is realization of every positive odd Hamiltonian path total except 7 and 21.
This packet continues the actual H13/binary-shell construction work and the
user's connections through 99·2^k+23, the first 25 primes, and 2,5,16.

## Positive results and their scopes

| Object | Exact result | Proof / data |
|---|---|---|
| Actual native return | `891·2^k+215 -> 2673·2^k+647`, k>=1; k1 is 1997->5993. Add two siblings; retain every old arc. | [Return proof](returns/PROOF.md) |
| Fixed nontransitive attachments | Exact multipliers 375 and 377 on h+1, for every Q_k; a necessary and sufficient rooted three-coordinate gate. | [One-exit construction](returns/ONE_EXIT.md) |
| Old prime square-sum kernel | Its five paths produce the two H9 cores and a one-arc 23/25 switch with identical E,U,d. | [Kernel switch](signed_arcs/KERNEL_SWITCH.md) |
| Paired binary-shell return | An infinite simultaneous H9/C3 rebase starts at 823/825; an exact two-coordinate recurrence proves closure. | [Kernel switch](signed_arcs/KERNEL_SWITCH.md) |
| Multiplicities | Complete spectra through order 8; exact ordinary, labelled, and cumulative counts for every H<=375 at every order. | [Compiler proof](multiplicity/PROOF.md), [query program](multiplicity/build.py) |
| Actual ready fibres | Every H>1 1-ready class is rigid, uniquely rooted, and a strong ready core followed by a transitive tail. Ready multiplicity stabilizes. | [Rigidity and constructors](multiplicity/RIGIDITY.md) |
| Growth and extension cost | At fixed H=3^a, the proportion with extension-gcd lifetime a tends to 1; exact first corrections differ for ordinary and labelled counting. | [Lifetime theorem](multiplicity/EXTENSION_LIFETIME.md) |
| Signed odd tables | The modulus product gives J(x,sigma)=(3x+sigma,-sigma); 7->21->63 becomes (2,-)->(5,+)->(16,-). | [Weighted atlas](WEIGHTED_ATLAS.md) |
| The two missing arrows | Missing moduli 7 and 21 mean the minus 2->1 and plus 5->16 edges are absent. A finite mixed-sign controller bypasses both on the candidate support. | [Punctured atlas](signed_arcs/PUNCTURED_ATLAS.md) |
| First 25 primes | Exactly 23 positive prime columns; all appear by order 7. Their ordinary sum is 172n-1073 for n>=8. | [Prime columns](multiplicity/prime_columns_lt100.csv) |

No support statement is inferred from a growing multiplicity or from the
signed controller. The full count, endpoint response, vertex order, and
automorphism sidecar are retained where the construction needs them.

## Reproduction

From the repository root:

```bash
python3 -B 05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/multiplicity/build.py
python3 -B 05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/multiplicity/build.py --n 100 --h 39
python3 -B 05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/multiplicity/rigidity.py
python3 -B 05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/signed_arcs/check.py
python3 -B 05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/signed_arcs/switch_check.py
python3 -B 05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/returns/native_return.py
python3 -B 05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/returns/one_exit_probe.py
python3 -B 05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/returns/rooted179_probe.py
python3 -B 05-knowledge/results/spectrum_multiplicity_signed_tables_20261111/seal_packet.py
```

The compiler reuses the pinned inherited ordinary and strong inventories;
it does not enumerate tournaments again. The literal controls test the new
constructors and the exact boundary coordinates. Each proof states the
unbounded conclusion separately from its finite reproduction universe.

## Inheritance and research board

Starting source commit: `63b86513f8be054cb9e2d4e1ab510df8712f5340`.
The first incoming update was `b305e1d03f`, including the all-k triangle
extension filtration and different finite defect-three returns. These are
integrated as proved inputs where cited, with their source scopes preserved.

Closest proved mechanism: native zero continuation and the actual ready
sources Q_k and R_k. Canonical hostile: 1997 has no 1-/2-ready strong
representative, although non-ready marked constructors exist. Corrected
near miss: an unweighted fibre or an odd-table route does not specify the
endpoint response needed to realize its next numerical target. Least-used
sidecars: automorphism orders, retained root ancestry, and the prime graph's
three independent cycle-parity bits.

| Live concept | Operation and preserved predicate | Missing coordinate / next gate |
|---|---|---|
| Anchor: exact target returns | Relocate a newborn predecessor; keep the old induced graph and exact H. | Further nonzero digits need more than the single transitive substituted block. |
| Anchor: paired shell switches | Retain E,U,d while changing H; rebase both members together. | Infinite families of overlaps do not give every address an overlap. |
| Niche: ordinary versus ready multiplicity | SCC decomposition, reciprocal automorphism weights, ready-root inverse. | Abundant unmarked witnesses can coexist with zero ready witnesses. |
| Wildcard: prime square-sum kernel | Suppress five paths, retain parity, then complete native comparisons. | Pairwise cycle overlap loses common ports and XOR dependencies. |
| Niche: signed table atlas | Transport product and multiplicity through m=4x+sigma. | Raw edge descent does not manufacture a tournament of that modulus. |
| Scale: finite control of infinite objects | Strong-order cutoff for fixed H; eight parity states; exact suffix recurrence. | Uniform positivity over all required prime lanes remains open. |

The productive loop was from an old prime graph to its minimal completion
cube, then to the two same-count cores, a native switch, and finally a
simultaneous binary-shell return. The multiplicity loop separates size
growth from usable marked sources and then recovers an injective constructor
from the correct ready subfibre. Both loops yield actual graph operations.

## What the multiplicities now decide

For H=221 and n>=9, the exact ordinary count is `2n^2+20n-262`, but
exactly one class has a 1-ready root. Its ready fraction therefore tends
to zero. For H=63, the ordinary count is exactly `2(n-7)` and the cumulative
count is `(n-7)(n-6)`, for n>=7; both strong cores are rigid. For an
attained product-primitive count, the second backward difference of N_n
recovers the exact number of new strong classes at order n.

Combining the incoming triangle filtration with the two H9 switches gives
a further all-size law. At fixed H=3^a, almost every large-order source
needs a new vertices before any induced extension can remove divisibility
by 3. For a>=2, its exact-lifetime-a proportion is

    ordinary: 1 - 3a(a-1)/n + O_a(n^-2),
    labelled: 1 - 19a(a-1)/n + O_a(n^-2).

The coefficients come from ordered gaps and the automorphism weight of
the two rigid H9 cores. This identifies which growing part of a count
fibre is hard to use for extensions, rather than treating every witness
as interchangeable.

The fixed one-exit multiplier-243 problem reduces to the rooted target
`(H,E,U)=(179,46,278)`. Prime factorization and the exact strong inventory
reduce every possible order to 39 strong classes. None has that response:
the only E=46 roots have U=232 or U=276. Thus this entire specified
uniform attachment route is closed, while the 375/377 routes are positive.
The next construction must change that interface, use another source,
or leave the fixed one-exit grammar. This is not an exclusion of the
requested Hamiltonian totals themselves.
