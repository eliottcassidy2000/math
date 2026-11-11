# Historical operations recovered for odd modular tables

Status: **PROVENANCE / TRANSFER MAP**. These notes recover usable operations from distant work; the proofs of the new statements are in [PROOF.md](PROOF.md). Old reflections supply direction, not new theorem dependencies.

The exact central-table map was already recorded in [October 3, tiling wedges and modular atoms, §3](../../tiling_modular_atoms_20261003.md):

> `m=4N-1 (B): [[N,3N-1],[3N-1,N]]`
>
> `m=4N+1 (A): [[3N+1,N],[N,3N+1]]`.

The table uses indices 1 through m−1. This source also retains the wedge addresses, sign reflections, and transpose; a residue alone loses its row/column address. Its odd-index multiplication law
\(a\star b=a+b+2ab\) satisfies \(2(a\star b)+1=(2a+1)(2b+1)\), giving an actual factor witness rather than merely a prime/composite label.

Four older lanes materially change the present operation.

| Source and date | Recovered operation | Map into the present problem | Predicate retained; information still needed |
|---|---|---|---|
| [July 21 CRT fibre-product reflection](../../../../07-reflections/crt-fiber-products-not-modular-cusps-the-exact-clock-carrier-codex-20260721.md), built around THM-2059 | Independent clock phases must lift to one integer; overlap is controlled by reduction modulo gcd | Each marked table contributes an order clock, and the simultaneous lift uses their least common multiple | All marked congruences are simultaneous; individual clock sizes alone would not specify nonzero phases |
| [July 26 THM-2422, operation fibres and twin-center ancestry](../../../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md) | Keep the actual factor fibre and its diagonal, rather than just an existential operation shadow | The new table carries its old modulus as a named proper divisor | Actual compositeness certificate survives; a prime/composite bit would lose the reusable divisor |
| [September 17 signed-cycle arithmetic, §§2,4,7](../../arithmetic_braids2_20260917_signed_cycles.md) | Gcd content, the inverse braid \(4n+b\), and annihilator ideals preserve different coordinates | The new observer period cancels exactly \(\gcd(Q_i,3h+\sigma)\) before taking the order of 4 | The correct factor clock survives; discarding the coefficient would overestimate the period or conflate different arithmetic strata |
| [October 3 route memory, §4](../../collatz_atom_route_memory_20261003.md) | Increase an internal inverse-fibre index by \(3^r\) to preserve r previous exact valuations | Refine that one clock to \(\operatorname{ord}_M(4)\), where M retains all marked observers | Exact valuation head, sign, selected tail, and table factors survive in one infinite affine recurrence |

The June 3 [two-block reflection](../../../../07-reflections/lrc-collatz-the-same-two-block-question-s596.md) supplied an earlier imaginative direction: keep the 2-power, 3-power, and ordered carry together. We use that representation-level pointer, rather than its broader historical analogy, in the exact formula
\[
3^r n+\sigma S=2^A h.
\]
The new operation is not an estimate on the difference of large powers. It chooses a lift exponent whose multiplicative-order constraints make that difference divisible by every required observer modulus.

The distinction between two named signed dynamics was already explicit in [September 17, signed cycles](../../arithmetic_braids2_20260917_signed_cycles.md):
\[
U_{-b}(-n)=-U_b(n),
\]
and in [October 3, paired prime rows](../../paired_prime_collatz_20261003.md). Negation transports forward arrows between signed domains; it does not reverse arrows. The three displayed negative cycles of the plus map correspond to the displayed positive minus cycles based at 1, 5, and 17. The present construction keeps the chosen cycle or tail rather than requiring a classification of all cycles.

The [current modular-center packet](../../modular_centers_collatz_pell_20261111/PROOFS.md) already contains the four mod-8 table classes, affine trajectory guards, row-map information loss, and prime-power decomposition. We do not repeat that work as a new result. Its table decoder, combined with the older head-preserving lift, gives the new infinite-family theorem:

* source: any actual finite signed odd trajectory with a selected tail;
* operation: a synchronized common-target lift at one marked internal state;
* target: an infinite affine family with the same exact valuation head;
* preserved predicate: all specified marked modular observations, hence all supplied factor witnesses;
* discarded information: no original source value is fixed by the lift;
* needed sidecars: sign, ordered valuation head, carry, observer at each mark, and common-tail pointer;
* decisive control: replay a lifted route and check the first smaller increment rejected by the order conditions.

The key positive consequence is stronger than finding more composite examples: **every finite marked trajectory has an infinite common-tail copy family in which every marked table is composite, with an explicit old-table divisor certificate**. The remaining question is how a family library can cover an arbitrary source; that question is not replaced by the existence of arbitrarily many copies.
