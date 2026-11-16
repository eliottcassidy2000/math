# Chan–Pak: where the gluing (J) can live, one shared vertex, and {7, 21} as conflict graphs

Session klein-2026-11-16-join (2026-11-16). The owner's prompt:
- 21 is h(10);
- search for the gluing (J), p + q − 1 for primes ≥ 11;
- the 24-vertex two-dial bank and the literal triangle shield;
- an analogy with {7, 21}.

The same prompt was worked concurrently by wave 43 of
[chanpak_waves_20261113](../chanpak_waves_20261113/README.md#wave-43-the-carry-odometer-tested-an-exact-triangle-shield-two-apex-gluing-rigidity-and-7-21)
and by [fresh_unit_prime_gluing_20261116](../fresh_unit_prime_gluing_20261116/README.md). They prove, and
this packet does not repeat:
- two-apex and transitive-overlap gluing gives only 3 + 3 − 1;
- no overlap of three or fewer vertices works for primes ≥ 11;
- (11, 11) fails, and exact deletion-form families exist;
- the reading of {7, 21} as {h(3), h(3 ∘ 1)}.

This packet adds the following.

| result | status | file |
|---|---|---|
| **Target map.** For every target N ≤ 123 and each row p + q − 1, p + q + 1, p + q + 3 (primes ≠ 7): which prime pairs can sit as induced subtournaments inside any tournament of count N. Complete, because the strong component carrying both cores has count N and order ≤ 9 (Busch floor). | FINITE-EXACT + CITED (Busch 2006); independent OCF replay for N < 45 | [TARGET_MAP.md](TARGET_MAP.md) |
| In the (J) row only 18 of 51 targets admit any pair. 23, 27, 29, 31 and 29 others admit none. 11 never sits as a (J) leg below 125. | FINITE-EXACT | [TARGET_MAP.md](TARGET_MAP.md) |
| Over all three rows, **27, 57, 67, 71, 91, 103, 105, 115** admit no induced pair. So CONDITIONAL.md's gluing must rebuild ("prepare") its witnesses at these counts. | FINITE-EXACT | [TARGET_MAP.md](TARGET_MAP.md) |
| **One shared vertex.** Prime witnesses (not both 3) meeting in at most one vertex force H ≥ p + q + 7 at every order. This sharpens CONDITIONAL.md's p + q + 3 and kills all three of its rows at one-vertex overlap. | PROVED; exhaustive replay at orders ≤ 9 | [ONE_VERTEX.md](ONE_VERTEX.md) |
| **HYP-9862.** H(A ∪ B) ≥ H(A)H(B) whenever \|A ∩ B\| ≤ 1. No failure among 218,925,648 order-9 instances; it fails at two shared vertices (5 < 9). | OPEN; FINITE-EXACT at orders ≤ 9 | [ONE_VERTEX.md](ONE_VERTEX.md) |
| **Correction (MISTAKE-775).** CONDITIONAL.md's quantitative overlap bound drops the factor r = H(W) in the intersection term. Counterexample at order 9: H = 99 < 107. The corrected bound holds over all 1.79·10⁹ order-9 covers. | PROVED; brute-force check | [ONE_VERTEX.md](ONE_VERTEX.md) |
| **All conflict graphs of value 7 and 21.** 7 = K_3 = 3 ∨ 5 is the only form of 7. 21 has six forms: the shield of K_3, P_4, and four joins. Every join form splits as 3 ∨ 19 and 5 ∨ 17; three of them (K_10 = h(10), K_8 − e, the matching form) also split as **11 ∨ 11**. | PROVED (finite enumeration) | [SEVEN_TWENTYONE.md](SEVEN_TWENTYONE.md) |
| **Two dials, one move.** The +6 dial is the +2 dial under one triangle shield. Over the empty base the two-dial rectangle is (K_j ⊔ K_1) ∨ K_i, values 3 + 2i + 6j, and its forbidden cells are exactly where K_3 sits at the top level or under the shield. "Exactly once" holds because the disjoint channel is the only exact ×3. | PROVED identities; ANALOGY reading | [SEVEN_TWENTYONE.md](SEVEN_TWENTYONE.md) |
| **Binary Goldbach is not needed.** A physical join of at most five prime legs, plus Helfgott's ternary theorem, implies Chan–Pak. A leg 7 splits into the dials 3 and 5. | PROVED reduction; the physical hypothesis is OPEN | [TERNARY.md](TERNARY.md) |

**Not claimed.** Chan–Pak (OPEN); any physical gluing rule; HYP-9862.

## Reproduce

```
cc -O3 jglue.c -o jglue; for n in 4 5 6 7 8 9; do gentourng -q $n | ./jglue $n; done   # cover census, CONDITIONAL bound check (~40 s)
cc -O3 jcontain.c -o jcontain; for n in 3 4 5 6 7 8 9; do gentourng -q $n | ROWS=1 ./jcontain $n; done | sort -u   # rows_census.out
for n in 5 6 7 8 9; do gentourng -q $n | ./jcontain $n; done   # jrow_forms.out (forms C/K/D)
python3 make_tables.py         # target_map_table.out
python3 verify_rows_ocf.py     # independent OCF replay, N < 45
cc -O3 onevx.c -o onevx; cc -O3 supmult.c -o supmult   # onevx.out, supmult.out
python3 omega_forms.py         # omega_forms.out (needs networkx)
```
