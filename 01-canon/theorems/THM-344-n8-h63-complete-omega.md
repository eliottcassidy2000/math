---
id: THM-344
name: n8-h63-complete-omega
status: PROVED (finite exhaustive isomorphism-class census)
date: 2026-05-29
session: opus-2026-05-29-S10
scripts:
  - 04-computation/h63_n8_isoclass_census_s10.py
  - 04-computation/h63_complete_omega_s10.py
results:
  - 05-knowledge/results/h63_n8_isoclass_census_s10.out
  - 05-knowledge/results/h63_complete_omega_s10.out
---

# THM-344: At n=8, H=63 Occurs in Exactly Two Classes, Both with Complete Ω

## Statement

Among tournaments on 8 vertices, exactly two isomorphism classes have
Hamiltonian path count H(T) = 63.

For both classes:

- |Aut(T)| = 1.
- Ω(T) is complete.
- Ω(T) has 31 vertices, with directed odd-cycle length distribution
  `{3: 8, 5: 17, 7: 6}`.
- Hence `I(Ω(T), 2) = I(K31, 2) = 1 + 2·31 = 63`.

Therefore every labeled n=8 tournament with H(T)=63 realizes 63 by the
complete-conflict mechanism, not by a disconnected K3 factor or by any
higher α-vector.

## Exact Census

Using nauty `gentourng 8`, there are 6880 non-isomorphic tournaments on
8 vertices. The S10 census computes H(T) for each representative.

Results:

| H value | n=8 isomorphism classes |
|---|---:|
| 7 | 0 |
| 21 | 0 |
| 63 | 2 |

The two H=63 classes are:

| gentourng cid | score sequence | |Aut| | labeled copies | Ω |
|---:|---|---:|---:|---|
| 2519 | (1,2,2,3,3,5,6,6) | 1 | 40320 | K31 |
| 3285 | (1,1,2,4,4,5,5,6) | 1 | 40320 | K31 |

Thus there are 80640 labeled n=8 tournaments with H=63, a labeled
frequency of `80640 / 2^28 = 0.030041%`.

## Proof

The proof is a finite exhaustive computation over isomorphism classes.

1. `gentourng 8` generates one representative for each of the 6880
   isomorphism classes of tournaments on 8 vertices.
2. For each representative, `h63_n8_isoclass_census_s10.py` computes H(T)
   by Held-Karp dynamic programming.
3. The only representatives with H=63 are cids 2519 and 3285.
4. For those two representatives, the script enumerates all directed odd
   cycles, constructs Ω(T), and checks that every pair of cycles intersects.
   Thus Ω(T)=K31.
5. Since both classes have trivial automorphism group, each contributes
   8! labeled tournaments.

The result is independent of random sampling. The companion script
`h63_complete_omega_s10.py` sampled 100000 labeled n=8 tournaments and
found 33 H=63 hits; all belonged to the same two score profiles and all
had Ω=K31, matching the exact census.

## Consequences

- H=63 is a temporary n≤7 gap, not a permanent forbidden value.
- The old disconnected obstruction `63 = I(K3,2)·I(2K1,2)` is irrelevant
  to the actual n=8 realization.
- The right structural question is now: which complete conflict graphs
  K_r can occur as Ω(T), and what is the minimal n for each r?

## Related

- MISTAKE-024, MISTAKE-050: H=63 was falsely treated as permanently
  forbidden.
- INV-191: H=63 unlocks at n=8 via complete conflict graph.
- HYP-1754: refuted universal H≠63 hypothesis.
- THM-343: H=7 remains impossible for all tournaments.
