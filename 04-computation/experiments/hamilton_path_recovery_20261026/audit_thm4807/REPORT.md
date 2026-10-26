# Independent audit of THM-4807

2026-10-26. Separate-agent audit. All audit code is new; it imports nothing from the lane's scripts or from
`../joint_port_bisimulation_20261026`. The audited file is the working-tree version of
`01-canon/theorems/THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation.md`
(203 lines, sha256 `baa0246a…`). The lane session applied all 15 edits in `required_edits.json`, built by
`build_required_edits.py`, before the first push. The agent could not write this file because of a tool hook; the
lane saved the report text verbatim.

## Summary

**Verdict.** Parts (a)–(e) are correct at every order and need only Rédei and Camion. Every lane number reproduced
exactly. The census was extended to all classes at n = 9 and to all strong classes at n = 10. The commentary had one
major error and twelve minor ones.

**Major**

- **M1.** Section (e) said that the traced core alone does not suffice. This is false: `HP(core) = HP(T)` and the core
  is a subgraph of T, so by (d) the core alone determines T. The core is distinct for all 33,867 labelled tournaments
  with n ≤ 6.

**Minor**

- **m1.** Reversing an untraced arc creates one new Hamiltonian cycle for each Hamilton path from u to v, not exactly
  one: 1, 2, 3 or 5 for n ≤ 8, and up to 45 at n = 10. It also creates new 3-cycles.
- **m2.** The figure 6,412,788 is the count of arcs checked at n = 9 alone. The total for n ≤ 9 is 6,589,019.
- **m3.** Case (γ) uses at most two skipping arcs, not exactly two. Both extremal choices are needed: the alternatives
  fail on 3,120 and 3,360 of the 8,160 γ-pairs at n = 5.
- **m4.** The sufficiency proof in (c) left implicit that A meets the components in nondecreasing order.
- **m5.** "The minimum is never above 1" holds at every order, but it rests on Busch's Theorem 4.7 / the BJR
  separator, not on the elementary inputs.
- **m6.** "Exactly the places" overstates. The converse fails: 3,600 of the 3,633 strong 8-classes with a cut vertex
  have no untraced arc.
- **m7.** Missing credit for the third concurrent proof (growing_sieve_path P1–P4, merged in sieve_path_merge). All
  packets already note single-guard separation; what is new here is the explicit word.
- **m8.** THM-1370 is not "driven" by the factorisation; the factorisation is only its step (1).
- **m9.** The description of Busch's Lemma 4.13 was garbled.
- **m10.** The conventions needed |C| ≥ 3 and a one-line proof of the factorisation.
- **m11.** The HYP-9418 link dangled; the lane creates that file.
- **m12.** The maximum number of untraced arcs is known (Busch's thesis Theorem 4.26; Busch–Jacobson, JCMCC 63, 2007).
- **m13.** In the C5 count, an absorbing failure state, if counted, adds one.

## Verdict per claim

| Claim | Verdict |
|---|---|
| Conventions (ordered-join factorisation) | CONFIRMED (m10) |
| (a), both halves including the dual | CONFIRMED |
| (b): chord path, \|C_k\| ≥ 3, no Hamilton path v→u, uP'u is a Hamiltonian cycle of T_D for every D | CONFIRMED |
| (b) "exactly one new cycle" | REFUTED → CORRECTED (m1) |
| (c), necessity and sufficiency | CONFIRMED (m4, m2) |
| (d) cases α/β/γ, the i/j choice, boundaries | CONFIRMED ("at most two", m3) |
| (d) consequence for C5 | CONFIRMED (m13) |
| (e) decoder, cut positions, one direction qualifies | CONFIRMED |
| (e) "neither half suffices alone" | REFUTED → CORRECTED (M1) |
| Endpoint relation collides from n = 4 (44 of 64; 15,343 of 32,768) | CONFIRMED |
| Three decoders agree and are exact on every class n ≤ 8 | CONFIRMED; also all classes n = 9 and strong n = 10 |
| Two-step values through n = 8 | CONFIRMED |
| "The minimum is never above 1" | TRUE, but not elementary → CORRECTED (m5) |
| Busch Prop 4.1, Lemma 4.13, Cor 4.8; no reconstruction statement in the thesis | CONFIRMED (m9) |
| Description of the concurrent packets | accurate but incomplete (m7) |
| "drives THM-1370" | overstated → CORRECTED (m8) |
| "untraced arcs are exactly the places…" | REFUTED → CORRECTED (m6) |
| Untraced-arc counts 1, 2, 11, 39, 200 | CONFIRMED; next term 1,207 at n = 10; the maximum is known (m12) |

## Independent computation (0 failures)

### Method

- Own DFS word enumerator.
- Subset DP for first and last vertices of every T[S].
- Own strong-component routine.
- Tournament classes from `gentourng`.
- DFS and DP agree on every tournament, and the C census matches the Python census at n = 5..8.

### Census of classes

| Quantity | n = 5 | 6 | 7 | 8 | 9 | 10 (strong) |
|---|---|---|---|---|---|---|
| Classes / strong | 12 / 6 | 56 / 35 | 456 / 353 | 6880 / 6008 | 191,536 / 178,133 | — / 9,355,949 |
| Untraced arcs in strong classes | 1 | 2 | 11 | 39 | 200 | 1,207 |
| Strong classes with an untraced arc | 1 | 2 | 9 | 33 | 167 | 1,072 |
| (c) exact, arcs | 60 | 525 | 7,413 | 168,224 | 6,412,788 | 421,017,705 |
| Untraced intra-block pairs (all classes) | 1 | 4 | 18 | 73 | 336 | 1,207 |
| d+(u) − d+(v) at untraced arcs | 2 | 3 | 3–4 | 3–5 | 3–6 | 3–7 |

In every column, all three decoders (E, G1, G4) are exact and agree, and the cut positions equal the component
boundaries.

### Two-step counts (true, false) at untraced pairs

The maximum always lies in the true direction, and the false-direction count never exceeds 1.

| n | Counts |
|---|---|
| 9 | (2,1):30, (3,0):5, (3,1):50, (4,0):9, (4,1):73, (5,0):12, (5,1):84, (6,1):73 |
| 10, strong classes | (2,1):20, (3,0):2, (3,1):78, (4,0):8, (4,1):178, (5,0):16, (5,1):224, (6,0):43, (6,1):212, (7,1):426 |

### Labelled tournaments, n = 1..6

- Strong labelled tournaments: 1, 0, 2, 24, 544, 22,320.
- Hamilton-path sets, trace supports and all three decoders are distinct or exact at every n.

### Separating word, all ordered labelled pairs

| n | Pairs | α | β | γ |
|---|---|---|---|---|
| 2 | 2 | 2 | 0 | 0 |
| 3 | 56 | 50 | 0 | 6 |
| 4 | 4,032 | 3,864 | 0 | 168 |
| 5 | 1,047,552 | 1,039,272 | 120 | 8,160 |

Biased random pairs at n = 6..12 (8,500 pairs) also pass.

### Stress test

Recursively planted Busch Corollary 4.14 structures, n = 9..14, 1,000 cases per n:

- 4,091 untraced arcs;
- 143,891 arcs checked for (c);
- 33,572 flip sets;
- all decoders exact.

## Literature and attribution

### Busch's thesis

- **Proposition 4.1** (p. 39) states that non-traceable arcs join two cut vertices.
- **Lemma 4.13** (pp. 47–48) is stated under the structure of Theorem 4.12. Its condition (ii) has the shape of
  THM-4807 (c), with z in place of u.
- **Corollary 4.8** (p. 43) is quoted correctly.
- The thesis contains no reconstruction statement.
- **Theorem 4.26** (p. 58) bounds the number of non-traceable arcs by (n² − 4n + 3)/8, with equality iff T ≅ T_max.
  The census maxima are 1, 1, 3, 3, 6, 6 for n = 5..10.

### BJR's published statement

Busch–Jacobson, JCMCC 63 (2007), Theorem 1.1, on page 4, attributes the result to BJR (JGT 53, 2006) and matches
thesis Theorem 4.12. The 2006 BJR text itself could not be read (403).

### Concurrent packets

- guard_reconstruction G1 and Session70 G4/G5 are described accurately.
- growing_sieve_path P1–P4 (commit `1b61c5c884`) and sieve_path_merge (`65a93d7747`) were missing from the
  attribution.
- None of these packets deletes the tail u, proves the second half of (a), or uses (b)'s chord path.

### Classical status

- The first half of (a) is Busch's Proposition 4.1.
- "No Hamilton path v→u" is a two-line corollary of BJR.
- The chord path, the tail form (c), the explicit word and the endpoint decoder were not found in the sources read.
- No source found states that Hamilton paths determine a labelled tournament.
- The new sequences (1, 2, 11, 39, 200, 1207; 1, 2, 9, 33, 167, 1072; 1, 2, 8, 44, 614, 15343) are not in OEIS.

## Reproduction

```
python3 -B audit_census.py 8 6
cc -O2 -o audit_census_c audit_census_c.c
gentourng -q 9 | ./audit_census_c 9
gentourng -q -c 10 | ./audit_census_c 10
python3 -B audit_sepword.py 5 1500 4807
python3 -B audit_stress.py 1000 20261026
python3 -B audit_examples.py
python3 -B audit_c_conditions.py
python3 -B build_required_edits.py
```
