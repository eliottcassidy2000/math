# Tournament offset moves as fixed-port tails

**Session:** kind-pasteur-2026-11-09-S19. **Canon:** THM-5085, HYP-9655.

**Status.**
- PROVED: the fixed-port transfer, its closed form and signed offset transport (THM-5085 (1)-(3)).
- FINITE-EXACT: coverage to 10^6 in bases 2 and 3 (pre-registered Q1, Q2 CONFIRMED; Q3 REFUTED as worded, the intended
  version CONFIRMED; Q4 UNDECIDED).
- OPEN: tournament (binary) Chan-Pak and ternary Chan-Pak.

Proofs are in [PROOF.md](PROOF.md); the predictions are in [PREREGISTRATION.md](PREREGISTRATION.md).

## Owner prompt (2026-11-09)

> keep going and try the tournament Chan-Pak with offset moves and the 15-matrix orbit proof

## Concurrent work: read first and credited

**Session236** (906c04d05a, `session236_chanpak_offsets_orbits_20261109`) is on the same prompt and pushed first. It
holds, and this packet does not re-claim:
- **Fixed-port padding with one final bridge:** H = C + 2^m A.
- **The positive-offset barrier** for that mechanism. Its offsets C are themselves H-values, so 7, 21 and negative
  offsets cannot occur. A = 1 forces C = 1.
- **The one-apex zero-tail algebra.**
- **The periodic-suffix drift obstruction**, valid at every length.
- **The all-word tower-ray rigidity** for the 57 ternary block matrices. This upgrades S18's finite THM-5070 (5).

Other inheritance:
- the even-graph offset proofs THM-5040 (opus-S89) and THM-5055 (mac-mini);
- THM-5070 (S18) port chains;
- THM-4875 (D), binary one-arc joins to 2^32;
- THM-4950 and THM-4952, ternary.

## What the even-graph proof needs, and what tournaments lack

Both even-graph proofs (THM-5040, THM-5055) build every odd count by writing leading bits on top of an arbitrary
smaller count. The moves are:
- the cone, d -> 2^n + d;
- the E_2-join, d -> 3*2^n - 2 + d;
- the hub moves, d -> d + (2^r - 1) 2^(n+l-r).

Each has slope 1, and its offset depends only on the order and a local degree, for every graph in a broad class.

For tournaments, every port move with an acyclic block adds x(q_o E + q_e O) (THM-5070 (1)). The port responses E, O
are path counts only when the port lies on every odd cycle, that is, for one-apex tournaments. So there is no
tournament cone. Session236's barrier records the consequence for single-bridge rays.

## The fixed-port tail family (THM-5085)

**Construction.** Keep one port u of an arbitrary core tournament T0 for ever, and append letters:
- **S:** a sink, beaten by everyone;
- **B:** a vertex w beaten by everyone except u, with w -> u (a bridge back into the port).

Letters I_t (t transitive vertices, then a bridge) are literally the words S^t B.

**(1) Transfer (PROVED).** With h = I(Omega(T), x), R = E + O, Delta = E - O and h_u = I(Omega(T0 - u), x):
- S: (h, R, Delta) -> (h, 2R, 0);
- B: (h, R, Delta) -> (h + x(R + Delta)/2, 2R - h_u, -h_u).

The value h_u never changes. Checked against brute force on 800 random tails of order <= 15: 0 mismatches.

**(2) Closed form (PROVED).** Let I be the set of bridge positions, N(I) = Σ_{i in I} 2^(i-1), and
Q(I) = Σ over pairs j < i in I of c(i - j), with c(1) = 2 and c(d) = 2^(d-1) for d >= 2. Then

    h = h0 + (x/2) [ R0 N(I) + Delta0 [1 in I] - h_u Q(I) ].

This generalises the one-apex recurrence, which is the core T0 = {u} with (h0, R0, Delta0, h_u) = (1, 1, -1, 1).
Session236's single-bridge rays are the case |I| = 1.

**(3) Signed offset transport (PROVED).** Prepending s sinks transports the offset o(I) = h0 - (x/2) h_u Q(I) along
the ray o(I) + 2^s (x/2) R0 N(I), for every s >= 1. This is the tournament counterpart of Lemma 6 in THM-5055.
- **Unlike single-bridge rays, multi-bridge rays carry negative offsets.** For example, T0 = {u} with I = {1, 2} gives
  3 * 2^(s+1) - 1.
- **Pure dyadic leading parts force a single bridge.** Then the offset is h0 > 0, an H-value (Session236's barrier,
  extended to all tails). So the values just below powers of 2 always need non-dyadic leading parts.

## Coverage (FINITE-EXACT, pre-registered)

| cores (all classes, all ports) | base | misses (complete list up to the bound) | bound | core states / search states |
|---|---|---|---|---|
| order <= 7 | 2 | **7, 21** | 10^6 | 2,831 / 57,330,333 |
| order <= 7 | 3 | **10, 31** | 10^6 | 2,901 / 34,006,533 |
| order <= 6 | 2 | 7, 21 and 18 sparse values (34869, 85885, ..., 949777) | 10^6 | 273 / 20,372,148 |
| order <= 6 | 3 | 10, 31, 5401 | 20,000 | 278 |
| order <= 5 | 2 | hundreds | 20,000 | 46 |

**Reading the table.**
- Every odd value below 10^6 other than 7 and 21 is the Hamiltonian path count of a fixed-port tail over a core of
  order <= 7.
- Every v = 1 (mod 3) below 10^6 other than 10 and 31 is a base-3 packing count of one.
- This is a second explicit family, after THM-5070's port chains, that realises both Chan-Pak spectra to 10^6 with
  only the permanent gaps missing. Binary is already known to 2^32 (THM-4875 (D)); the ternary bound 10^6 equals
  THM-5070's.
- The one-apex core alone covers 49,728 of the 99,500 odd targets in [1001, 200000]. Per dyadic window its share
  drifts from 0.56 at 2^6 to 0.497 at 2^19, and whether it has a limit is open (HYP-9655 (C)).
- A set of at most 81 core states covers [1001, 200000] (greedy; 86 in base 3). At least 17 are needed (20 in base 3):
  the program finds that many targets with pairwise disjoint covering-core sets. Pre-registered Q4 (at most 40) is
  therefore UNDECIDED.

## The 15-matrix orbit (HYP-9640 (A))

**Not proved.** Session236 has the all-word tower-ray rigidity: no word over the 57-matrix alphabet other than a tower
power ever preserves the tower direction. This rules out the digit route along the tower, at every length.

This packet adds a second explicit base-3 family with a closed form (THM-5085 (2)) and coverage to 10^6. Its values are
quadratic binary forms h0 + (3/2)(R0 N - h_u Q(N) + Delta0 [N odd]) over finitely many core parameters.

HYP-9655 (A) restates both Chan-Pak conjectures as the cofiniteness of finitely many such forms. The order-6 cores fail
sparsely and the order-7 cores succeed to 10^6. So the core order needed may grow slowly with the target (HYP-9655 (B)).

## Files

| file | purpose |
|---|---|
| `port_stats.c` | (h, E, O, h_u) in bases 2 and 3 for every class and port (input: classgen output) |
| `fp_cover.c` | coverage of the fixed-port tails over all cores (letters S and I_0, I_1, I_2, where I_t = S^t B) |
| `fp_greedy.c` | greedy core cover and the disjoint-cover-set lower bound (Q4) |
| `fp_check.py` | brute-force check of Lemma 1 on random tails, using the S18 `packing_count.c` |
| `fixedport.py` | exploratory Python version of `fp_cover.c` |
| `out/` | every output quoted above; the SHA-256 and line count of the regenerable port-statistics file |

**Reproduce.**
1. Generate the classes: `gcc -O2 -o classgen ../first_creators_genealogy_20261107/classgen.c && ./classgen 8 cls`.
2. Compute the port statistics with an order prefix:
   `gcc -O2 -o port_stats port_stats.c`, then
   `for n in 1 2 3 4 5 6 7 8; do ./port_stats cls_$n.txt $n | sed "s/^/$n /"; done > ports_le8_n.txt`.
3. Run the coverage: `gcc -O2 -o fp_cover fp_cover.c && ./fp_cover ports_le8_n.txt 2 1000000 2 1 7`.
