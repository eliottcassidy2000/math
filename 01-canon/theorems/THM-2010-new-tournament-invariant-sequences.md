---
id: THM-2010
title: "A CATALOG OF SEQUENCE-VALUED TOURNAMENT INVARIANTS, with OEIS identifications — ten candidate-NEW integer sequences and a clean skew/adjacency asymmetry. Computed exhaustively over iso classes n=3..6 (exact integer arithmetic), the corpus' invariants generate integer sequences of three kinds: (a) INVARIANT-VALUE-COUNTS |X|(n) = #distinct values invariant X takes; (b) EXTREMAL maxX(n); (c) STRUCTURAL counts. KNOWN (identified): iso=A000568; |score|=A000571; #strongly-connected=A051337; #regular=A007079; max H (Hamiltonian-path count) = 3,5,15,45 = A038375 (Moon, 'max spanning paths', hard); #self-converse=2,2,8,12=A002785; and |specS| (distinct SKEW/Seidel char polynomials) = 1,2,2,6,11,50 = the Breen-Stover-Yates Seidel-tournament-matrix count (arXiv:2406.09697). CANDIDATE-NEW (no OEIS match on 4 terms; the four marked [0] returned ZERO OEIS hits): the ADJACENCY char-poly count |specA| = 2,3,9,28 (the A-side twin of the KNOWN skew count |specS| — a clean gap); the Redei-spectrum size |H| = 2,3,7,19 (#distinct Hamiltonian-path-counts); |cyc| = 2,3,9,32; |R| = 2,2,6,8; |disc| = |var(lambda^2)| = 1,2,2,5 (equal, both skew-spectral); |arb_inv| = 2,4,12,55 [0]; max-total-arborescences = 3,10,55,333 [0]; max|R| = 3,3,15,15,147 [0]; metagraph edges |E(G_n)| = 1,5,30,290 [0]; metagraph 1-WL color classes = 1,2,10,34. THE ASYMMETRY: the SKEW characteristic-polynomial count is a studied sequence but the ADJACENCY one is not — |specA| is strictly finer (28 vs 6 at n=6, since the adjacency matrix is not antisymmetric and keeps the score/hierarchy data the skew form discards). |arb_inv| = 2,4,12,55 is a NEAR-COMPLETE invariant (55 of 56 iso classes at n=6). |disc| and |var(lambda^2)| taking the SAME count (1,2,2,5) reflects both being functions of the skew spectrum."
status: >
  All first terms VERIFIED-EXACT over every iso class n=3..6 (bit-packed canon, exact integer
  Faddeev-LeVerrier / Bareiss / Held-Karp; script + .out).  OEIS identifications done by direct
  API lookup (curl): the KNOWN seven match named sequences; the ten CANDIDATE-NEW returned no
  tournament-relevant OEIS match, and the four marked [0] (arb_inv, max-total-arb, max|R|,
  metagraph-edges) returned ZERO hits.  CAVEAT: four terms is short — a candidate could match a
  longer OEIS entry my query did not surface, or sit in OEIS under a different offset; the [0] four
  are the confidently-new ones.  The specA/specS asymmetry and the arb_inv near-completeness are
  exact n<=6 facts.  This is a REFERENCE CATALOG (verified data + external identifications), not a
  proof; its value is surfacing ~10 uncataloged tournament-invariant sequences and the skew/adjacency
  gap, and giving each a home for OEIS submission and reciprocal-signature study (THM-1990).
source: kind-pasteur-2026-07-21-S128c145 (owner: look for more sequence-like invariants)
depends_on:
  - THM-1965    # the invariant lattice (these are the value-counts of its invariants)
  - THM-1990    # the reciprocal-sum lens (each new sequence gets a signature)
related: [THM-1980, THM-1936, THM-1930]
external:
  - "OEIS: A000568 (tournaments), A000571 (score seqs), A051337 (strong tournaments), A007079 (regular tournaments), A038375 (max Hamiltonian paths, Moon 1968), A002785 (self-complementary oriented graphs)."
  - "Breen, Stover, Yates, 'Determinants of Seidel Tournament Matrices' (arXiv:2406.09697): distinct Seidel char-poly count 1,1,1,2,2,6,11,50 = |specS| here."
script: 04-computation/new_sequence_invariants_kps_S128c145.py (+ .out)
---

# THM-2010 — a catalog of sequence-valued tournament invariants

Every invariant `X` of tournaments generates an integer sequence in `n`. Beyond the census sequences
(A000568 etc.), three under-explored families — **value-counts** `|X|(n)`, **extremal** `maxX(n)`,
and **structural** counts — yield a batch of sequences, most of them uncataloged. Exhaustive over iso
classes `n = 3..6` (exact arithmetic).

## Known (identified)

| sequence | `n=3..6` | OEIS / source |
|---|---|---|
| iso classes `V(G_n)` | 2, 4, 12, 56 | **A000568** (tournaments) |
| `|score|` distinct | 2, 4, 9, 22 | **A000571** (score sequences) |
| `|specS|` distinct (skew char polys) | 1, 2, 2, 6 | **Breen–Stover–Yates** 1,2,2,6,11,50 (arXiv:2406.09697) |
| max H (Hamiltonian paths) | 3, 5, 15, 45 | **A038375** (Moon, "max spanning paths"; hard) |
| # strongly connected | 1, 1, 6, 35 | **A051337** |
| # regular | 1, 0, 1, 0 | **A007079** (n odd) |
| # self-converse | 2, 2, 8, 12 | **A002785** (self-complementary oriented graphs) |

## Candidate-new (no OEIS tournament match at 4 terms; `[0]` = zero hits)

| sequence | `n=3..6` | note |
|---|---|---|
| `|specA|` distinct (adjacency char polys) | **2, 3, 9, 28** | the A-side twin of the known skew `|specS|` — a gap |
| `|H|` distinct (Rédei spectrum size) | **2, 3, 7, 19** | # distinct Hamiltonian-path counts |
| `|cyc|` distinct | 2, 3, 9, 32 | # distinct cycle-count vectors |
| `|R|` distinct | 2, 2, 6, 8 | # distinct signed-Rédei magnitudes |
| `|disc| = |var(λ²)|` | 1, 2, 2, 5 | equal — both skew-spectral functions |
| `|arb_inv|` distinct | **2, 4, 12, 55** `[0]` | near-complete (55 of 56 classes at n=6) |
| max total arborescences | **3, 10, 55, 333** `[0]` | max spanning out-trees over all roots |
| max `|R|` | **3, 3, 15, 15, 147** `[0]` | (mac-mini THM-1936) |
| metagraph edges `|E(G_n)|` (wiggly) | **1, 5, 30, 290** `[0]` | edges of the flip metagraph |
| metagraph 1-WL color classes | 1, 2, 10, 34 | THM-1965's canonical vertex coloring |

## Two structural observations

- **The skew/adjacency char-poly asymmetry.** The number of distinct **skew** (Seidel) characteristic
  polynomials `|specS| = 1,2,2,6,11,50` is a studied sequence; the number of distinct **adjacency**
  characteristic polynomials `|specA| = 2,3,9,28` is not, and is strictly larger (28 vs 6 at n=6).
  The adjacency matrix is not antisymmetric, so it retains the score/hierarchy data the skew form
  discards — `specA` is the finer spectral invariant. `|specA|` is the obvious missing A-side companion
  to a published sequence.
- **`|disc| = |var(λ²)| = 1,2,2,5`.** Both are functions of the skew spectrum, so they take equally
  many values — a small witness that the "skew-spectral" invariants form one coarse layer (below the
  full `specS` count 1,2,2,6).

## Named next

- **Extend the `[0]` sequences to n=7** (max total-arb, `|arb_inv|`, metagraph edges, max`|R|` — the
  latter is `…,147` at n=7 already) and **submit to OEIS**; four terms is thin, seven would settle
  novelty.
- **Submit `|specA| = 2,3,9,28,…`** as the adjacency companion of the Seidel count, with a cross-ref.
- **Reciprocal signatures (THM-1990).** Each new sequence has a `σ = Σ 1/a_n`; the fast-growing ones
  (max H = A038375 grows like `n!/2ⁿ`) converge — fold into the reciprocal-signature table.
- **A closed form for metagraph edges `1,5,30,290`** — the flip-metagraph edge count should have a
  generating expression (per-class wiggly degrees summed, THM-1965/CLAUDE.md waggly layers).
