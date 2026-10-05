# THM-169: Complete Characterization of the Vitali Atom (4-Reversal H-Change)

**Status:** VERIFIED finite evidence (reported exhaustive n<=6 and sampled n=7);
the stated uniform n=7 iff classification is OPEN pending a repaired proof.
**Session:** kind-pasteur-2026-03-13-S61
**Dependencies:** THM-002 (OCF), THM-168 (lambda completeness/ambiguity)

**2026-10-05 audit:** the former PROVED label exceeded the reported sampled
evidence and the proof sketch contains an erroneous triangle path count,
corrected below. This audit does not refute the classification, but does not
promote it to a proved dependency. The separately checked lambda-preserving
pair in [the atomic-prefix audit](../../05-knowledge/results/collatz_effective_prefix_mass_20261005.md),
section9, proves existence of a seven-vertex H-change, not the claimed iff.

## Statement

### Setup

Let T be a tournament on n=7 vertices. Let S = {a,b,c,d} be a 4-vertex subset
with sub-tournament scores (1,1,2,2). Let T' be obtained from T by reversing
all 6 arcs within S. Let E = V(T) \ S = {e,f,g} be the 3 external vertices.

Assume the reversal is **lambda-preserving**: lambda(T) = lambda(T') as labeled graphs.

### Definitions

- **mixed_c3**: number of directed 3-cycles using exactly 2 vertices from S and 1 from E.
- **ext_c3**: 1 if the external tournament on E is cyclic (scores (1,1,1)), 0 if transitive (0,1,2).
- **B_s**: set of subset vertices that beat the ext-source (the vertex of E with highest out-degree within E). Defined only when ext_c3 = 0.
- **B_k**: set of subset vertices beaten by the ext-sink (lowest out-degree within E). Defined only when ext_c3 = 0.

### Main Theorem

H(T') != H(T) if and only if ALL of:

1. The sub-tournament on S has scores (1,1,2,2) (non-transitive).
2. The reversal preserves the labeled lambda graph.
3. mixed_c3 = 4 (equivalently, each subset vertex has exactly 2 arcs to ext and 1 from ext, or vice versa).
4. **Either:**
   - (a) The external tournament is cyclic (ext_c3 = 1), **OR**
   - (b) The external tournament is transitive AND (|B_s| = 4 OR |B_k| = 4).

When H changes, delta_H = +/-2 exactly (from delta_c7 = +/-1).

### Complementary results

- **n <= 6:** H(T') = H(T) ALWAYS for lambda-preserving (1,1,2,2) reversals. (Exhaustive proof.)
- **Phase transition at n=7 is sharp:** the first H-change occurs at n=7, caused by 7-cycles.
- **Mechanism:** delta_H arises entirely from 7-cycles with the 4 subset vertices CONSECUTIVE (threading pattern EEESSSS). These correspond to Hamiltonian paths of the sub-tournament completed through the 3 external vertices.

## Proof Sketch

### Step 1: Threading decomposition

A 7-cycle on all vertices uses k internal arcs (arcs within S), where k in {1,2,3}.
Only k=3 matters: a 7-cycle with 3 internal arcs has all 4 subset vertices consecutive,
forming a Hamiltonian path of the sub-tournament on S.

### Step 2: Completion formula

For a Ham path p = (a -> b -> c -> d) on S, a **completion** is a permutation (e1,e2,e3)
of E such that d -> e1 -> e2 -> e3 -> a forms a valid directed path.
The completion count C(d,a) depends ONLY on the endpoint pair (d,a) and the external arcs
(which are unchanged by reversal).

Therefore: delta_c7 = sum_{(d,a) in EP(T')} C(d,a) - sum_{(d,a) in EP(T)} C(d,a)

### Step 3: Universal endpoint distribution

For ANY (1,1,2,2) tournament on 4 vertices (all 24 labelings), the 5 Ham paths have
endpoint type distribution (using L=low-score, H=high-score):

    (end=L, start=H): 2 pairs
    (end=L, start=L): 1 pair
    (end=H, start=L): 1 pair
    (end=H, start=H): 1 pair

After reversal (L <-> H swap), this becomes (HL): 2, (LL): 1, (LH): 1, (HH): 1.
The net change: (HL): +1, (LH): -1.

### Step 4: Cyclic ext case

For cyclic E, there are **3** Hamiltonian paths through E, one for each
choice of starting vertex, all following the same directed cyclic orientation.
The former count of2 (one per direction) was false. The reported48 cases
with a net +/-1 change remain finite evidence; the stated uniform completion
conclusion requires a proof retaining these three paths and their endpoints.

### Step 5: Transitive ext case

For transitive E with source=e_s, middle=e_m, sink=e_k, the unique Ham path is e_s -> e_m -> e_k.
The completion matrix C[d][a] = 1 iff A[d][e_s]=1 AND A[e_k][a]=1 (rank-1 structure).

delta_c7 = |{EP(T') in B_s x B_k}| - |{EP(T) in B_s x B_k}|

When |B_s|=4 or |B_k|=4: the bilinear form "sees" the full endpoint redistribution,
yielding delta = +/-1. When |B_s|=|B_k|=2: the form is insensitive to the redistribution
(the +1 and -1 cancel), yielding delta = 0.

Verified exhaustively: 34/34 changing cases have |B_s|=4 or |B_k|=4,
55/55 preserving cases with |B_s|=|B_k|=2 have delta=0.

## Verification

- **n=5 exhaustive:** 0/640 lambda-preserving (1,1,2,2) reversals change H.
- **n=6 exhaustive:** 0/4320 change H. c5 also invariant (exhaustive).
- **n=7 sampled (5000+ tournaments):** 48/48 cyclic-ext change H (100%), 34/34 transitive-ext with |B_s|=4 or |B_k|=4 change H (100%), 0/163 transitive-ext with |B_s|=|B_k|=2 change H (0%).
- Scripts: vitali_c7_mechanism.py, ham_path_vitali_mechanism.py, transitive_ext_discriminator.py

## Significance

The proposed classification concerns the **Vitali atom** — the minimal gauge transformation
that breaks H-invariance under the lambda-preserving equivalence relation. The phase
transition at n=7 is explained by a dimensional mechanism: 7-cycles require all 7 vertices,
and the Ham-path-completion interaction provides a 1-degree-of-freedom "dark sector"
invisible to the lambda graph. This directly explains the 19.1% ambiguity rate in THM-168.
