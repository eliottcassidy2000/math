# egcp_20261109: sharp even-graph Chan-Pak (THM-5040)

These scripts accompany [THM-5040](../../../01-canon/theorems/THM-5040-sharp-even-graph-chan-pak.md). They are
pure Python 3 with no dependencies. Run them from this directory.

| Script | What it does | Output |
| --- | --- | --- |
| `base_extract.py` | Extracts the census base from `even511.txt` (see below). Writes `base_witnesses.txt` (248 minimal-order realisers W and 20 leaf graphs L at orders 8 and 9) and `base_sets.json` (E, E0, C and L for n <= 9). | committed outputs |
| `verify_base.py` | Re-verifies every base graph by brute force (d, evenness, order, leaf, no universal vertex), checks E(n) = O(n) for n <= 9, and checks that L covers T(8), T(9) and the offsets. | `verify_base.out` |
| `lemma_check.py` | Random tests of the complement formula and Lemmas 5, 7 and 8, including the evenness criterion of Lemma 8 in both directions. | `lemma_check.out` |
| `proof_check.py N` | Walks every x through the case analysis of THM-5040 Section 4 for 10 <= n <= N and asserts each side condition. | `proof_check_24.out` |
| `construct.py N S` | Builds explicit witnesses with the proof's own recursion and brute-forces d, evenness and order. It covers all delicate values plus S random values at each order 10..N, and all top-band graphs. | `construct_12.out` |
| `cover.py N` | Independent closure under padding, cones, E_2/K_1/E_3 multi-joins, two-hub graphs, the top band and the leaf recursion. Compares E(n) with O(n), and E0(n) with the conjectured exact no-universal-vertex set. | `cover_25.out` |
| `egcp_lib.py` | Shared helpers: exact automorphism enumeration and the operations cone, E2join, SR, jump, cocktail and frontier. | |

**Regenerating the census base** (optional; the extracted files are committed):

```
gcc -O3 -o dom_even ../paleyeven_20261109/dom_even.c
./dom_even 511 16 > dom511.out
awk '$4=="E"' dom511.out > even511.txt
python base_extract.py even511.txt
```

The census reproduces THM-4996 (3): connected graphs with d <= 511 number 1, 1, 2, 6, 21, 112, 853, 11,117, 261,080,
199,053 and 1,855 on 1 to 11 vertices.
