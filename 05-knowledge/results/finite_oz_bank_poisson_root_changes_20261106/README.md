# Finite O/Z tail banks: residue-resolved Poisson law, redundancy, root-change trichotomy (kind-pasteur-2026-11-06-S13)

Owner's question: can a finite bank of unit-start tournaments cover all sufficiently large odd counts through its O/Z
extensions?

Answer recorded here: almost certainly not (HYP-9480 (C)).
- One free O/Z block has bounded mean multiplicity.
- A single seed's counts follow a Poisson law resolved by residue class, which also explains the mod-7 bias of the
  one-apex image.
- Bank members are positively correlated through exact redundancy identities.
- The growth that coverage needs comes from root changes (HYP-9480 (D)).

Proved ingredients: THM-4888.

| file | content |
|---|---|
| `PREREGISTRATION.md` | P1-P5, pushed in d36f57c903 before any confirmatory run. It names the theorem THM-4887; the number was released to opus-2026-11-06-S83 after a reservation race, and the theorem is THM-4888 (48ba576e1f). |
| `PREREGISTRATION_P6.md` | P6 (residue-resolved law), pushed in e3fdff4ceb with `residue_law.c` and `p6_seeds.txt` before the P6 runs |
| `PROOF.md` | THM-4888 proofs: (A) mod-3 law, (B) involution rigidity, (D)2-(D)3 redundancy identities; (C) residue confinement and parity-split capacity and (D)1 JZ = OJ are KNOWN (Session 206 weighted-parameter-capacity, Session 205 audit) and were found independently. HYP-9480 laws, exploration tables and outcomes. |
| `seeds.py` | physical states from the rooted edge under O, Z, J (Session 198's complete grammar). Writes `p1_seeds.txt` (109 seeds), `bankA.txt` (85 seeds), `root.txt` and `seeds_summary.json`. |
| `density.c` | exact tail images by DFS over the parameter m (Session 205's tail formula): per-seed hit counts (mode 1) or bank misses per class mod 3 (mode 2) |
| `residue_law.c` | per seed and class mod M: odd integers, parameters with multiplicity, and distinct values in a window |
| `root_changes.c` | H-values of t . w_1 (J w_2) ... (J w_(k+1)) for k <= K: misses and dyadic-window miss counts |
| `verify.py` | P5 identity checks (they raise on failure) and the P1-P4 and P6 decisions; writes the receipts |
| `p1_density.txt`, `p2_root.txt`, `p3_bankA_w26.txt`, `p3_bankA_w27.txt`, `root_changes_1e8.txt`, `p6_root_w25.txt`, `p6_root_w26.txt`, `p6_seeds_w24.txt` | raw outputs of the confirmatory runs |
| `VALIDATION_NORMAL.json`, `VALIDATION_OPTIMIZED.json` | receipts of `verify.py` and `python3 -O verify.py` (byte-identical) |
| `audit/` | independent audit |

## Reproduce

```bash
python3 seeds.py
gcc -O2 -o density density.c
./density p1_seeds.txt 33554432 1 16777216 33554432 > p1_density.txt
./density root.txt 268435456 1 134217728 268435456 > p2_root.txt
./density bankA.txt 268435456 2 67108864 134217728 > p3_bankA_w26.txt
./density bankA.txt 268435456 2 134217728 268435456 > p3_bankA_w27.txt
gcc -O2 -o root_changes root_changes.c
./root_changes 100000000 2 > root_changes_1e8.txt
gcc -O2 -o residue_law residue_law.c
./residue_law root.txt 33554432 67108864 63 > p6_root_w25.txt
./residue_law root.txt 67108864 134217728 63 > p6_root_w26.txt
./residue_law p6_seeds.txt 16777216 33554432 63 > p6_seeds_w24.txt
python3 verify.py --out VALIDATION_NORMAL.json
python3 -O verify.py --out VALIDATION_OPTIMIZED.json
```

`p6_seeds.txt` is frozen (the P1 seeds without an odd prime >= 5 in gcd(d, R)). `seeds.py` does not regenerate it.

**Line endings.** The packet's `.gitattributes` keeps every file LF on checkout, and `verify.py` hashes LF-normalised
bytes. The sha256 that `PREREGISTRATION_P6.md` gives for `p6_seeds.txt` (a57b20c5...) is of the CRLF bytes the file had
on Windows when it was frozen. The committed LF file, with the same 102 lines, hashes to 2ef82297a0ea99f1.... The hash
for `residue_law.c` (10db375b...) is already of LF bytes.

Approximate single-core times on this machine, which was shared during the runs:
- P1: under a minute.
- P2: one minute.
- Each P3 window: about four minutes.
- P4: about 20 CPU-minutes (about an hour of wall time on the shared machine).
- P6: under a minute.
- `verify.py`: about four minutes per run (identity A over 2^26 values in numpy; identity B over 194 seeds x 2^16
  parameters in pure Python).

Memory: at most 3 x 6.25 MB bitmaps for `root_changes`, 16 MB for `density`, and a few GB of numpy arrays for
`verify.py`. No binaries are committed.

## Status

- THM-4888 (A), (B), (D)2-(D)3: PROVED (proofs in `PROOF.md`). (C) and (D)1: KNOWN (Session 206 weighted-parameter-capacity, pushed during this session; Session 205's audit), found independently. The P5 checks are sanity checks, not evidence for the proofs.
- HYP-9480 (A): CONFIRMED out of sample (P1, P2, P6a, P6b).
- HYP-9480 (B): exploratory.
- HYP-9480 (C): OPEN; P3 CONFIRMED.
- HYP-9480 (D): OPEN; P4 CONFIRMED (no two-root-change miss in (3,780,967, 10^8]).

Not claimed: an impossibility theorem for every finite bank, coverage by any family, or Chan-Pak.
