# Session133: eight fixed compiled counters and two small direct counters

**Actual finite ordered FRACTRAN lists, with all-input proofs and exact replays.**
The unchanged A000568 source shrinks from **5,543 to 4,294 fractions**.

| Sequence | Fractions | Input | Output |
|---|---:|---|---|
| tournament classes A_n | 4,294 | `149*5^n` | terminal `v_3=A_n` |
| self-converse S_n | 4,744 | `163*5^n` | terminal `v_3=S_n` |
| converse-merged M_n | 5,152 | `179*5^n` | terminal `v_3=M_n` |
| strongly connected C_n | 5,663 | `199*5^n` | terminal `v_3=C_n` |
| ordered tournament pairs | 4,296 | `149*5^n` | terminal `v_3=F_n(4,1)` |
| 2-colored pairs | 4,437 | `157*5^n` | terminal `v_3=F_n(4,2)` |
| 2-colored tournaments | 4,436 | `157*5^n` | terminal `v_3=F_n(2,2)` |
| minimum-core extremizers e(N), N>=6 | 829 | `53*5^(N-6)` | terminal `v_3=e(N)` |
| minimum core H, g(N), N>=3 | **20** | `43*5^(N-3)` | exactly `41*2^g(N)` |
| same g(N), tail N>=6 | **14** | `11*5^(N-6)` | exactly `41*2^g(N)` |

All ordinary n inputs are nonnegative. C_0=0 and C_1=1. The compiled programs
retain work registers; pure terminal forms are claimed only for core20/core14.
No program consumes a sequence lookup table. Every file is below 100KB.

The core-floor and extremizer meanings use
[the separately audited core theorem](../core_reduction/THEORY.md).
Their scalar outputs sharpen the size bounds relevant to H=7 or21; they do
**not** replace the local H-spectrum certificates or prove open completeness.
See [THEORY.md](THEORY.md) for contracts, all-input termination, exact runtime
of the small programs, and the information-loss boundaries.

```bash
python3 verify.py
cp VERIFIED.json VERIFIED.normal.json
python3 -O verify.py
cmp VERIFIED.normal.json VERIFIED.json
sha256sum -c SHA256SUMS

# Read delivered fractions and literally scan actual integers; no compilation:
python3 run_program.py self_converse 3   # 2
python3 run_program.py converse_merged 3 # 2
python3 run_program.py strong 3           # 1
python3 run_program.py extremizer_tail 5 # e(11)=3
```

`VERIFIED.json` separates native arithmetic, primitive counter, prime-exponent
fraction and literal integer coverage. Higher source-native terms are not
misrepresented as completed FRACTRAN runs. Checks remain active under `-O`;
independent signed-permutation, whole-orbit and pair-orbit censuses plus hostile
controls corroborate the proofs. The two core lists are in `core_full20` and
`core_tail14` fraction files, with independent direct runners in `core_counter.py`.

## Exact execution coverage

All eight source-native programs are checked at inputs 0..12. The following
ranges are **actual counter/fraction replays**, not merely formula evaluations.

| Program | Primitive IR inputs | Exponent-fraction inputs | Literal-integer inputs | Largest recorded IR steps |
|---|---|---|---|---:|
| `a000568` | 0..6 | 0..5 | 0..3 | 22,926,170 |
| `colored_pair_tournaments` | 0..3 | 0..3 | 0..2 | 446,933 |
| `colored_tournaments` | 0..4 | 0..4 | 0..2 | 579,378 |
| `converse_merged` | 0..5 | 0..4 | 0..2 | 491,598 |
| `extremizer_tail` | 0..12 | 0..12 | 0..2 | 4,563 |
| `pair_tournaments` | 0..4 | 0..3 | 0..2 | 6,276,086 |
| `self_converse` | 0..8 | 0..7 | 0..3 | 6,338,264 |
| `strong` | 0..5 | 0..4 | 0..3 | 993,943 |

Core20 is replayed at N=3..24 and core14 at N=6..24, in both actual integers
and exponent vectors. Every run checks the exact terminal shape and step formula.
The verifier uses a 100,000,000-step guard for each IR/exponent run and a
1,000,000-step guard for each on-disk literal run; these are test safety guards,
not asserted all-input runtime bounds. The proofs establish all-input termination.
