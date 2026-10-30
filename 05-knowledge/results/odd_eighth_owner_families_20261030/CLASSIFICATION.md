# Four odd-eighth centers: fourteen complete infinite functions

**PROVED uniform tail / FINITE-EXACT complete finite head.**

Let V(E,m)={1,...,7} union E union {8m}, with E any subset of
{9,...,15} and m>=1. Every row gives the complete function M(V(E,m)).
For EVERY m>=5 its value is m/(8m+s); the first four values follow.
The all-size proof is in [PROOFS, O1-O3](PROOFS.md).

Pattern positions are speeds 9,10,11,12,13,14,15 in that order:
`1` means present, `0` absent, and `*` arbitrary. Each row is exactly
one Boolean cube; the disjoint rows cover all 128 choices.

| Pattern | Sources | m=1 | m=2 | m=3 | m=4 | s for m>=5 |
|---|---:|---:|---:|---:|---:|---:|
| `0******` | 64 | 1/9 | 2/17 | 3/25 | 4/33 | 1 |
| `100****` | 16 | 1/10 | 2/19 | 1/9 | 4/35 | 3 |
| `110****` | 16 | 1/11 | 2/19 | 1/9 | 4/35 | 3 |
| `101*0**` | 8 | 1/10 | 1/10 | 3/29 | 4/37 | 5 |
| `11100**` | 4 | 1/12 | 2/21 | 3/29 | 4/37 | 5 |
| `11110**` | 4 | 1/13 | 2/21 | 3/29 | 4/37 | 5 |
| `101*1*0` | 4 | 1/10 | 1/10 | 1/10 | 4/39 | 7 |
| `11101*0` | 2 | 1/12 | 2/23 | 3/31 | 4/39 | 7 |
| `1111100` | 1 | 1/14 | 2/23 | 3/31 | 4/39 | 7 |
| `1111110` | 1 | 1/15 | 2/23 | 3/31 | 4/39 | 7 |
| `101*1*1` | 4 | 1/10 | 1/10 | 1/10 | 1/10 | 9 |
| `11101*1` | 2 | 1/12 | 1/12 | 1/11 | 4/41 | 9 |
| `1111101` | 1 | 1/14 | 2/25 | 1/11 | 4/41 | 9 |
| `1111111` | 1 | 1/16 | 2/25 | 1/11 | 4/41 | 9 |

No equality is inferred merely from a finite prefix: matching s proves
the entire remaining tail agrees. Distinct rows already differ for m<=5.

| Consumer and allowed continuation | Exact number of states |
|---|---:|
| Present answer at any fixed m>=5 | 5 |
| Entire present function, all m>=1 | 14 |
| Future band additions, queries m>=5 | 16 |
| Future band additions, queries all m>=1 | 128 |

The last two are minimal deterministic carriers, with constructive
separating continuations for every retained bit ([O4](PROOFS.md)).
These state counts concern this declared family and operation menu.
In particular the number 14 here does not settle LRC(14).

Reproduce with `python -X utf8 -B 05-knowledge/results/odd_eighth_owner_families_20261030/verify.py`,
then run `render_table.py` in the same directory. The full exact times,
finite overlap and native hostile controls are in [certificate.json](certificate.json).
