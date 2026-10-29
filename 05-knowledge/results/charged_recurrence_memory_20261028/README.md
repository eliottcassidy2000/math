# Charged recurrences, sharp port thresholds and growing linear memory

**2026-10-28. PROVED scoped statements in the notes / FINITE-EXACT receipts /
CITED classical and repo mechanisms / OPEN general attachment closure.
Same-agent proof review and alternate computations; not independently audited.**

The main result is an exact hierarchy of observations: imposing the multiplier
constraint at tree depth j first changes an ordered port count at arity p^j.
The complete charged cycle recurrence can therefore be reduced, for fixed
point counts, to character blocks indexed by their deepest nontrivial depth.
For bounded arity, deeper blocks contribute exactly zero. This compression
does not classify joint transporters or arbitrary complement/re-marking
contexts.

This session takes the existing loop kernel as an input and changes the
counting direction. It also recovers two sparsely cited mechanisms from
earlier work, then tests an odd-field Frobenius extension.

| result | statement and boundary | proof / checks |
|---|---|---|
| Charged cycle recurrence | full multiplier sectors; Hol counts from the identity sector | [FOURIER F1](FOURIER.md), [charged.py](charged.py) |
| Sharp visibility | successive constraints first separate at 1,p,p^2,...; exact leading gaps | [FOURIER F2--F3](FOURIER.md), [jets.py](jets.py) |
| Low-order Hol counts | p+2k-2 and 8k^2+(10p-30)k+p^2-10p+22 | [FOURIER F4](FOURIER.md) |
| Minimal arity order | p^k+1-(p-1)(p-2)/2 for affine wreath actions, k>=2 | [MEMORY M1--M2](MEMORY.md), [memory.py](memory.py) |
| Ternary recurrence | explicit characteristic roots 3,9,...,3^k for centralizer ports | [MEMORY M3](MEMORY.md) |
| Actual complement memory | L_j=(-M)^j(I-M^2); exact complement-only set and linear minima | [COMPLEMENT N1--N3](COMPLEMENT.md), [renewal.py](renewal.py) |
| Frobenius boundary | F_27 semilinear offspring law; arity-rank deficit 109 versus affine 325 | [FROBENIUS S1--S3](FROBENIUS.md), [frobenius.py](frobenius.py) |
| Cycle-free rook recurrence | actual native partner H counts with injectivity and closed-cycle exclusion | [ROOKS R1--R3](ROOKS.md), [rooks.py](rooks.py) |

Two concrete controls: at p=7, depth two, C has order 83,349 and Hol has
order 9,261; their two-port counts are 49 and 113. The root-child stabilizer
and Hol agree through six ports, then differ by 1,290,240 at seven ports.
The 81/27 tournaments still have identical labelled cyclic triples, and
their complete automorphism counts are retested.

[RECOVERY](RECOVERY.md) records the inheritance pass, concept board, three
recovery waves, the sparse-citation selection method, typed connections,
and three concrete open problems. [AUDIT](AUDIT.md) states exact universes,
source laws, independent computational routes, and limitations.

Reproduce from the repository root:

```text
python -X utf8 -B 05-knowledge/results/charged_recurrence_memory_20261028/replay.py
python -X utf8 -B -O 05-knowledge/results/charged_recurrence_memory_20261028/replay.py
python -X utf8 -B 05-knowledge/results/charged_recurrence_memory_20261028/seal.py
```

No A000568 prefix beyond the existing work is claimed. The new recurrences
count specified marked-atom actions and actual finite complement observations.
The smallest general carrier closed under re-marking, complement and strong
substitution remains open.
