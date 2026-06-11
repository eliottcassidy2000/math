# THM-471: the B3 general-shape enumerator — m=2 Schipperus-cutoff probes and the first m=3 (open-case) data

**Status:** STUB — CLAIMED by kind-pasteur-2026-06-11-S1. Implements THM-460 B3's
general-shape recursive enumerator (POKE Task 1.2). Computations this session; do
not build on this file until status changes.

## Claimed content

Tuple-form of the B3 grammar (exponents δ ∈ N^m, position significance δ[m−1] high):
Bin(0)=point; finite peel (least nonzero index 0): two Bin(δ−e₀) blocks, common
cross-split position, internals strictly lex-below; limit peel (least nonzero index
i ≥ 1): M order-separated parts Bin((δ−e_i) + (j,…,j) in positions 0..i−1), j=1..M.
BT(m,M) = parts Bin((j,…,j)), j=1..M. Size analysis: at m=3, s=2 the j=1..M march
makes BT(3,2) have 1873 leaves > 256 ambient points (vacuous — the m=3 analogue of
THM-460 D's miniature-design guard); feasible probes are M=1 (16-leaf shape, a
legitimate rung: B3(i) holds per M including M=1) and the truncated j-from-0 M=2
variant (43 leaves). Plan: validate at m=1 against THM-460 D (R(1,2)=3 shadow at
M=1; M=2 SAT region), then m=2 sweeps (Schipperus positivity FORCES cutoffs via
THM-460 C2 — first such number ever computed if reached) [HYP-2394], then the m=3
probe at (s,c)=(2,2),(2,3) [HYP-2395].

Script: `04-computation/erdos592_shape_miniatures_kp0611.py` (+ .out).
