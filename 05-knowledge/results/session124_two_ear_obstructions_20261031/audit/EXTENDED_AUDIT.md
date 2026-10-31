# Extended audit: occupancy, modular lower mode, radix27/81

**Verdict:** the new results pass independently, with the stateless/fixed-assigned-digit scope preserved. Nothing here restricts the spectrum of arbitrary tournaments.

## Exact independent replay

`extended_audit.py` uses the independently audited tensor-sum matrices from `independent_audit.py`, not producer imports, saved Krylov rows, numerical LP, or cached candidate solutions. For each certificate it reconstructs the word operator, recurrence-defect row and all required Krylov iterates. The multipliers are used solely as exact rational witness data.

| Certificate family | Exact exclusions | Coordinate inequalities | Eligible words remaining |
|---|---:|---:|---|
| length3, conserved units |48|1,728|000|
| length4, conserved units |220|7,920|0000,0303,1212,2121,3030|
| length2, occupancy + units |8|288|00|
| length4, occupancy + units |224|8,064|0000|

Every certified universe is checked explicitly, including distinctness, completeness, and exclusion of slow words. Every separating constant is strictly negative; every coordinate of its separating row is nonnegative. The three legacy digit/H-bound multipliers in all48 radix27 certificates are zero. Thus these individual-template exclusions require neither a bound on the constant digit nor H>1, but they apply only to the stated eligible words.

The occupancy relations have rank8. Their rational nullspace construction uses unit pivots and has integral free-variable columns, so it is a complete **integer** kernel basis of rank28. The exact reachable linear span also has dimension28. All four letter operators preserve this kernel, and the initial vector lies in it.

For all four complementary words03,12,21,30, the independently reconstructed row q is integral, has value−4 on the initial state, and satisfies448 total modulo3 transition tests on the full integer kernel. The scalar annihilator rows vanish identically. The derivative at1 of `z(z-1)(z-4)^2(z-9)/8` equals−9. Therefore the linear1-mode coefficient is exactly `D=-q(v)/9`; its universal nonvanishing follows from `q(v)=2 mod3`, not from an unjustified numerical fit.

Normal and optimized independent outputs are `../receipts/extended_audit.normal.out` and `../receipts/extended_audit.optimized.out`; they include source-certificate hashes. The original certificate verifiers, including their mutations, are replayed by `../verify_all.py`, with individual outputs in `../receipts/`.

## Proof and quantifier checks

1. **Occupancy swap:** ignoring last-vertex good/bad flags, exchange Q and R while reversing ear order. The first-vertex requirements match exactly; P stays fixed. This proves the eight cross-order identities. The relation is not a claim that accepted xy/yx path counts agree: acceptance still distinguishes the last flags and the direct ear edge.
2. **Separating identities:** any constant-digit recurrence implies `c(M-I)(M-BI)M^j v=0` for every j≥0. Occupancy relations vanish and conserved units equal1. A nonnegative separating row then has negative value on an alleged nonnegative actual state, impossible. No numerical solver's feasibility status is used.
3. **Forced mode:** the scalar annihilator permits9^k,4^k,k4^k,1,k (after the zero transient). Its polynomial filter annihilates all but the k term. This is an exact Jordan-mode argument, including k≥1 and the vanishing zero transient; it does not assume diagonalizability.
4. **Radix81:** even the weaker220-certificate family leaves at most five eligible templates, fewer than81. The stronger224-family leaves only0000. Either suffices for a complete stateless length4 bank obstruction.
5. **Combined2≤B<243 theorem:** a complete stateless bank assigns one fixed physical word to each centered digit and the same scalar action holds at every descendant of a seed H>1. An empty word is impossible since it would require `|d|=(B-1)H>B-1`. Every nonempty repeatable word has strictly positive B^k coefficient, so the spectral filtration gives `B=2^L` or `B=3^L`. Unique prime factorization forces the same length for every bank word, even if unequal lengths were initially permitted. Binary powers are excluded by the previous theorem; the remaining ternary powers below243 are exactly3,9,27,81. The fixed-bank obstructions cover them. This is not a claim about arbitrary changing physical banks, adaptive templates, or phase-dependent digit assignments for the individual-template certificates.
6. The all-zero eligible exceptions are **not** certified centered-digit constructions. Conversely, slow words omitted from the eligible universe cannot be declared impossible for unrestricted constant d: universal-sink words can maintain a constant count with a sufficiently negative d. Their exclusion in the complete-bank theorem specifically uses H>1 and the centered digit bound.

## Publication hygiene coordination

The producer has been asked to remove or explicitly label the obsolete fixed-z calculations in `check.py` main, and to mark the old “length≥3/4 remains open” statements as superseded by the subsequent finite-length results. These are presentation issues; the corrected affine conditions, exact new certificates, and combined theorem do not depend on the obsolete calculations.


## Curated-copy hygiene resolution

The publication copy removes the obsolete fixed-z computation from `code/check.py`
and marks earlier open-length statements as superseded by the audited length3/4
results. The original research/audit history is preserved outside this package.
Run `../verify_all.py` to regenerate all included proof-checker and audit outputs;
no workspace paths or numerical LP dependencies are needed.
