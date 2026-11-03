# Verification, priority, and limits

Status: **PASS, same-agent review with independent exact algorithms**.
Incoming Session144/145 have their own independent audits; those are cited,
not represented as sub-agent work requested by this session.

## Exact universes and receipts

| Check | Universe | Independent route / result |
|---|---|---|
| Uniform path-cover transfer | all1099 labelled tournaments n=1..5 and5405 births | subset DP plus block decomposition versus permutations/cuts; agrees |
| Uniform birth-word moments | every word through n=8 | direct counts versus all-order first-moment recurrence; agrees |
| Aggregate polynomial roots | n=1..12 | exact Sturm counts; general proof uses D+positive preservation |
| Parity carrier closure | all12412 forests on7 labels and all64 cuts to6 |170752 basis pullbacks; coefficients and missing terms checked |
| Carrier dimensions | n=0..20 | finite sum versus separate EGF expansion; agrees |
| All-context mean hostile | all32768 six-label sources under89 fixed contexts |2916352 literal Boolean counts; equal means, cut17 separates second moments |
| Short actual query bank |2979 discovery words,571 retained | signed block contraction gives rank571 modulo101,103 |
| Independent minimum certificate |571 actual words x768 actual sources |438528 Boolean subset-DP counts; full evaluation matrix rank571 modulo101,103 |
| Keller map | explicit rational polynomial K and quotient | determinant2, collision, composed twisted equation constant4; exact SymPy |
| Orbit degree recurrence | proof for every k; specialization controls k<=5 | degrees1,7,37,187,937 modulo1009 match proved upper bound |
| Prime table | seven specified rows | sieve/list versus totient, factorization and eligible-prime counts |
| Normalized prime records | every even4<=m<=20000 | complete sieve scan; each of15 records replayed with independent primality/primepi |
| Retained source controls | actual81/27 pair and three arithmetic orbits | automorphism enumeration, Hamilton counts, exact valuations; agrees |

The six-vertex lower bound is particularly strong: its independent verifier
does not import or compute a forest coefficient. The integer evaluation
matrix has rank571 modulo an odd prime, so rank over Q is at least571.
The proved closed-carrier dimension supplies the matching upper bound.
The768 selected sources suffice for this lower-bound certificate; they
are not described as an exhaustive census of32768 orientations.

Finite tests do not prove the all-order assertions by extrapolation. Their
proofs are the deletion/reconstruction bijection, conditional variance,
signed block replacement, closure under grafting, composition chain rule,
leading-degree induction, cyclic-gap sum, and elementary CRT construction.

## Reproduction

From the worktree root, let P=05-knowledge/results/positive_source_jacobian_prime_clocks_20261103.
Run the following commands with P expanded (PowerShell uses `$p` for it):

```powershell
$p = '05-knowledge/results/positive_source_jacobian_prime_clocks_20261103'
python -B "$p/verify.py" > "$p/verification.json"
python -O -B "$p/verify.py" > "$p/verification_optimized.json"
python -B "$p/carrier.py" > "$p/carrier.json"
python -O -B "$p/carrier.py" > "$p/carrier_optimized.json"
python -B "$p/sidecars.py" > "$p/sidecars.json"
python -B "$p/rank_input.py"
g++ -O3 -std=c++17 "$p/carrier_literal.cpp" -o "$env:TEMP/positive_jc_carrier_literal.exe"
& "$env:TEMP/positive_jc_carrier_literal.exe" > "$p/carrier_literal.json"
g++ -O3 -std=c++17 "$p/carrier_rank.cpp" -o "$env:TEMP/positive_jc_carrier_rank.exe"
Get-Content "$p/rank_input.txt" | & "$env:TEMP/positive_jc_carrier_rank.exe" > "$p/carrier_rank.json"
g++ -O3 -std=c++17 "$p/rank_literal.cpp" -o "$env:TEMP/positive_jc_rank_literal.exe"
Get-Content "$p/rank_contexts.txt" | & "$env:TEMP/positive_jc_rank_literal.exe" > "$p/rank_literal.json"
python -B "$p/seal.py"
```

The saved rank_contexts.txt is the ordered571-word list in carrier_rank.json,
preceded by571; each subsequent line begins with the word length. The seal
checks this equality. Runtime executables are outside the repository.
Python normal/-O main receipts and carrier receipts are byte-identical;
checks use explicit failures, not asserts. g++ and SymPy provide the local
exact arithmetic implementations; no interval or floating-point rank is used.

## Proof audit and corrections

- Internal parity means the number of internal vertices mod2; reversal
  sign is(-1)^(internal count+1). The total internal pool remains named.
- Closure is for fixed source-oblivious attachments, ordinary converse and
  fixed relabelling. Native path complement and adaptive marking are outside it.
- Modular rank is used as a LOWER bound. Equality at571 follows from a
  proved upper bound. The lower ranks at shorter word lengths have no
  separately proved rational upper bounds here.
- The positive-law certificate needs an actual source law. Uniform births
  fail it at every order>=5; this does not rule out other laws.
- K is the same-weight output reorder, not the unpermuted published map.
  Quotient degree, geometric degree and determinant are distinct clocks.
- A raw nonlinear orbit degree is gauge-dependent. The incoming shear
  settles every degree>=2 without creating new equivalence classes; earlier
  draft next-steps asking whether degrees2..6 occur were removed before push.
- The general fold count cannot assume no sheets escape at a target
  coordinate hyperplane. Incoming corrected THM-4727 is retained intact.
- HC(4)=>JC(2) is the cited direction; a converse is not supplied.
- phi(m)/2 is an iid model mean, not a proven first-prime waiting law.
  CRT controls the multiple clock, not a normalized prime-index extremum.

## Concurrent priority and additional contribution

Inherited base was5e1cedf452c7. A later fetch brought in b687cd28ca and
7b76044ca1 (Session144),5e14c36ccf (weighted-family fibre/inertia extension),
and2a836157af (Session145 parity carrier and exact571 minimum).
Their overlapping results have priority and are explicitly credited.

This packet adds a short bank of571 ACTUAL S/B words of depth at most5,
verified by an independent full evaluation-matrix rank. The incoming
Session145 certificate instead used exterior-coefficient contexts, each
representable by a potentially large signed sum of actual queries. The
uniform path-cover recurrence/real-rooted aggregate/variance obstruction,
one-cut second-moment witness, and finite cyclic waiting identity are
separate scoped contributions. The Jacobian tower, CRT delay mechanism,
prime-table correction and MRB scope overlap are consolidated confirmations.

No complete Keller classification, JC(2), LRC(14), Collatz convergence,
Fermat-prime finiteness, MRB irrationality, or all-order7/21 spectrum proof
is claimed. ETH is a primary-source scope citation; no Lean build was run.

Artifact and dependency SHA-256 pins are in manifest.json. These
pins normalize CRLF to LF so they survive ordinary Git platform conversion.
The full81/27 cover profiles were independently recomputed after reading
the incoming uniform-prime transport result; the all-order averaged-query
blindness follows from the new recurrence, not a finite word scan.
The maintained
documentation check passed after mainline integration, using the Windows
extended-path prefix because some unrelated repository filenames exceed
the ordinary Windows path length.
