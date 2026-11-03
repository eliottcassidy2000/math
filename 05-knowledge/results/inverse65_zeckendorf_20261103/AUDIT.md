# Verification and scope audit

Status: same-agent proof audit plus independent exact algorithms; no independent
human or separate-agent review is claimed. All new Python validity checks use
explicit exceptions and survive `-O`.

## Reproduction and input identity

Run `python -X utf8 -B reproduce.py` from this packet or by its repository path.
It checks68 inherited file hashes, then runs each of the five probes normally
and with `-O`. Each pair must produce identical LF-normalized JSON bytes.
`reproduction.json` records commands, timings, outputs and digests. Hashes
normalize CRLF to LF to allow ordinary Windows/Git line-ending conversion.
The inherited atlas is pinned at the source commit, rather than silently
reselected from a future search.

The independent C++ zero-cycle check is:

```text
g++ -O2 -std=c++17 zero_cycle.cpp -o zero_cycle
./zero_cycle units.txt
```

`zero_cycle_cpp.out` records960799 words, one trace candidate, zero matches.
The Python proof also checks the64-bit matrix-entry range for this experiment;
the second characteristic coefficient uses128-bit intermediate arithmetic.
The pruned SymPy certificate checks2,8,53 admissible words independently.
No executable is committed.

## Declared universes and controls

| Probe | Universe and independent path | Result |
|---|---|---|
| Inverse | All747 reachable states h<=151; compare bounded closure with adjugate recursion |All decode|
| Negative image | All physical-cone triples with odd h<=15 and x,y<=2h |942 absent states reject|
| Norm / cone | All65 matrices and stated extreme rays / positive vectors |9/5 norm growth, monotone h,16h bound|
| aa offsets | All states h<=192 plus all278915 words of lengths1..3 |41 states, maximum h155|
| H55 sources | All334 words h<=55, literal subset DP on every cut of all six H55 graphs |768 failures to reach111; separate known source succeeds|
| Charge classification | Both signs, rational nullspaces of all65 matrices |Seven letters, five charge lines|
| Prime and7p examples | Independently reconstruct and count literal target rows |499,1009,19861,21001 and91,147007,139027 agree|
| Quadratic clocks | Complete pair orbits for listed moduli,13^1..13^3 and mixed CRT controls |All-e proof agrees with every finite check|
| Literal ear lift | All30 code-8 ears, every Hamilton path, exact full-row equality after relabelling |28 lift; two H41 fail; one extra atom repairs both|
| Zero-cycle bound | All seven-unit-letter words length<=7, independently norm-pruned exact characteristic polynomials |No q3,6,9 cycle|
|81/27 | Actual12-vertex graphs, automorphism enumeration and named-port path DP |Same67 cyclic triples and H; port response separates|

No program samples arbitrary real weights as a substitute for Boolean edges.
There is no assumption that equal H, equal cyclic triples, or a positive cone
state gives the same source. Source reconstruction is tested at the adjacency
level in the ear lift and on the large quadratic examples.

## Proof audit and boundaries

1. **I2 inverse:** root is a formal unit; actual scalar H1 uses singleton a.
   Every inverse branch is checked, with a decreasing positive integer norm.
   No least-predecessor heuristic is treated as complete. Finite runtime follows
   mathematically; large scalar searches can still be expensive.
2. **I3 branch deficit:** inequalities hold only for the stated same-end bank
   operation. Re-marking, a second end, arbitrary ears, and replacing the parent
   are outside that theorem. The inherited9-dimensional both-end observer is
   a candidate enlargement, not already a solution.
3. **I4 exceptions:** long words have a nonempty prefix before the qualifying
   suffix. Short words are exhausted separately; a finite cutoff is derived
   before the search. The maximum155 is therefore all-order for aa, not merely
   the largest observed example.
4. **C1 charges:** classified single letters only. Products may possess other
   charges. The sign phase is retained; no scalar second-order recurrence is
   asserted when its forcing has not been removed.
5. **C4 clocks:** the norm13 recurrence is singular only at13. Pair equality
   is tested before projecting to one residue; this is essential to the exact
   transient proof. No numerical coincidence is used to infer primality.
6. **Z3-Z4:** the finite transducer is deterministic, appends literal words,
   and allows finite initial/terminal output. Its zero cycle must realize a
   three-root scalar sequence, so the characteristic-polynomial restriction is
   necessary. Passing that test would not itself give a compiler.
7. **E2 lift:** every Hamilton marking is tried, so failed children are truly
   outside the declared bank on their existing vertices. The repair is minimal
   for these two children only; no minimal all-ear-closed alphabet is claimed.
8. **E3 series:** coefficients count presentations. Neither weighted H totals
   nor presentation counts are substituted for isomorphism classes or scalar
   support. A000568 recreation is not advanced by relabelling this series.

The full required prime/7p inverse, unrestricted Hamilton spectrum, arbitrary
re-marking closure and all-order smallest attachment carrier remain **OPEN**.
