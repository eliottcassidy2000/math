# Evidence, reproduction, and validity boundaries

This packet uses separate implementations and exact arithmetic controls.
It does not claim an independent human or subagent audit. No external theorem
ID is newly allocated. The general proofs are stated in COLLATZ.md and
INTERSECTIONS.md, with their classical and inherited inputs named there.

## 1. Potential certificate

Run from this directory with a C++17 compiler and Python3 plus NumPy:

```sh
c++ -std=c++17 -O2 -Wall -Wextra -pedantic lift_queue.cpp -o /tmp/lift-queue
/tmp/lift-queue 1 /tmp/lift-queue-1.bin
/tmp/lift-queue 2 /tmp/lift-queue-2.bin
/tmp/lift-queue 3 /tmp/lift-queue-3.bin
/tmp/lift-queue 4 /tmp/lift-queue-4.bin
/tmp/lift-queue 5 /tmp/lift-queue-5.bin
python3 audit_potential.py 5 /tmp/lift-queue-5.bin --output /tmp/potential5_audit.json
python3 audit_lifts.py --directory /tmp --output /tmp/lift_profiles.json
```

The sparse monotone queue is only a producer. Its termination is not a proof
assumption: the final full inequality scan is the certificate. The independent
auditor uses direct modular exponentiation, integer vector transitions and
divisibility tests rather than the producer's inverse-power recurrence or
digit-peeling reward. Compact state indexing is specified in both sources.
All states, including unreachable concerns one might otherwise filter away,
are checked. The unbounded-letter tail is proved by the integer cost bound.

Receipts: [queue4.json](queue4.json), [queue5.json](queue5.json),
[potential5_audit.json](potential5_audit.json), [lift_profiles.json](lift_profiles.json).
The E=5 byte hash is
`f2522660b63908ee3b9fe414d2930efcba2c0c1cb7dd9a54c1d8314c6806dcca`.
The O0 and O2 builds produce the same31,250,000 bytes. The E=4 vector agrees
entrywise with the preceding packet's certificate. Files are explicit
little-endian signed32-bit values in unit-multiplier row order, then b order.

Positive control: all1250 letters1 attain equality, with reward781. Negative
controls: multiplier-only word(3,1^16) has spurious truncated reward12; the
actual affine content reward is0. Deleting old slack from the chosen lift
fails at the explicit state in COLLATZ.md§2.

## 2. Elementary theorem controls

```sh
python3 defect_controls.py
python3 periodic_controls.py
python3 cyclic_intersections.py
python3 epoch_budget.py
python3 -O defect_controls.py --output /tmp/defects-O.json
python3 -O periodic_controls.py --output /tmp/periodic-O.json
python3 -O cyclic_intersections.py --output /tmp/intersections-O.json
python3 -O epoch_budget.py --output /tmp/epoch-budget-O.json
cmp defects.json /tmp/defects-O.json
cmp periodic.json /tmp/periodic-O.json
cmp intersections.json /tmp/intersections-O.json
cmp epoch_budget.json /tmp/epoch-budget-O.json
```

- Defect controls: every nonempty positive composition with total clock<=16,
  65,535 words, with1,148,392 checks. Complete integer coefficients verify the
  untruncated inequality. The dark-height family is checked integrally
  through requested height5 and modularly through7. No infinite word is
  inferred from the finite samples; the lifting argument is in the proof.
- Periodic controls: all1,336 primitive words of lengths1–5 over{1,2,3,4},
  all4,932 proper phases, repeats0–20, and boundary valuations at5,7,11.
  In addition2,520 marked-entry cases include96 exact resonances; repeats
  through5 test the nonresonant divisibility. There are367,923 checks.
  Full direct integer iteration supplies the comparator for the symbolic
  geometric-sum and commutator formulas.
- Tournament controls: every cyclic subgroup at3,5,7, all constant-step
  charts there, and every legal path at5; composite hosts9,15,21,27.
  There are35,483 checks, including independent adjacency invariance and
  edge-orbit counting. The prime theorem uses the cited group theorem;
  the composite all-size claim uses the explicit fibre translations.
- The refined epoch-budget hypothesis is tested on all262,143 nonempty
  positive words with clock<=18. The separate finite searches truncate
  height at E=1,2,3 and allow letters<=12; their precision is part of the
  result. This is evidence for an OPEN sufficient target, not an unbounded
  certificate. Explicit runtime guards raise errors rather than certify
  interrupted searches.

Python checks use explicit exceptions, so optimization does not erase them.
No floating-point estimates occur. The mathematical review checked the
chronological order in `vu` versus `uv`, nonprimitive resonance, equal ports
in the composite hostile, owner versus pair multiplicity, the strict range
in the CRT zero test, and necessary-versus-sufficient directions.

## 3. Source and concurrent-work review

Started from clean b51ef5cd7393 and fetched the three incoming Session112
commits through b7285c3d9a. Read their maintained README, canonical-count
and bank-ranking interfaces, paired aliquot decoder, and paper-scope note;
the new intersections specialize their counting interface without claiming
an unrestricted running-time bound. The paper was read as text and rendered
pages1–6. The prime-degree input was checked directly in Dobson–Witte
Theorem1. No factoring/ETH Lean build was run here.

The known unsupported equivalence in the previous Collatz packet was repaired
in place and logged in01-canon/MISTAKES.md. No finite certificate was retracted.
The broad number connections are classified in CONNECTIONS.md, with stronger
older results cited before the weaker new specializations.

## 4. Remaining obligations

Neither a fifth layer nor a finite-budget search proves all-layer feasibility.
Fixed-controller content divisors do not bound a growing entry word for free.
Prime-parent singleton intersections do not bound the number of parents or
native charts needed for every tournament, and fail at composite order by
an exponential family. These are explicit open obligations, not hidden
steps in a universal Collatz or unrestricted enumeration proof.
