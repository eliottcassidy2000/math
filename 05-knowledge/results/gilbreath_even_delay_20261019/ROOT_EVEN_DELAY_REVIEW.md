# Independent review: all finite size-six delays are even

Status: PASS, 2026-10-19.

The root reviewer read the complete pair-language proof. All seven input
pairs and three incoming values agree with the displayed table. The three
forbidden-output arguments each retain the necessary entering state and
correctly force a forbidden input pattern. The first-one formula includes
the corrected index z+j+1; its even/odd splitting proves entry to the
language after any finite initial zero run. Inputs with no one have delay
zero. Thus the proof covers every binary left word, every width, and every
age, not a bounded initial segment.

An independently written monitor/transducer product has exactly seven
reachable states and44 edges for each initial boundary0 or1, and no rejected
output transition. These are exhaustive finite-state checks proving
preservation of the regular language for arbitrary word lengths. The normal
and optimized receipts match byte-for-byte.

The invariant is stronger than the requested no-odd-first-defect statement,
but does not count word multiplicities, bound the first finite delay, or
prove that the alternate branch kills every later defect. Paired halving,
the sharp dyadic horizon, and the exact asymptotic remain open in this packet.

Reproduce:

    python3 verify_pair_language_independent.py > pair_language_independent.json
    python3 -O verify_pair_language_independent.py > pair_language_independent_optimized.json
