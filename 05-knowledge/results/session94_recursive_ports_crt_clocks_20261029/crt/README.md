# Session94 CRT / reflection carrier lane

Source repository pin: `88c969273f2b6bc4ea2a8c5d29bb36c558eb59d2`.

Read:
1. `COMPATIBLE_CRT_CARRIERS.md`: all-odd ring normal form; minimal dihedral
   reflected double; exact addition obstructions; translation subgroup law.
2. `NINE_CARRIES_AND_ALTERNATING_GROUPS.md`: exact nine-carry inverse reader;
   sparse primitivity certificates at periods25/55; CFSG-dependent alternating
   group conclusions. The proof explicitly charges algebraic operations,
   not constant-cost big-integer arithmetic.

Standard-library reproduction, no network or external package required:

```
python3 check_crt.py
python3 verify_sparse_certificate.py
python3 fast_reader.py
python3 verify_period55.py
python3 check_jones_degrees.py
```

Every command was rerun with `python3 -O`; paired JSON receipts match.
`check_crt.py` exhausts words only at9,15,21,25. The period55 certificate
verifier does not import the decoder and checks just44 word occurrences.
`check_jones_degrees.py` verifies degree exclusions and signs, not Jones's
classification, whose primary source and CFSG dependency are in the proof.

The general compatible-chart and nine-carry theorems have proofs. The
specific group identifications depend on complete finite certificates,
the inherited all-odd reader theorem, and cited classification. No claim
of universal composite primitivity or resolution of Collatz is made.
