# Automata / Ising / numerical-semigroup lane

Read [PROOF.md](PROOF.md). The complete32-attachment central fiber of A^k
is a pinned two-state spin chain. Its normalized H spectrum has an exact
<8,17> description with56 persistent boundary holes from k17; actual H
spacing is2*13^(k−2). No Chan–Pak coverage claim.

Run `python3 verify.py` and `python3 -O verify.py` from any working directory.
The verifier uses only the standard library, recomputes all32 packets by two
methods, reconstructs and trims the carry automaton, checks the all-length
formulas on a declared finite control range, and counts selected literal graphs.
