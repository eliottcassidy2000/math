#!/usr/bin/env python3
"""export_automaton.py -- write the condensation (strong-connectivity) automaton of span w for the C
enumerator: nstates naccept / accept list / transition table rows (state) with one entry per letter
in the order of ribbon.letters(w) (-1 = dead)."""
import sys
from condensation import automaton
from ribbon import letters
w = int(sys.argv[1])
order, states, trans, accept = automaton(w)
L = letters(w)
with open(f"span{w}_condensation.txt", "w") as f:
    f.write(f"{len(order)} {len(accept)}\n")
    f.write(" ".join(map(str, accept)) + "\n")
    for i in range(len(order)):
        f.write(" ".join(str(trans[(i, a)]) for a in L) + "\n")
print(f"span {w}: {len(order)} condensation states written")
