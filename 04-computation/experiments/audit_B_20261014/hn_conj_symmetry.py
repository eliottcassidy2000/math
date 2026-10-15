#!/usr/bin/env python3
"""Audit B: M3 and P1 are invariant under complex conjugation (exact, in the zeta_24 representation), hence
P1 u psi M3 = conj(psi^-1 ... ) is congruent to M3 u psi P1:  conj(psi^-1 (P1 u psi M3)) = psi conj(P1) u conj(M3)
= M3 u psi P1 (|psi| = 1).  So the second 'mixed pair' is the same graph up to an isometry."""
import sys
from importlib import util
sys.argv = ['x', '/dev/null', 'none']
spec = util.spec_from_file_location('rb', 'hn_spindle_free_rebuild_zeta24.py')
rb = util.module_from_spec(spec); spec.loader.exec_module(rb)
for nm, A in (('M3', rb.M3), ('P1', rb.P1), ('M1+M1', rb.M11)):
    C = rb.kconj_arr(A)
    print(f"{nm}: conj(A) == A as sets: {set(map(tuple, C.tolist())) == set(map(tuple, A.tolist()))}", flush=True)
