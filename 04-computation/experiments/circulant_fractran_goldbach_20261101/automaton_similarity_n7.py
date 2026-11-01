#!/usr/bin/env python3
"""mac-mini-2026-11-01-lenses: GF(2)-similarity types of 7-tournament adjacency matrices versus the 272 merged nodes
(THM-643 census: 456 classes, 88 self-converse, so 88 + 368/2 = 272 merged nodes).

Every class has a Hamiltonian path, so the 2^15 tilings of the fixed path realise every class; similarity types are
counted over all tilings (no canonicalisation needed, since similarity is an isomorphism and converse invariant).
PRE-REGISTERED (before the first run): S4  fewer than 272 similarity types at n = 7 (the automaton does not classify
the merged nodes at n = 7).
usage: python3 automaton_similarity_n7.py
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('ast', os.path.join(HERE, 'automaton_similarity_types.py'))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def main():
    n = 7
    IRR = m.irreducibles(7)
    types = set()
    for A in m.tilings(n):
        rows = [sum(1 << v for v in range(n) if A[u][v]) for u in range(n)]
        types.add(m.similarity_type(rows, n, IRR))
    ok = len(types) < 272
    print(f'  n=7: {len(types)} GF(2)-similarity types over all 32768 tilings; merged nodes 272')
    print(f'  [{"PASS" if ok else "FAIL"}] S4: fewer than 272 similarity types at n = 7')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
