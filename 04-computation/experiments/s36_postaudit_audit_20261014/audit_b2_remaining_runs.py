"""Referee check (B, completion): the session runs of the near-transitive sweep that audit_b_near_transitive_dfs.py did
not reproduce: (8,2), (9,2), (15,2), (8,3), (9,3), (13,3), (14,3), (8,4), (9,4), (12,4).  Same independent code (explicit
DFS, Landau + Warshall strongness, triangle test for M) imported from audit_b_near_transitive_dfs.py; nothing new.
Usage: python audit_b2_remaining_runs.py n k [n k ...]   (each invocation kept under ~3 minutes)
"""
import sys
from audit_b_near_transitive_dfs import run

if __name__ == '__main__':
    a = [int(x) for x in sys.argv[1:]]
    for n, k in zip(a[0::2], a[1::2]):
        run(n, k)
