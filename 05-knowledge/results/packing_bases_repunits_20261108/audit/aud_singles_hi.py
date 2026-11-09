"""aud_singles_hi.py - print the carry-free singletons a with LO <= a < HI from an aud_fewcycle intercyclic inventory.
usage: python aud_singles_hi.py ic8191.txt LO HI"""
import sys
from collections import Counter

c = Counter(int(line.split()[1]) for line in open(sys.argv[1]))
lo, hi = int(sys.argv[2]), int(sys.argv[3])
for a in sorted(c):
    if c[a] == 1 and lo <= a < hi:
        print(a)
