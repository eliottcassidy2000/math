"""xcheck.py -- the C checker (edimsearch d 0 -x: incremental-free keys, hash test + sort test) against the
pure-Python checker (pycheck.py: histograms as tuples, straight from the definition) on random and on
near-resolving sets.  Compares the defect (#edges - #distinct histograms) of every set.

usage: python3 xcheck.py [nrandom_per_d] [seed]
"""
import random, subprocess, sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import pycheck

def main():
    nr = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 7)
    S19 = [5, 57, 54, 104, 35, 109, 115, 49, 6, 39, 55, 102, 15, 32, 85, 97, 21, 113, 41]
    P6 = [v for v in range(64) if (0x02283022a042a00a >> v) & 1]
    allok = True
    for d in (4, 5, 6, 7):
        sets = []
        if d == 7: sets += [S19] + [S19[:j] + S19[j + 1:] for j in range(19)]          # resolving, and 19 near misses
        if d == 6: sets += [P6] + [P6[:j] + P6[j + 1:] for j in range(15)]             # resolving, and 15 near misses
        for t in range(nr):
            k = rng.randint(1, min(31, 1 << d)); sets.append(rng.sample(range(1 << d), k))
        inp = '\n'.join(' '.join(map(str, S)) for S in sets) + '\n'
        out = subprocess.run([os.path.join(ROOT, 'src', 'edimsearch'), str(d), '0', '-x'], input=inp, capture_output=True, text=True, check=True).stdout.split('\n')
        nres = 0
        for S, line in zip(sets, out):
            dd = pycheck.defect(d, S)
            cdef = int(line.split('defect=')[1])
            if dd != cdef or (line.startswith('RESOLVING') != (dd == 0)): allok = False; print('MISMATCH', d, S, line, dd)
            nres += dd == 0
        print('d=%d: %d sets compared (C check mode vs pure Python), %d resolving' % (d, len(sets), nres), flush=True)
    print('XCHECK ' + ('ALL AGREE' if allok else 'DISAGREEMENT'))

if __name__ == '__main__':
    main()
