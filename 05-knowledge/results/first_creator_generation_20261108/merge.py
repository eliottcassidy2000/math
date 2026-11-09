import sys, struct
from array import array
HCAP = 100000
def load(prefixes):
    A = [0] * HCAP; F = [0] * HCAP
    for p in prefixes:
        a = array('Q'); 
        with open(p + '.cnt', 'rb') as fh: a.fromfile(fh, 2 * HCAP)
        for i in range(HCAP):
            if a[i]: A[i] += a[i]
            if a[HCAP + i]: F[i] += a[HCAP + i]
    return A, F
def readset(fn): return set(int(x) for x in open(fn).read().split())
if __name__ == '__main__':
    n = int(sys.argv[1]); Sn = readset(sys.argv[2]); Snext = readset(sys.argv[3]) if sys.argv[3] != '-' else None
    A, F = load(sys.argv[4:])
    got = {h for h in range(HCAP) if A[h]}
    print('order %d -> %d: values attained by extensions %d (max %d)' % (n, n + 1, len(got), max(got)))
    if Snext is not None:
        print('  matches S(%d) from the class file: %s (missing %d, extra %d)' % (n + 1, got == Snext, len(Snext - got), len(got - Snext)))
    new = sorted(got - Sn)
    unc = [h for h in new if F[h] == 0]
    print('  new values at order %d: %d; not reached from a first creator of order %d: %d %s' % (n + 1, len(new), n, len(unc), unc[:40]))
    old_only_fc = [h for h in new if F[h] == A[h]]
    print('  new values reached ONLY from first-creator parents: %d' % len(old_only_fc))
