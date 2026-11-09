"""Decomposability defects of a tournament (genesis lane).
  lam(T)     = fewest arc reversals making T non-strong (arc-strong connectivity).
  defsub(T)  = fewest arc reversals creating a nontrivial module M (2 <= |M| <= n-1):
               min over M of sum_{x not in M} min(|N+(x) & M|, |N-(x) & M|)  (reversals inside M or outside never help).
  twin(T)    = defsub restricted to |M| = 2 (fewest distinguishers of a pair).
  dec(T)     = min(lam, defsub) = fewest reversals reaching a decomposable (non-strong or non-prime) tournament.
               For n >= 3, defsub <= lam always: a non-strong tournament is non-prime, because the side of
               size >= 2 of a dominating cut is a module. So dec = defsub (audit G1, m1).
"""
import sys
from verify_escapers import adj_from, lam

def defsub(out, n, maxsize=None):
    full = (1 << n) - 1
    best = None; bestM = None
    for M in range(1, full):
        sz = bin(M).count('1')
        if sz < 2 or sz > n - 1:
            continue
        c = 0
        for x in range(n):
            if M >> x & 1:
                continue
            k = bin(out[x] & M).count('1')
            c += min(k, sz - k)
            if best is not None and c >= best:
                break
        if best is None or c < best:
            best, bestM = c, M
    return best, bestM

def twin(out, n):
    return min(sum(1 for x in range(n) if x not in (u, w) and ((out[x] >> u) ^ (out[x] >> w)) & 1)
               for u in range(n) for w in range(u + 1, n))

if __name__ == '__main__':
    for line in sys.stdin:
        p = line.split()
        n, code = int(p[0]), p[1]
        out = adj_from(code, n)
        d, M = defsub(out, n)
        L = lam(out)
        print(n, 'lam', L, 'defsub', d, '|M|', bin(M).count('1'), 'twin', twin(out, n), 'dec', min(L, d), ' '.join(p[2:]))
