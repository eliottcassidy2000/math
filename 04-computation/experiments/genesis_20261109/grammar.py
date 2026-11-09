"""Classify tournaments under THM-4942's five-operation grammar (kind-pasteur S17):
   FC   first creator            H(T) not in S(n-1)
   FCV  first creator + vertex   some T - v is a first creator: H(T-v) not in S(n-2)
   SUM  ordinal sum              T not strong
   SUB  substitution             T has a nontrivial module (not prime)
   SUM2 ordinal sum re-connected by <= 2 arcs   lambda(T) <= 2
A tournament matching none of these is an ESCAPER (HYP-9532 (F)).
Uses tprop (exact H, lambda, primality) and the exact spectra S(n), n <= 11.
"""
import subprocess
import sys
import os
from spectra import in_S

HERE = os.path.dirname(os.path.abspath(__file__))


def code_of(n, arc):
    """arc(i, j) -> True if i -> j (i < j)."""
    return ''.join('1' if arc(i, j) else '0' for i in range(n) for j in range(i + 1, n))


def profile(items):
    """items: list of (n, code). Returns list of dicts."""
    inp = ''.join('%d %s\n' % (n, c) for n, c in items)
    r = subprocess.run([os.path.join(HERE, 'tprop')], input=inp, capture_output=True, text=True, check=True)
    res = []
    for line in r.stdout.strip().split('\n'):
        a, b, c = line.split('|')
        n, H, strong, lam, prime, msd = map(int, a.split())
        dels = list(map(int, b.split()))
        res.append(dict(n=n, H=H, strong=strong, lam=lam, prime=prime, msd=msd, dels=dels, code=c.strip()))
    return res


def classify(p):
    """Return (labels, status) where status in {'ESCAPER', 'generated', 'unknown'}."""
    n = p['n']
    labels = []
    unknown = False
    s = in_S(n - 1, p['H'])
    if s is None:
        unknown = True
    elif s is False:
        labels.append('FC')
    fcv = False
    for h in p['dels']:
        t = in_S(n - 2, h)
        if t is None:
            unknown = True
        elif t is False:
            fcv = True
    if fcv:
        labels.append('FCV')
    if not p['strong']:
        labels.append('SUM')
    if not p['prime']:
        labels.append('SUB')
    if p['lam'] <= 2:
        labels.append('SUM2')
    if labels:
        return labels, 'generated'
    return labels, ('unknown' if unknown else 'ESCAPER')


def band(n, s):
    """back arcs j -> i exactly when 1 <= j - i <= s; all longer arcs forward."""
    return code_of(n, lambda i, j: (j - i) > s)


if __name__ == '__main__':
    items = []
    for s in (2, 3, 4, 5):
        for n in range(2 * s + 2, 15):
            items.append((n, band(n, s)))
    for p in profile(items):
        labels, st = classify(p)
        s = None
        print('n=%2d H=%7d strong=%d lam=%d prime=%d msd=%d  dels=[%d..%d]  %-14s %s' % (
            p['n'], p['H'], p['strong'], p['lam'], p['prime'], p['msd'], min(p['dels']), max(p['dels']),
            ','.join(labels) or '-', st))
