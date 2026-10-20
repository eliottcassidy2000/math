"""Affine GF(2) compatibility compiler for blue paths and original arc guards.

This solves exact finite labelled orientation families, not tournament classes.
An inconsistent system certifies universal blackness within its stated guards.
A consistent system returns full rows as a hostile, not an all-order theorem.
"""
from itertools import combinations
import probe as P


class Affine:
    def __init__(self):
        self.basis = {}
        self.contradiction = None

    def add(self, mask, value, reason):
        provenance = {reason}
        while mask:
            pivot = mask.bit_length()-1
            if pivot not in self.basis:
                self.basis[pivot] = (mask, value, provenance)
                return
            previous, rhs, reasons = self.basis[pivot]
            mask ^= previous
            value ^= rhs
            provenance ^= reasons
        if value and self.contradiction is None:
            self.contradiction = sorted(provenance)

    def solution(self):
        P.need(self.contradiction is None, 'inconsistent affine system')
        assignment = 0
        for pivot, (mask, value, _) in sorted(self.basis.items()):
            if value ^ ((mask & assignment).bit_count() & 1):
                assignment |= 1 << pivot
        return assignment


def blue_compatible(n, exceptional, words, step, extra_arcs=()):
    """Can a chosen step be blue with all ledger words valid and stated guards?

step is zero based. extra_arcs contains triples(u,v,bit) for oriented u->v.
The returned contradiction reasons use equation labels; a solution includes
the entire labelled source tournament and a literal sequential replay.
"""
    pairs = list(combinations(range(n),2))
    indices = {pair:i for i,pair in enumerate(pairs)}
    affine = Affine()
    def fixed(u,v,bit,reason):
        pair = P.edge(u,v)
        affine.add(1 << indices[pair], bit ^ int(u>v), reason)
    for i in range(n-1):
        fixed(i,i+1,1,f'path:{i}')
    for u,v in pairs:
        reflected = (n-1-v,n-1-u)
        affine.add((1<<indices[(u,v)]) ^ (1<<indices[reflected]),0,f'blueP:{u},{v}')
    y = n-1-exceptional
    for v in range(1,n):
        fixed(0,v,int(v not in {exceptional,n-1}),f'first:{v}')
    for v in range(n-1):
        fixed(n-1,v,int(v in {0,y}),f'last:{v}')
    for i,(u,v,bit) in enumerate(extra_arcs):
        fixed(u,v,bit,f'extra:{i}:{u},{v}')
    compiled = P.compile_ledger(n,words)
    for u,v,bit in compiled['guards']:
        fixed(u,v,bit,f'direction:{u},{v}')
    if compiled['guard_conflicts']:
        return {'compatible':False,'reason':'internally inconsistent ledger directions'}
    prior_edges = set()
    for word in words[:step]:
        prior_edges ^= P.edges(word)
    q = words[step]
    for i,j in combinations(range(n),2):
        u,v = q[i],q[j]
        a,b = q[n-1-j],q[n-1-i]
        pair, image = P.edge(u,v),P.edge(a,b)
        rhs = int(u>v) ^ int(a>b) ^ int(pair in prior_edges) ^ int(image in prior_edges)
        affine.add((1<<indices[pair]) ^ (1<<indices[image]),rhs,f'blueQ:{i},{j}')
    if affine.contradiction is not None:
        return {'compatible':False,'contradiction_equations':affine.contradiction}
    assignment = affine.solution()
    rows = [0]*n
    for i,(u,v) in enumerate(pairs):
        if assignment>>i&1:
            rows[u] |= 1<<v
        else:
            rows[v] |= 1<<u
    rows = tuple(rows)
    P.need(P.path(rows,tuple(range(n))) and P.blue(rows,tuple(range(n))),'affine source model')
    now = rows
    for i,word in enumerate(words):
        P.need(P.path(now,word),'affine ledger direction')
        if i==step:
            P.need(P.blue(now,word),'affine blue model')
        now = P.complement(now,word)
    return {'compatible':True,'rows':rows,'scores':[row.bit_count() for row in rows],
            'strong_partner':P.strong(P.complement(rows,tuple(range(n))))}


def B_words(n, cut=4):
    a,b,x,z = 0,1,2,n-1
    d = tuple(range(3,cut))
    r = tuple(range(cut,n-1))
    P.need(d and r,'B blocks nonempty')
    return ((a,b)+r+(z,x)+d,(b,)+r+(a,z,x)+d,(b,x)+d+r+(a,z))


if __name__ == '__main__':
    import json
    for n in range(6,17):
        words=B_words(n)
        result=blue_compatible(n,n-3,words,1,[(1,n-2,1)])
        print(n,json.dumps(result),flush=True)
