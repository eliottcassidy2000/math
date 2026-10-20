#!/usr/bin/env python3
"""Exact log intervals and all neighbor-disjoint candidates below weight 8000."""
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path
import json

checks = 0


def need(test, message):
    global checks
    checks += 1
    if not test:
        raise ArithmeticError(message)


def logarithm(t, terms=80):
    lower = 2*sum((t**(2*j+1)/F(2*j+1)
                   for j in range(terms)), F())
    upper = lower+2*t**(2*terms+1)/((2*terms+1)*(1-t*t))
    return lower, upper


L2, U2 = logarithm(F(1, 3))
L3, U3 = logarithm(F(1, 2))
need(L2 > F(2, 3), 'log2 strict lower bound')
N = 2**26
need(N == 8192**2 and 7999**2 < N, 'Legendre magnitude/weight relation')
need(F(3, 2*(3*N+2)) < F(1, 2*N), 'negative integer-minimum bound')
lo, hi = L3/U2, U3/L2
hprev, h, kprev, k = 0, 1, 1, 0
convergents = []
prefix = []
for step in range(50):
    a = lo.numerator//lo.denominator
    need(a == hi.numerator//hi.denominator, 'ambiguous partial quotient')
    prefix.append(a)
    hprev, h = h, a*h+hprev
    kprev, k = k, a*k+kprev
    q, p = h, k
    if p >= 8000:
        terminal = {'q': q, 'p': p}
        break
    if p >= 2 and p < q < 2*p:
        convergents.append((p, q))
    lo, hi = 1/(hi-a), 1/(lo-a)
else:
    raise ArithmeticError('no terminal denominator')
need(convergents == [(2, 3), (5, 8), (12, 19), (41, 65),
                     (53, 84), (306, 485), (665, 1054)],
     'small convergent universe')
need(terminal == {'q': 24727, 'p': 15601}, 'terminal convergent')


def direct_carry(word):
    b = 0
    for i, bit in enumerate(word):
        b = (3 if bit else 1)*b+bit*2**i
    positions = [i for i, bit in enumerate(word) if bit]
    positional = sum(2**i*3**(len(positions)-j-1)
                     for j, i in enumerate(positions))
    need(b == positional, 'carry algorithms disagree')
    return b


def path(source, edges):
    word = source[:]
    for i in edges:
        j = (i+1) % len(word)
        need(word[i] != word[j], 'nonactual edge')
        word[i], word[j] = word[j], word[i]
    return word


def all_actual_neighbors(source):
    result=set(); q=len(source)
    positions=[i for i,bit in enumerate(source) if bit]
    for label,i in enumerate(positions):
        j=positions[(label+1)%len(positions)]
        for direction in (-1,1):
            edges=[(i if direction==1 else i-1)%q,(j if direction==1 else j-1)%q]
            if {edges[0],(edges[0]+1)%q}&{edges[1],(edges[1]+1)%q}: continue
            if all(source[e]!=source[(e+1)%q] for e in edges):
                result.add(tuple(edges)); result.add(tuple(edges[::-1]))
    return result


rows = []
hits = []
trace = sha256()
# The known -17 cycle is not subject to the new-cycle minimum bound.
shapes = sorted(set(convergents+[(7, 11)]))
for p, q in shapes:
    r = q-p
    need(gcd(p, r) == 1 and 1 <= r < p, 'candidate shape hypothesis')
    D = 2**q-3**p
    source = [((i+1)*p)//q-(i*p)//q for i in range(q)]
    labels = {(-j*r) % p: (j+(j*r+p-1)//p-1) % q for j in range(p)}
    predicted = set()
    count = 0
    maximum_content = 0
    n=p-r
    for family,residues in [('plus',range(max(0,r-n))),('minus',range(n,r))]:
        for u in residues:
            base_edges=[labels[u],labels[u+n]]
            if family=='minus': base_edges=[(edge-1)%q for edge in base_edges]
            for edges in (base_edges,base_edges[::-1]):
                predicted.add(tuple(edges))
                word=path(source,edges)
                B=direct_carry(word); content=gcd(abs(D),B)
                maximum_content=max(maximum_content,content)
                record={'p':p,'q':q,'r':r,'family':family,'u':u,'edges':edges,'B':B,'D':D,'content':content,'word':''.join(map(str,word))}
                trace.update(json.dumps(record,sort_keys=True,separators=(',',':')).encode()+b'\n')
                if content==abs(D): hits.append(record)
                need(content<abs(D),'integral neighbor-disjoint hostile')
                count+=1
    need(predicted == all_actual_neighbors(source),
         'candidate legality not exhaustive')
    rows.append({'p': p, 'q': q, 'r': r, 'clock': D, 'paths': count,
                 'maximum_content': maximum_content,
                 'reason_included': 'known negative 17-cycle' if (p, q) == (7, 11)
                                    else 'small ordinary convergent'})
total_paths=sum(row['paths'] for row in rows)
need(not hits,'unexpected integer output')
# Positive rational Taylor terms certify the transfer to the audited p/9 gate.
from math import factorial
need(sum((F(29,6)**j/factorial(j) for j in range(21)),F())>121,'exp(29/6)>121')
need(L3>F(13,12),'log3 lower bound')
need(F(13,12)-F(4*29,6*20)==F(7,60)>F(1,9),'large weight gap')
need(8000==20**3,'large weight threshold')
print(json.dumps({
    'status': 'FINITE-EXACT complete candidate check after signed descent reduction',
    'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    'checks': checks, 'log_series_terms': 80,
    'continued_fraction_prefix': prefix, 'terminal_convergent': terminal,
    'small_convergent_shapes': [{'p': p, 'q': q} for p, q in convergents],
    'candidate_shapes': rows, 'total_paths': total_paths, 'integer_hits': hits,
    'all_case_trace_sha256': trace.hexdigest(),
    'dependency': 'The separate complete signed-descent certificate through 2^26 '
                  'is required to make this reduced candidate universe complete.',
    'scope': 'Primitive k=1 same-direction disjoint swaps of consecutive source ones, not arbitrary two edits.'
}, indent=2))
