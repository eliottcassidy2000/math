"""Independent standard-library controls for the session198 source guardrail.

Uses no inherited source oracle, scalar implementation, subset counter, or
external tournament bank. Finite controls do not replace the written proofs.
"""
from functools import lru_cache
from itertools import permutations
import json
from pathlib import Path


def require(test, message):
    if not test:
        raise RuntimeError(message)


@lru_cache(None)
def constants(L):
    if L == 1:
        return (1, 1, 1)
    ell = (L+1)//2
    g, q, b = constants(ell)
    if L % 2 == 0:
        return (g, q, 2*b)
    return (g+ell-ell%2, q+b, 2*b)


def sources(cap):
    origin = (1, 1, 1, 1)
    result = {origin: ()}
    pending = [origin]
    for v in pending:
        D,E,U,d = v
        moves = [(('O', None), (D,0,2*U,d))]
        for L in range(1, cap//U+1):
            g,Q,B = constants(L)
            moves.append((('M',L), (g*D-(L%2)*E+L*U,
                                       D-int(L==1)*E,Q*D+B*U,D)))
        for move, w in moves:
            if w[2] <= cap and w not in result:
                result[w] = result[v] + (move,)
                pending.append(w)
    return result


def N(v,w):
    D,E,U,d = v
    g,e,L,f = w
    return (D-E+U)*g-(D-E+d)*e+(D+3*U-d)*L-(U+E)*f


def M(v,w):
    D,E,U,d = v
    g,e,L,f = w
    return (D+2*U)*g-(U+d)*e+(2*D-E+6*U-d)*L-(U+E)*f


def paths(rows, vertices):
    return [p for p in permutations(vertices)
            if all(rows[a]>>b&1 for a,b in zip(p,p[1:]))]


def rooted_profile(rows, root):
    vertices = tuple(range(len(rows)))
    hp = paths(rows, vertices)
    D,E = len(hp),sum(p[-1]==root for p in hp)
    s = sum(p[0]==root for p in hp)
    d = len(paths(rows, tuple(v for v in vertices if v!=root)))
    U = 0
    for mask in range(1 << len(rows)):
        if mask>>root&1:
            qv = tuple(v for v in vertices if mask>>v&1)
            pv = tuple(v for v in vertices if not mask>>v&1)
            U += len(paths(rows,pv))*sum(q[0]==root for q in paths(rows,qv))
    return (D,E,U,d),s


def audit():
    origin = (1,1,1,1)
    v,w = (3,1,3,1),(1,0,6,1)
    for D,E,U,d in (v,w):
        require(1<=d<=D<=U and 0<=E<=D, 'basic inequalities')
        require(D-E>=d and U>=D+d-E and U>=2*d, 'strong inequalities')
        require(E<=d, 'endpoint deletion inequality')
        require(D%2==1 and d%2==1 and (U-E)%2==0, 'parity controls')
    require(N(v,w)==67, 'formal target67')
    bank = sources(11)
    require(len(bank)==15 and v in bank and w not in bank,'physical source bank')
    actual = sorted(set([15]+[form(a,b) for a in bank for b in bank
        if a!=origin and b!=origin for form in (N,M) if form(a,b)<=67]))
    require(67 not in actual, 'all-order C2 inverse cutoff rejects67')
    triangle_profile,s = rooted_profile([2,4,1],0)
    require(triangle_profile==v and s==1,'literal triangle source')
    for n in range(2,9):
        rows = [sum(1<<j for j in range(i+1,n)) for i in range(n)]
        profile,s = rooted_profile(rows,0)
        require(profile==(1,0,2**(n-1),1) and s==1,'transitive source ray')
    rows = [70,124,120,97,73,17,32]
    n = len(rows)
    require(all(not(row>>i&1) for i,row in enumerate(rows)), 'no loops')
    require(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1
                for i in range(n) for j in range(i)), 'tournament orientation')
    H = len(paths(rows,tuple(range(n))))
    C = sum(all(rows[a]>>b&1 for a,b in zip((0,)+p,p+(0,)))
            for p in permutations(range(1,n)))
    require((H,C)==(67,3), 'literal independent C3 witness')
    return dict(status='PASS',source_cap=11,source_states=len(bank),
                formal_C2_target=67,C2_values_through67=actual,
                actual_C3_rows=rows,actual_H=H,actual_C=C,
                transitive_source_orders=list(range(2,9)),
                independence='No imported packet functions; direct permutation counters.')


if __name__=='__main__':
    report=audit()
    Path(__file__).with_name('GUARDRAIL_AUDIT.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,sort_keys=True))
