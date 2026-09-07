"""Exact two-diagonal deletion potential; no repository imports or float gates."""
from collections import Counter, deque
from fractions import Fraction
from itertools import combinations, product, permutations
from math import factorial, comb
from pathlib import Path
import json
import sys
sys.stdout.reconfigure(encoding='utf-8', newline='\n')

gates = 0
def check(condition, label):
    global gates
    gates += 1
    if not condition:
        raise RuntimeError(label)

def boards(n):
    pairs = tuple(combinations(range(n), 2))
    degrees = [0]*n
    rows = []
    def visit(r):
        if r == n:
            if all(v == 2 for v in degrees):
                yield tuple(rows)
            return
        for a,b in pairs:
            if degrees[a] == 2 or degrees[b] == 2:
                continue
            degrees[a] += 1; degrees[b] += 1
            if all(2-v <= n-r-1 for v in degrees):
                rows.append((a,b))
                yield from visit(r+1)
                rows.pop()
            degrees[a] -= 1; degrees[b] -= 1
    yield from visit(0)

def potential(points):
    ds = sorted({c-r for r,c in points})
    ss = sorted({c+r for r,c in points})
    di = {d:i+1 for i,d in enumerate(ds)}
    si = {s:i+1+len(ds) for i,s in enumerate(ss)}
    target = len(ds)+len(ss)+1
    adj = [[] for _ in range(target+1)]
    def edge(u,v,cap):
        adj[u].append([v,cap,len(adj[v])])
        adj[v].append([u,0,len(adj[u])-1])
    for d in ds: edge(0,di[d],2)
    for s in ss: edge(si[s],target,2)
    for r,c in points: edge(di[c-r],si[c+r],1)
    flow = 0
    while True:
        parent = [None]*(target+1); parent[0] = (0,0)
        todo = deque([0])
        while todo and parent[target] is None:
            u = todo.popleft()
            for i,(v,cap,_) in enumerate(adj[u]):
                if cap and parent[v] is None:
                    parent[v] = (u,i); todo.append(v)
        if parent[target] is None: break
        v = target
        while v:
            u,i = parent[v]
            e = adj[u][i]; e[1] -= 1
            adj[v][e[2]][1] += 1
            v = u
        flow += 1
    return len(points)-flow

def certificates(points):
    d = Counter(c-r for r,c in points)
    s = Counter(c+r for r,c in points)
    z = Counter(c+r for r,c in points if d[c-r] <= 2)
    iso = Counter(c+r for r,c in points if d[c-r] == 1)
    f = sum(max(0,v-2) for v in d.values())
    g = sum(max(0,v-2) for v in z.values())
    i = sum(v == 3 and iso[k] == 3 for k,v in s.items())
    return f,g,i

def brute(points):
    for k in range(len(points)+1):
        for lost in combinations(range(len(points)),k):
            removed = set(lost)
            d = Counter(c-r for j,(r,c) in enumerate(points) if j not in removed)
            s = Counter(c+r for j,(r,c) in enumerate(points) if j not in removed)
            if max(d.values(),default=0)<=2 and max(s.values(),default=0)<=2:
                return k
    raise RuntimeError('empty subset must work')

def eligible(rows, selected):
    n = len(rows)
    other = [[] for _ in range(n)]
    for r,pair in enumerate(rows):
        for c in pair: other[c].append(r)
    for r in selected:
        for c in rows[r]:
            t = next(v for v in other[c] if v != r)
            if any(j != r and (j == t or 2*j == r+t) for j in selected):
                return False
    return True

def conditional_control(rows, s, selected):
    n = len(rows)
    ds = {s-2*r for r in selected}
    targets = [s-r for r in selected]
    check(len(ds)==3 and len(set(targets))==3, 'three distinct directions/labels')
    assignments = 0
    for chosen in product(*(rows[r] for r in selected)):
        check(len(set(chosen))==3, 'eligible source columns distinct')
        fixed = dict(zip(chosen,targets))
        forced = [(r,fixed[c]) for r,pair in enumerate(rows) for c in pair if c in fixed]
        inside = [(r,c) for r,c in forced if r+c==s or c-r in ds]
        check(set(inside)==set(zip(selected,targets)), 'all forced companions outside union')
        free = [c for c in range(n) if c not in chosen]
        labels = [c for c in range(n) if c not in targets]
        forbidden = {(c,t) for c in free for t in labels
                     if any(c in pair and (r+t==s or t-r in ds)
                            for r,pair in enumerate(rows))}
        check(max(Counter(c for c,t in forbidden).values(),default=0)<=8,
              'source forbidden degree')
        check(max(Counter(t for c,t in forbidden).values(),default=0)<=8,
              'target forbidden degree')
        lengths = sum(1 for r in range(n) if 0<=s-r<n)+sum(n-abs(d) for d in ds)
        check(len(forbidden)<=2*lengths, 'retained four-line count bound')
        assignments += 1
    check(assignments==8, 'all eight assignments')

def avoidance_controls():
    count = 0
    for n in range(2,8):
        masks = [set(), {(r,r) for r in range(n)},
                 {(r,c) for r in range(n) for c in range(n) if (c-r)%n in (0,1)},
                 {(r,c) for r in range(n) for c in range(n) if r==0 or c==0},
                 {(r,c) for r in range(n) for c in range(n)}]
        for forbidden in masks:
            m = len(forbidden)
            delta = max([0]+list(Counter(r for r,c in forbidden).values())+
                        list(Counter(c for r,c in forbidden).values()))
            hits = Counter(sum((r,p[r]) in forbidden for r in range(n))
                           for p in permutations(range(n)))
            moments = []
            for k in range(n+1):
                moment = Fraction(sum(comb(w,k)*ct for w,ct in hits.items()),factorial(n))
                ordered = moment*factorial(k)*(factorial(n)//factorial(n-k))
                check(ordered.denominator==1, 'integer ordered rook count')
                if k>=2:
                    check(0<=m**k-ordered<=k*(k-1)*delta*m**(k-1), 'uniform collision bound')
                moments.append(moment)
                partial = sum((-1)**j*moments[j] for j in range(k+1))
                exact = Fraction(hits[0],factorial(n))
                check(partial<=exact if k%2 else partial>=exact, 'Bonferroni parity')
            count += 1
    return count

def main():
    census = []
    positive = None; strict = None; full_hostile = None
    for n,expected in zip(range(2,7),(1,6,90,2040,67950)):
        total = 0; sums = [0,0,0,0]; hist = Counter(); controls = 0
        for rows in boards(n):
            points = tuple((r,c) for r,pair in enumerate(rows) for c in pair)
            tau = potential(points); f,g,i = certificates(points)
            check(tau>=f+g>=f+i, 'disjoint dual deletion certificate')
            if n<=5:
                check(tau==brute(points), 'independent exhaustive subset deletion')
            if i and positive is None: positive = dict(n=n,rows=rows,tau=tau,F=f,G=g,I=i)
            if tau>f+g and strict is None: strict = dict(n=n,rows=rows,tau=tau,F=f,G=g,I=i)
            if tau==0 and full_hostile is None:
                bad = next(((a,b,c) for a,b,c in combinations(points,3)
                            if (b[0]-a[0])*(c[1]-a[1])==(c[0]-a[0])*(b[1]-a[1])),None)
                if bad: full_hostile = dict(n=n,rows=rows,collinear_triple=bad)
            total += 1; hist[tau] += 1
            sums = [a+b for a,b in zip(sums,(tau,f,g,i))]
            # Complete one-transposition universe through n=4; deterministic sample at n=5,6.
            if n<=4 or total%997==1:
                for a,b in combinations(range(n),2):
                    moved = tuple((r,b if c==a else a if c==b else c) for r,c in points)
                    t2 = potential(moved)
                    k = len(set(points)-set(moved))
                    check(abs(tau-t2)<=k<=4, 'retained-point and column-transposition Lipschitz')
                    controls += 1
        check(total==expected, 'complete saturated universe')
        row = dict(n=n,boards=total,mean_tau=str(Fraction(sums[0],total)),
                   mean_F=str(Fraction(sums[1],total)),mean_G=str(Fraction(sums[2],total)),
                   mean_I=str(Fraction(sums[3],total)),tau_histogram=dict(sorted(hist.items())),
                   transposition_controls=controls)
        census.append(row)
        print(json.dumps(row,sort_keys=True))
    check(positive is not None, 'isolated-line positive control')
    check(strict is not None, 'dual certificate need not be equality')
    check(full_hostile is not None, 'two directions do not imply full success')
    check(potential(((0,0),))==0, 'one-cell capacity-one hostile')
    # Complete eligible-row-triple test on the cyclic skeletons n=6..12.
    eligible_counts = []
    for n in range(6,13):
        rows = tuple((r,(r+1)%n) for r in range(n))
        count = 0
        for s in range(2*n-1):
            rs = [r for r in range(n) if 0<=s-r<n]
            bad = 0
            for triple in combinations(rs,3):
                if eligible(rows,triple):
                    conditional_control(rows,s,triple); count += 1
                else: bad += 1
            check(bad<=4*len(rs)*max(0,len(rs)-2), 'uniform ineligible triple bound')
        eligible_counts.append([n,count])
    check(all(c>0 for n,c in eligible_counts), 'nonvacuous conditional controls')
    avoidance = avoidance_controls()
    # Exact exp(2) enclosure; signs in the new coefficient are retained.
    term = Fraction(1); lo = term
    for k in range(1,61):
        term *= Fraction(2,k); lo += term
    hi = lo + term*Fraction(2,61)/(1-Fraction(2,62))
    qlo,qhi = 1/hi,1/lo
    elo = qlo/12-qhi**2/2+5*qlo**3/4+qlo**4/6
    ehi = qhi/12-qlo**2/2+5*qhi**3/4+qhi**4/6
    glo,ghi = 1-5*qhi+elo,1-5*qlo+ehi
    check(elo>Fraction(1,200), 'strict improved constant')
    check(glo>Fraction(328,1000) and ghi<Fraction(329,1000), 'gamma rational enclosure')
    # Formal exponential antiderivative coefficients at 1 minus at 0.
    anti = {4:Fraction(1,4),2:Fraction(-3,2),0:Fraction(15,4),-2:Fraction(1,2)}
    expected = {-2:Fraction(1,12),-4:Fraction(-1,2),-6:Fraction(5,4),-8:Fraction(1,6)}
    check({k-6:v/3 for k,v in anti.items()}==expected, 'exact integral coefficients')
    result = dict(census=census,positive=positive,strict_dual_hostile=strict,
                  two_direction_zero_hostile=full_hostile,
                  eligible_cyclic_triples=eligible_counts,avoidance_matrices=avoidance,
                  epsilon_interval=[str(elo),str(ehi)],gamma_interval=[str(glo),str(ghi)],
                  status='FINITE-EXACT; asymptotic theorem is proved in companion prose')
    path = Path(__file__).with_name('continuing10_20260907_no3_repair_certificate.json')
    path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print('Positive isolated-line control:',json.dumps(positive,sort_keys=True))
    print('Strict dual hostile:',json.dumps(strict,sort_keys=True))
    print('Two-direction zero hostile:',json.dumps(full_hostile,sort_keys=True))
    print('Eligible cyclic triples:',json.dumps(eligible_counts))
    print('Avoidance matrices:',avoidance)
    print('Exact interval checks: epsilon>1/200; 328/1000<gamma<329/1000')
    print('Always-active exact gates:',gates)

if __name__=='__main__': main()
