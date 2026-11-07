#!/usr/bin/env python3
"""Exact controls for loopless digraph cover algebra. Standard library only."""
from itertools import permutations
from math import factorial, isqrt
from collections import Counter
import hashlib, json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def graph(n, mask):
    a = [0]*n
    k = 0
    for u in range(n):
        for v in range(n):
            if u != v:
                if mask >> k & 1:
                    a[u] |= 1 << v
                k += 1
    return a


def complement(a):
    full = (1 << len(a))-1
    return [full ^ (1 << i) ^ a[i] for i in range(len(a))]


def gate(a):
    n = len(a)
    return [x | (1 << n) for x in a] + [1 << (n+1), (1 << n)-1]


def defects(a):
    c = [0, 0, 0]
    for i in range(len(a)):
        for j in range(i+1, len(a)):
            c[((a[i] >> j) & 1) + ((a[j] >> i) & 1)] += 1
    return c


def paths_on_subsets(a):
    n = len(a)
    dp = [[0]*n for _ in range(1 << n)]
    h = [0]*(1 << n)
    h[0] = 1
    for v in range(n):
        dp[1 << v][v] = 1
    for s in range(1, 1 << n):
        for v in range(n):
            count = dp[s][v]
            if count:
                nxt = a[v] & ~s
                while nxt:
                    bit = nxt & -nxt
                    dp[s | bit][bit.bit_length()-1] += count
                    nxt ^= bit
        h[s] = sum(dp[s])
    return h


def covers(h, qmax):
    full = len(h)-1
    f = [0]*len(h)
    f[0] = 1
    out = [f[full]]
    for _ in range(qmax):
        new = [0]*len(h)
        for s in range(len(h)):
            t = s
            while True:
                new[s] += f[s ^ t]*h[t]
                if not t:
                    break
                t = (t-1) & s
        f = new
        out.append(f[full])
    return out


def histogram(a):
    return Counter(sum((a[u] >> v) & 1 for u, v in zip(p, p[1:]))
                   for p in permutations(range(len(a))))


def choose(z, n):
    v = 1
    for j in range(n):
        v *= z-j
    return v // factorial(n)


def evaluate(hist, n, q):
    return sum(c*choose(q+f, n) for f, c in hist.items())


def hc(a):
    n = len(a)
    if n < 3:
        return 0
    return sum(all((a[u] >> v) & 1 for u, v in zip(p, p[1:]+p[:1]))
               for tail in permutations(range(1, n)) for p in [(0,)+tail])


def strong(a):
    n = len(a)
    for root in range(n):
        seen = 1 << root
        while True:
            new = seen
            for u in range(n):
                if seen >> u & 1:
                    new |= a[u]
            if new == seen:
                break
            seen = new
        if seen != (1 << n)-1:
            return False
    return True


def arithmetic(n):
    a = [0]*n
    for x in range(1, n+1):
        for y in (x+2, 2*x):
            if y <= n:
                a[x-1] |= 1 << (y-1)
    return a


def symmetric(a):
    return [row | sum(((a[j] >> i) & 1) << j for j in range(len(a)))
            for i, row in enumerate(a)]


def transitive(n):
    return [sum(1 << j for j in range(i+1, n)) for i in range(n)]


def main():
    digest = hashlib.sha256()
    census = []
    for n in range(5):
        total = 1 << (n*(n-1))
        for mask in range(total):
            a = graph(n, mask)
            h = paths_on_subsets(a)
            values = covers(h, 4)
            hist = histogram(a)
            ca = complement(a)
            cv = covers(paths_on_subsets(ca), 4)
            b = gate(a)
            bh = paths_on_subsets(b)
            bv = covers(bh, 3)
            for q in range(5):
                require(evaluate(hist,n,q) == values[q], ('word-cut',n,mask,q))
                require(evaluate(hist,n,-q) == (-1)**n*cv[q], ('reciprocity',n,mask,q))
            for q in range(1,4):
                require(bv[q] == q*(values[q+1]+q*values[q]-values[q-1]), ('gate',n,mask,q))
            d, e = defects(a), defects(b)
            require((d[0], d[2]) == (e[0], e[2]), ('defect',n,mask))
            if n:
                require(strong(b), ('strong',n,mask))
                require(hc(b) == h[-1], ('cycle-count',n,mask))
            else:
                require(hc(b) == 0 and h[-1] == 1 and bv == [0,1,4,9], 'empty-boundary')
            digest.update(json.dumps([n,mask,values,cv,bv,d],separators=(',',':')).encode()+b'\n')
        census.append({'n':n,'all_labelled_loopless_digraphs':total})
    controls = []
    for n in range(1,11):
        a = arithmetic(n)
        b = symmetric(a)
        t = transitive(n)
        ha,hb,ht = [paths_on_subsets(x) for x in (a,b,t)]
        va,vb,vt = [covers(h,2) for h in (ha,hb,ht)]
        word = list(range(n if n%2 else n-1,0,-2))+list(range(2,n+1,2))
        require(len(word)==n and len(set(word))==n and all((b[x-1]>>(y-1))&1 for x,y in zip(word,word[1:])),('arithmetic-word',n))
        require(ha[-1] == (1 if n<=2 else 0), ('directed-H',n))
        require(n==1 or va[2]==4, ('constant-two-cover',n))
        require(next(q for q,v in enumerate(va) if v) == (1 if n<=2 else 2), ('directed-pc',n))
        require(next(q for q,v in enumerate(vb) if v) == 1, ('symmetric-pc',n))
        require(hb[-1] >= 1 and (n==1 or hb[-1]%2==0), ('symmetric-H',n))
        require(ht[-1] == 1, ('completion-H',n))
        edge_count=max(0,n-2)+n//2-(n>=4)
        require(sum(x.bit_count() for x in a)==edge_count,('edges',n))
        controls.append({'n':n,'directed_H':ha[-1],'symmetric_H':hb[-1],
                         'transitive_H':ht[-1],'P_directed_at_2':va[2],
                         'defect_counts_0_1_2':[defects(x) for x in (a,b,t)],
                         'symmetric_hamilton_word':word})
    hostile=[]
    for name,a in [('incomplete_with_two_cycle',[2,5,0]),('transitive',transitive(3))]:
        hist=histogram(a)
        require([hist[f] for f in range(3)]==[1,4,1],('hostile-hist',name))
        vals=covers(paths_on_subsets(a),5)
        require(vals==[q**3 for q in range(6)],('hostile-P',name))
        b=a[:];response=[];ds=[]
        for r in range(4):
            response.append(paths_on_subsets(b)[-1]);ds.append(defects(b));b=gate(b)
        hostile.append({'name':name,'rows':a,'P_coefficients_ascending':[0,0,0,1],
                        'H_nested_r0_to_3':response,'defects_nested_r0_to_3':ds})
    require(hostile[0]['H_nested_r0_to_3']==hostile[1]['H_nested_r0_to_3'],'observer-blind')
    # Distinguish all order-two comparison types by actual cover evaluations.
    order2=[{'rows':a,'P_at_0_1_2':covers(paths_on_subsets(a),2)} for a in ([0,0],[2,0],[2,1])]
    require([r['P_at_0_1_2'] for r in order2]==[[0,0,2],[0,1,4],[0,2,6]],'minimal-boundary')
    for x in range(1,101):
        require(2*(x+2)==(2*x)+4,'literal-affine')
        h=(x&-x).bit_length()-1
        horizontal=lambda y:y+2*(y&-y)
        require(2*horizontal(x)==horizontal(2*x),'scaled-grid-commutes')
    print(json.dumps({'status':'PASS','scope':'All labelled loopless digraphs n=0..4, arithmetic graph controls n=1..10, explicit order3 hostile; no spectrum completeness claim.',
        'census':census,'exhaustive_rows_sha256':digest.hexdigest(),
        'arithmetic_controls':controls,'same_cover_different_comparison_alphabet':hostile,
        'order2_minimality_controls':order2,'affine_and_row_scaled_controls':'x=1..100'},indent=2,sort_keys=True))

if __name__=='__main__':
    main()
