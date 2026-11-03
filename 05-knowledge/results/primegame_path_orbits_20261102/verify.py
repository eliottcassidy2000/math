"""Independent literal-arc replay of the depth-two hostile and arithmetic census.

Standard library only. Explicit checks also run under python -O.
"""
from collections import Counter, deque
from itertools import combinations, permutations
from math import factorial
from fractions import Fraction
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def rows_from_chart(n, x):
    rows = [0] * n
    k = 0
    for a, b in combinations(range(n), 2):
        forward = b == a+1 or bool(x >> k & 1)
        k += b > a+1
        rows[a if forward else b] |= 1 << (b if forward else a)
    return tuple(rows)


def count(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for s in range(1, 1 << n):
        for v, ways in enumerate(dp[s]):
            left = rows[v] & ~s
            while ways and left:
                b = left & -left
                dp[s | b][b.bit_length()-1] += ways
                left -= b
    return sum(dp[-1])


def cmds(n):
    return [('C',)] + [(a,b) for i in range(n) for j in range(i+2,n)
                     for a,b in [(i,j),(j,i)]]


def apply(rows, cmd):
    n = len(rows)
    if cmd == ('C',):
        target = list(rows)
        for a, b in combinations(range(n), 2):
            if b > a+1:
                target[a] ^= 1 << b
                target[b] ^= 1 << a
        return tuple(target)
    a, b = cmd
    p = list(range(n))
    p.insert(b, p.pop(a))
    if not all(rows[u] >> v & 1 for u, v in zip(p, p[1:])):
        return None
    return tuple(sum(((rows[p[i]] >> p[j]) & 1) << j for j in range(n))
                 for i in range(n))


def signature(rows, depth):
    if depth == 0:
        return ()
    return tuple(None if (child := apply(rows,c)) is None else signature(child,depth-1)
                 for c in cmds(len(rows)))


def separator(a, b, max_depth=3):
    queue = deque([(a,b,())])
    seen = {(a,b)}
    while queue:
        x,y,w = queue.popleft()
        if len(w) >= max_depth:
            continue
        for c in cmds(len(x)):
            xx,yy = apply(x,c),apply(y,c)
            if (xx is None) != (yy is None):
                return w+(c,)
            if xx is not None and (xx,yy) not in seen:
                seen.add((xx,yy))
                queue.append((xx,yy,w+(c,)))
    return None


def ordinal(*blocks):
    answer=[]
    offset=0
    total=sum(map(len,blocks))
    for block in blocks:
        later=((1<<total)-1)^((1<<(offset+len(block)))-1)
        answer.extend((row<<offset)|later for row in block)
        offset+=len(block)
    return tuple(answer)


def trace(rows, word):
    answer=[]
    for c in word:
        rows=apply(rows,c)
        answer.append(None if rows is None else list(rows))
        if rows is None:
            break
    return answer


def succeeds(rows, word):
    for c in word:
        rows=apply(rows,c)
        if rows is None:
            return False
    return True


def probe(block):
    a=8*block
    return (('C',),(a+3,a+1),(a+5,a+2))


def inverse(word):
    return tuple(c if c==('C',) else (c[1],c[0]) for c in reversed(word))


def conjunction_word(indices):
    word=()
    for i in indices[:-1]:
        w=probe(i)
        word+=w+inverse(w)
    if indices:
        word+=probe(indices[-1])
    return word


def main():
    a,b=rows_from_chart(8,509),rows_from_chart(8,253)
    require(signature(a,2)==signature(b,2),'depth-two hostile')
    word=separator(a,b)
    require(word is not None and len(word)==3,'sharp depth')
    brute=[]
    adjacency=[]
    for rows in [a,b]:
        paths=[p for p in permutations(range(8))
               if all(rows[u]>>v&1 for u,v in zip(p,p[1:]))]
        brute.append(len(paths))
        adjacency.append(sum(any((u,v) in [(1,5),(5,1)] for u,v in zip(p,p[1:])) for p in paths))
    require(brute==[171,173]==[count(a),count(b)],'two independent counts')
    for r in (a,b):
        reach=list(r)
        for k in range(8):
            for i in range(8):
                if reach[i]>>k&1:
                    reach[i]|=reach[k]
        require(all((s|1<<i)==255 for i,s in enumerate(reach)),'strong core')
        autos=sum(all(((r[i]>>j)&1)==((r[p[i]]>>p[j])&1)
                      for i in range(8) for j in range(8))
                  for p in permutations(range(8)))
        require(autos==1,'asymmetric core')
    # Every guard which mentions the differing long chord uses only vertices 0..6.
    affected=[]
    for c in cmds(8)[1:]:
        p=list(range(8));p.insert(c[1],p.pop(c[0]))
        arcs=list(zip(p,p[1:]))
        if any(set(e)=={1,5} for e in arcs):
            altered=[e for e in arcs if e[1]!=e[0]+1]
            require(all(0<=v<=6 for e in altered for v in e),'guard localization')
            affected.append({'command':c,'new_arcs':altered})
    require(len(affected)==4,'four affected guards')
    triangle=(2,4,1)
    contexts=[]
    # Full signatures checked on arbitrary nontransitive context examples as well.
    for left,right in [((),()),((0,),()),((),(0,)),((0,),(0,)),
                       (triangle,()),((),triangle),(triangle,(0,))]:
        aa,bb=ordinal(left,a,right),ordinal(left,b,right)
        require(signature(aa,2)==signature(bb,2),'context depth-two identity')
        ha,hb=count(aa),count(bb)
        scale=(count(left) if left else 1)*(count(right) if right else 1)
        require([ha,hb]==[171*scale,173*scale],'context count')
        contexts.append({'left':list(left),'right':list(right),'counts':[ha,hb]})
    cubes=[]
    for r in range(1,5):
        reference=None
        all_h=[]
        for bits in range(1<<r):
            rows=ordinal(*(a if bits>>i&1 else b for i in range(r)))
            sig=signature(rows,2)
            if reference is None:
                reference=sig
            require(reference==sig,'whole word-cube depth-two equality')
            h=171**bits.bit_count()*173**(r-bits.bit_count())
            all_h.append(h)
            for i in range(r):
                if bits>>i&1:
                    reset=probe(i)+inverse(probe(i))
                    require(trace(rows,reset)[-1]==list(rows),'literal source reset')
            for subset in range(1<<r):
                indices=[i for i in range(r) if subset>>i&1]
                w=conjunction_word(indices)
                require(succeeds(rows,w)==(bits&subset==subset),'physical conjunction')
                require(len(w)==(6*len(indices)-3 if indices else 0),'word length')
        means=[sum(Fraction(all_h[z],1<<(r-1)) for z in range(1<<r)
                   if z.bit_count()%2==parity) for parity in [0,1]]
        require(means==[172**r+1,172**r-1],'parity-law means')
        for subset in range((1<<r)-1):
            moments=[sum(1 for z in range(1<<r) if z.bit_count()%2==p
                         and z&subset==subset) for p in [0,1]]
            require(moments[0]==moments[1],'proper conjunction moments')
        # No claim these means themselves are Hamilton counts of the source cube.
        cubes.append({'blocks':r,'states':1<<r,'conjunction_checks':4**r,
                      'even_mean_H':int(means[0]),'odd_mean_H':int(means[1]),
                      'number_of_H_values':len(set(all_h))})
    fractions=[(17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),
               (77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1)]
    nums={u for u,v in fractions};dens={v for u,v in fractions}
    require((len(nums),len(dens),len(nums|dens),len(nums&dens))==(12,14,19,7),'incidence counts')
    # Bare denominator->numerator incidence differs from actual divisibility dynamics.
    parent={v:v for v in nums|dens}
    def root(v):
        while parent[v]!=v:
            v=parent[v]
        return v
    for u,v in fractions:
        parent[root(v)]=root(u)
    components=Counter(root(v) for v in nums|dens)
    primes=[2,3,5,7,11,13,17,19,23,29]
    weights=[-4,0,-4,2,2,2,2,3,-1,3]
    block=[1,3,2,3,2,2,3,2,1,0,0,3,1,1]
    def valuation(x,p):
        e=0
        while x%p==0:
            x//=p;e+=1
        return e
    delta=[sum(w*(valuation(u,p)-valuation(v,p)) for p,w in zip(primes,weights)) for u,v in fractions]
    require(delta==[[0,-2,1,0][k] for k in block],'inherited block character')
    cycle=[4301]
    fired=[]
    for _ in range(6):
        for k,(u,v) in enumerate(fractions):
            if cycle[-1]%v==0:
                cycle.append(cycle[-1]//v*u);fired.append(chr(65+k));break
    require(cycle==[4301,17765,16302,9867,8671,23023,4301]
            and ''.join(fired)=='GBDEFA','inherited priority six-cycle')
    # At five vertices reconstruct blue-class split from literal path charts.
    # Classification uses minimum full adjacency mask over all vertex permutations.
    def canon(r):
        n=len(r)
        return min(sum(((r[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate(combinations(range(n),2)))
                   for p in permutations(range(n)))
    census={}
    for x in range(64):
        r=rows_from_chart(5,x);c=canon(r)
        reflected=tuple(sum(((r[4-j]>>(4-i))&1)<<j for j in range(5)) for i in range(5))
        row=census.setdefault(c,{'H':count(r),'blue':0,'black':0,'converse':canon(reflected)})
        row['blue' if reflected==r else 'black']+=1
    sc=[v for c,v in census.items() if c==v['converse']]
    pure=[v for v in sc if not v['black']]
    mixed=[v for v in sc if v['black']]
    require((len(census),len(sc),len(pure),len(mixed))==(12,8,3,5),'five-class split')
    result={
        'status':'FINITE-EXACT literal replay; all-size context extension proved in THEORY.md',
        'hostile':{'chart_masks':[509,253],'rows':[a,b],'H':brute,
                   'paths_using_changed_edge':adjacency,'shortest_legality_separator':word,
                   'traces':[trace(a,word),trace(b,word)],'affected_guards':affected},
        'contexts':contexts,
        'word_cubes':cubes,
        'primegame':{'numerators':sorted(nums),'denominators':sorted(dens),
                     'shared':sorted(nums&dens),'numerator_only':sorted(nums-dens),
                     'denominator_only':sorted(dens-nums),'union':sorted(nums|dens),
                     'incidence_component_sizes':sorted(components.values()),
                     'instruction_cycle_rank':len(fractions)-len(nums|dens)+len(components),
                     'height_increments':delta,'priority_cycle':cycle,'priority_word':fired},
        'traps_inherited_orbit_arithmetic':{'480':2*factorial(6)//3,'288':2*factorial(6)//5},
        'five_classes':{'ordinary':12,'self_converse':8,'merged':10,
                        'pure_blue_H':sorted(v['H'] for v in pure),
                        'mixed_self_converse_H':sorted(v['H'] for v in mixed),
                        'details':census},
    }
    (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ['five_classes','primegame','contexts']},indent=2))


if __name__=='__main__':
    main()
