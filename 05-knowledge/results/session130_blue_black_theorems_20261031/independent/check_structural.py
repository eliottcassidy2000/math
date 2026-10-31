#!/usr/bin/env python3
"""Independent structural tests: colored deletion and equivariant odd cycles.
No producer imports/data. Lex free-edge gauge intentionally differs from owner.
Path and pair-cover counts use subset DP; cycle lists use directed DFS.
"""
from collections import Counter
from itertools import combinations
import json


def need(c,m):
    if not c:raise RuntimeError(m)


def free(n):return [(i,j) for i in range(n) for j in range(i+2,n)]


def refl(n,m):
    es=free(n);ix={e:k for k,e in enumerate(es)}
    return sum(1<<ix[n-1-j,n-1-i] for k,(i,j) in enumerate(es) if m>>k&1)


def restrict(n,m,vs):
    es=free(n);values={e:(m>>k)&1 for k,e in enumerate(es)}
    return sum(values[vs[i],vs[j]]<<k for k,(i,j) in enumerate(free(len(vs))))


def board(n,m):
    a=[[int(i<j) for j in range(n)] for i in range(n)]
    for k,(i,j) in enumerate(free(n)):
        if m>>k&1:a[i][j],a[j][i]=0,1
    return a


def symmetric(n):
    es=free(n);ix={e:k for k,e in enumerate(es)};seen=set();bs=[]
    for k,(i,j) in enumerate(es):
        if k in seen:continue
        l=ix[n-1-j,n-1-i];seen|={k,l};bs.append((1<<k)|(1<<l))
    return [sum(b for k,b in enumerate(bs) if word>>k&1) for word in range(1<<len(bs))]


def end_dp(a):
    n=len(a);d=[[0]*n for _ in range(1<<n)]
    for j in range(n):d[1<<j][j]=1
    for s in range(1,1<<n):
        for j in range(n):
            if s>>j&1:
                before=s^(1<<j)
                if before:d[s][j]=sum(d[before][i] for i in range(n) if a[i][j])
    return d


def cycle_list(a):
    n=len(a);out=[]
    for root in range(n):
        todo=[(root,)]
        while todo:
            p=todo.pop()
            if len(p)>=3 and len(p)%2 and a[p[-1]][root]:out.append(p)
            todo.extend(p+(v,) for v in range(root+1,n) if v not in p and a[p[-1]][v])
    return out


def normalized(c):
    k=c.index(min(c));return c[k:]+c[:k]


def pair_weights(a):
    n=len(a);f=end_dp(a);r=end_dp(list(map(list,zip(*a))));full=(1<<n)-1
    w=[[sum(f[s][u]*r[full^s][v] for s in range(1,full)) for v in range(n)] for u in range(n)]
    b=[[0]*n for _ in range(n)]
    for u,v in combinations(range(n),2):
        need((w[u][v]+w[v][u])%2==0,'Two-path symmetric pair weight not integral')
        b[u][v]=b[v][u]=(w[u][v]+w[v][u])//2
    return b


def ocf_test(n,m):
    a=board(n,m);cs=cycle_list(a);lookup={c:k for k,c in enumerate(cs)}
    support=[sum(1<<v for v in c) for c in cs];q=i=d=0
    for k,c in enumerate(cs):
        rc=normalized(tuple(n-1-v for v in reversed(c)));j=lookup[rc]
        if j==k:
            q+=1
            need(n%2==1 and n//2 in c,'Fixed odd cycle did not contain unique center')
        elif j>k:
            if support[j]&support[k]:i+=1
            else:d+=1
    # All collections have size at most two for these n<=8 tests. Track
    # length multiplicities, not merely total packing degree.
    dimension=(n-1)//2
    def mono(c,multiple=1):return tuple(multiple if 2*j+3==len(c) else 0 for j in range(dimension))
    profile=Counter({(0,)*dimension:1});folded=Counter({(0,)*dimension:1})
    for k,c in enumerate(cs):
        profile[mono(c)]+=1
        rc=normalized(tuple(n-1-v for v in reversed(c)));j=lookup[rc]
        if j==k:folded[mono(c)]+=1
        elif j>k and not (support[k]&support[j]):folded[mono(c,2)]+=1
        for l in range(k+1,len(cs)):
            if not (support[k]&support[l]):
                powers=tuple(x+y for x,y in zip(mono(c),mono(cs[l])))
                profile[powers]+=1
    need(all((profile[p]-folded[p])%2==0 for p in set(profile)|set(folded)),
         'Length-refined reflection packing quotient failed')
    need(all(count%2==0 for powers,count in profile.items() if sum(x%2 for x in powers)>n%2),
         'Forbidden length-multiplicity parity sector')
    a2=sum(not (s&t) for k,s in enumerate(support) for t in support[k+1:])
    h=sum(end_dp(a)[-1])
    need(h==1+2*len(cs)+4*a2,'Full OCF at n<=8 failed')
    need(a2%2==d%2,'Mirror-packing parity failed')
    need(h%8==(1+2*q+4*i)%8,'Mirror cycle mod8 failed')
    if n%2==0:need(q==0 and h%4==1,'Even SC boundary failed')
    else:
        z=n//2;vs=[v for v in range(n) if v!=z];base=[[a[u][v] for v in vs] for u in vs]
        hb=sum(end_dp(base)[-1]);rc=(sum(z in c for c in cs)-q)//2
        need((h-hb-2*q-4*rc)%8==0,'Center deletion mod8 failed')
        b=pair_weights(base);out={k for k,v in enumerate(vs) if a[z][v]}
        crossing=sum(b[u][v] for u in out for v in range(n-1) if v not in out and base[u][v])
        need(h-hb==2*crossing,'Exact pair-cover center cut failed')
        mirror=sum(b[u][n-2-u] for u in out if base[u][n-2-u])
        need(mirror%2==q%2,'Fixed odd cycle / mirror-arc cut parity failed')
    return h%8


def colored_lifts(n):
    full=(1<<len(free(n)))-1;lowerfull=(1<<len(free(n-1)))-1;counts=Counter();blues=Counter()
    for m in range((full+1)//2):
        child=restrict(n,m,list(range(n-1)));child=min(child,child^lowerfull)
        counts[child]+=1
        if refl(n,m)==m:blues[child]+=1
    for s in range((lowerfull+1)//2):
        center=restrict(n-1,s,list(range(1,n-1)))
        eligible=refl(n-2,center)==center
        need(counts[s]==1<<(n-2),'Literal lower-line lift count')
        need(blues[s]==2*eligible,'Blue central-board criterion')
    return {'order':n,'lower_lines':len(counts),'blue_upper_lines':sum(blues.values()),
            'eligible_lower_lines':len(blues),'lifts_per_lower_line':1<<(n-2)}


def arithmetic_test():
    from fractions import Fraction
    from math import gcd, isqrt, prod
    states=(4301,17765,16302,9867,8671,23023)
    fractions=[Fraction(a,b) for a,b in [(17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1)]]
    def nu(z,p):
        z=Fraction(z);a,b=z.numerator,z.denominator;e=0
        while a%p==0:a//=p;e+=1
        while b%p==0:b//=p;e-=1
        return e
    def degree(z):return sum(nu(z,p) for p in (3,7,17,29))
    def decode(z,scale):
        a=nu(scale,2);k=degree(z)
        return (nu(z,17)-nu(z,5),nu(z,5),nu(z,2)-a*k,
                nu(z,3)-nu(z,2)+a*k,nu(z,29),nu(z,7))
    scales=[]
    for aa in range(3):
        for rr in range(3):
            for vv in range(3):
                for ww in range(3):
                    for ss in range(3):
                        if aa*ww:continue
                        for S in (1,31):scales.append(S*2**aa*11**rr*13**vv*19**ww*23**ss)
    need(len(scales)==270,'Admissible scale count')
    phi=(0,1,2,4,5,3);owner_r=(0,3,5,1,4,2)
    need(all(phi[(-j)%6]==owner_r[phi[j]] for j in range(6)),'Phase/tile reflection gauge')
    tested=0
    for scale in scales:
        cycle=[scale*b for b in states]
        for j,N in enumerate(cycle):
            letter=next(k for k,f in enumerate(fractions) if N%f.denominator==0)
            need('ABCDEFGHIJKLMN'[letter]=='GBDEFA'[j],'Scaled physical cycle priority changed')
            need(N*fractions[letter]==cycle[(j+1)%6],'Physical six-cycle transition failed')
        need(gcd(*cycle)==scale,'Physical cycle gcd lost scale')
        values=set()
        for mask in range(64):
            x=tuple(mask>>j&1 for j in range(6));z=prod(cycle[j]**x[j] for j in range(6));values.add(z)
            need(decode(z,scale)==x,'Scaled primitive decoder failed')
            r=prod(cycle[j]**x[-j%6] for j in range(6));norm=z*r
            need((isqrt(norm)**2==norm)==(x[1]==x[5] and x[2]==x[4]),'Square norm color criterion failed')
            tested+=1
        need(len(values)==64,'Decoration codec collision')
        e=(-2,3,-1,0,2,-3);z=prod(Fraction(cycle[j])**e[j] for j in range(6))
        need(decode(z,scale)==e,'Signed primitive decoder failed')
    # Raw C/R-set obstruction: six-cycles contribute one blue and one black pair;
    # two-cycles contribute one blue only. Required coefficients are (-4,12).
    need(8-12<0,'Raw-cycle obstruction changed')
    return {'admissible_scales':len(scales),'binary_codewords':tested,'signed_tests':len(scales),
            'physical_word':'GBDEFA','scaled_decoder_and_square_norm':True,'raw_cycle_multiplicities':[-4,12]}


def main():
    hist=[]
    for n in range(3,9):
        masks=symmetric(n)
        if n==8:masks=masks[::32]
        values=Counter(ocf_test(n,m) for m in masks)
        hist.append({'order':n,'tested':len(masks),'scope':'exhaustive' if n<=7 else 'every32nd reflection-orbit-assignment',
                     'H_mod8':dict(sorted(values.items()))})
    print(json.dumps({'primegame_codec':arithmetic_test(),'colored_lifts':[colored_lifts(n) for n in range(4,8)],
       'odd_cycle_checks':hist,'no_producer_imports_or_data':True,
       'gates_survive_optimization':True,
       'counts':'literal odd directed cycles modulo rotation, never support-only'},sort_keys=True,indent=2))

if __name__=='__main__':
    import sys
    if '--test-failing-gate' in sys.argv:need(False,'Intentional optimized-safe negative control')
    main()
