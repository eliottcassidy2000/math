#!/usr/bin/env python3
"""Reconstruct endpoint transport in a separate ascending classifier, then transpose."""
import itertools as it
from collections import Counter
from fractions import Fraction
import json, math
import sympy as sp
from pathlib import Path
from independent_classifier import classify
ROOT=Path(__file__).resolve().parent
OUT=Path(__file__).resolve().parent

def fold(scores):
    n=len(scores)
    return min(tuple(sorted(scores)),tuple(sorted(n-1-s for s in scores)))

def main(n):
    free,codes,up_pi,sc,upA=classify(n);m=len(free);N=2**m
    downfree=[(i,j) for i in range(n) for j in range(i-1)]
    to_up=[free.index((j,i)) for i,j in downfree]
    upword=lambda x:sum(((x>>k)&1)<<to_up[k] for k in range(m))
    As=[list(map(list,zip(*upA[upword(x)]))) for x in range(N)]
    ownpi=[up_pi[upword(x)] for x in range(N)]
    # Construct all labels/data independently in first-discovered-node order.
    discovered=list(dict.fromkeys(ownpi));newpi=[discovered.index(u) for u in ownpi]
    reps=[newpi.index(u) for u in range(len(discovered))]
    q=len(discovered)
    d={'n':n,'q':q,'pi':newpi,'reps':reps,
       'D':[newpi.count(u) for u in range(q)],
       'c3':[math.comb(n,3)-sum(math.comb(sum(row),2) for row in As[x]) for x in reps],
       'scores':[sorted(map(sum,As[x])) for x in reps],
       'selfconverse':[sc[upword(x)] for x in reps]}
    candidate=ROOT/f'n{n}.json'
    producer=json.loads(candidate.read_text()) if candidate.exists() else None
    if producer is not None:
        assert all(d[k]==producer[k] for k in d)
    mapping={ownpi[x]:d['pi'][x] for x in range(N)}
    assert len(mapping)==d['q'] and len(set(mapping.values()))==d['q']
    assert all(mapping[ownpi[x]]==d['pi'][x] for x in range(N))
    assert all(sc[upword(x)]==d['selfconverse'][d['pi'][x]] for x in range(N))
    assert all(d['D'][u]==d['pi'].count(u) for u in range(d['q']))
    h=[math.comb(n,3)-sum(math.comb(sum(r),2) for r in a) for a in As]
    scores=[list(map(sum,a)) for a in As]
    rperm=[downfree.index((n-1-j,n-1-i)) for i,j in downfree]
    refl=lambda x:sum(((x>>k)&1)<<rperm[k] for k in range(m))
    q=d['q']; B=[[0]*q for _ in range(q)];K=[[0]*q for _ in range(q)]
    for x in range(N):
        u=d['pi'][x];v=d['pi'][x^(N-1)]
        (B if refl(x)==x else K)[u][v]+=1
        assert h[x]==d['c3'][u]
    if producer is not None:
        assert B==producer['Ablue'] and K==producer['Ablack']
    rowproof=[]
    for u,word in enumerate(d['reps']):
        a=As[word];s=scores[word];mu=1 if d['selfconverse'][u] else 2
        perms=list(it.permutations(range(n)))
        autom=sum(all(a[i][j]==a[p[i]][p[j]] for i in range(n) for j in range(n)) for p in perms)
        hist=[Counter(),Counter()];target_score=Counter();hp=0
        for p in perms:
            if not all(a[x][y] for x,y in zip(p,p[1:])):continue
            hp+=1
            blue=all(a[p[i]][p[j]]==a[p[n-1-j]][p[n-1-i]] for i in range(n) for j in range(n))
            delta=s[p[0]]-s[p[-1]]-1
            hist[int(blue)][delta]+=Fraction(mu,autom)
            ss=[n-1-v for v in s];ss[p[0]]+=1;ss[p[-1]]-=1
            target_score[fold(ss)]+=Fraction(mu,autom)
        for color,M in enumerate([K,B]):
            true=Counter()
            for v in range(q):
                true[d['c3'][v]-d['c3'][u]]+=M[u][v]
            assert +true==+hist[color]
        true_score=Counter()
        for v in range(q):true_score[fold(d['scores'][v])]+=B[u][v]+K[u][v]
        assert +true_score==+target_score
        # Independently obtain SCCs from transitive closure and compare zero current variance.
        reach=[r[:] for r in a]
        for i in range(n):reach[i][i]=1
        for k in range(n):
            for i in range(n):
                for j in range(n):reach[i][j]|=reach[i][k]&reach[k][j]
        components=[];unused=set(range(n))
        while unused:
            v=min(unused);C={j for j in unused if reach[v][j] and reach[j][v]}
            components.append(C);unused-=C
        source=next(C for C in components if all(a[x][y] for x in C for y in range(n) if y not in C))
        sink=next(C for C in components if all(a[y][x] for x in C for y in range(n) if y not in C))
        regular_ends=all(len({sum(a[v][w] for w in C) for v in C})==1 for C in [source,sink])
        zero_variance=len(set(hist[0])|set(hist[1]))==1
        assert zero_variance==regular_ends
        rowproof.append({'zero_current_variance':zero_variance,'source_size':len(source),'sink_size':len(sink),'node':u,'HP':hp,'Aut':autom,'merged_sheets':mu,'blue_delta':dict(hist[1]),'black_delta':dict(hist[0])})
    # Verify the color-conditioned law-of-total-variance with exact rationals.
    BM=sp.Matrix(B);KM=sp.Matrix(K);AM=BM+KM
    bs=[sum(r) for r in B];ks=[sum(r) for r in K];ds=[b+k for b,k in zip(bs,ks)]
    DB=sp.diag(*bs);DK=sp.diag(*ks);DD=sp.diag(*ds)
    invB=sp.diag(*[sp.Rational(1,b) if b else 0 for b in bs])
    invK=sp.diag(*[sp.Rational(1,k) if k else 0 for k in ks])
    LL=DD-AM*DD.inv()*AM
    within=DB-BM*invB*BM+DK-KM*invK*KM
    between=sp.zeros(q)
    for u,(b,k,dv) in enumerate(zip(bs,ks,ds)):
        if not b*k:continue
        v=BM[u,:]/b-KM[u,:]/k
        between+=sp.Rational(b*k,dv)*(v.T*v)
    assert LL==within+between
    result={'zero_current_variance_nodes':[u for u,r in enumerate(rowproof) if r['zero_current_variance']], 'colored_covariance_identity_verified':True,'n':n,'all_tilings':N,'nodes':q,'all_color_row_polynomials_verified':True,'all_score_profile_rows_verified':True,'rows':rowproof}
    if n==5:
        M=[row[:] for row in B]
        for u,v,change in [(0,2,-1),(3,6,-1),(0,3,1),(2,6,1)]:M[u][v]+=change;M[v][u]+=change
        assert all(min(row)>=0 for row in M)
        assert all(sum(a)==sum(b) for a,b in zip(B,M))
        assert all(M[u][u]==B[u][u] for u in range(q))
        # All modified contacts were already allowed; support is NOT preserved.
        assert all(B[u][v]>0 for u,v in [(0,2),(3,6),(0,3),(2,6)])
        def global_joint(A,label):
            ans=Counter()
            for u in range(q):
                for v in range(u,q):ans[tuple(sorted((label[u],label[v])))]+=A[u][v]//(2 if u==v else 1)
            return +ans
        assert global_joint(M,d['c3'])==global_joint(B,d['c3'])
        assert global_joint(M,list(map(fold,d['scores'])))==global_joint(B,list(map(fold,d['scores'])))
        mismatch=[]
        for u in range(q):
            old=Counter();new=Counter()
            for v in range(q):
                old[d['c3'][v]-d['c3'][u]]+=B[u][v]
                new[d['c3'][v]-d['c3'][u]]+=M[u][v]
            if +old!=+new:mismatch.append({'node':u,'before':dict(+old),'after':dict(+new)})
        assert [x['node'] for x in mismatch]==[2,3]
        result['hostile_mismatched_rows']=mismatch
        result['hostile_preserves_support']=False
    return result

if __name__=='__main__':
    data=[main(n) for n in [3,4,5,6]]
    (OUT/'independent_transport.json').write_text(json.dumps(data,indent=2,default=str)+'\n')
    for d in data:print(json.dumps({k:v for k,v in d.items() if k!='rows'},default=str))
