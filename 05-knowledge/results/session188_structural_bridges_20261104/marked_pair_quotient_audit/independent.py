#!/usr/bin/env python3
"""Independent grouped-quartic audit, using literal word cuts and predecessor HP DP."""
from itertools import combinations,permutations,product
from functools import lru_cache
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent.parent

def need(ok,why):
    if not ok:raise RuntimeError(why)

def rows(n,code):
    a=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if code>>k&1:a[i]|=1<<j
        else:a[j]|=1<<i
    return tuple(a)

def opposite(a):return tuple(((1<<len(a))-1)^r^(1<<i) for i,r in enumerate(a))
def induced(a,names):return tuple(sum(1<<j for j,v in enumerate(names) if a[u]>>v&1) for u in names)
def codeof(a):return sum(1<<k for k,(i,j) in enumerate(combinations(range(len(a)),2)) if a[i]>>j&1)

def blow(Q,blocks):
    groups=[];off=0
    for b in blocks:groups.append(list(range(off,off+len(b))));off+=len(b)
    a=[0]*off
    for i,b in enumerate(blocks):
        for u,row in enumerate(b):a[groups[i][u]]=sum(1<<groups[i][v] for v in range(len(b)) if row>>v&1)
        for j in range(len(Q)):
            if Q[i]>>j&1:
                for u in groups[i]:a[u]|=sum(1<<v for v in groups[j])
    return tuple(a),groups

def extend(a,x,c):
    n=len(a);r=list(a)+[1<<(n+1),0]
    for i in range(n):
        if x>>i&1:r[i]|=1<<n
        else:r[n]|=1<<i
        if c>>i&1:r[i]|=1<<(n+1)
        else:r[n+1]|=1<<i
    return tuple(r)

@lru_cache(None)
def hp(a):
    n=len(a);d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,1<<n):
        if s.bit_count()==1:continue
        for j in range(n):
            if s>>j&1:
                t=s^(1<<j);d[s][j]=sum(d[t][i] for i in range(n) if t>>i&1 and a[i]>>j&1)
    return sum(d[-1]) if n else 1

@lru_cache(None)
def covers_by_words(a):
    n=len(a);C=[[0]*n for _ in range(n)]
    for w in permutations(range(n)):
        bad=[k for k in range(1,n) if not a[w[k-1]]>>w[k]&1]
        if len(bad)>1:continue
        for k in bad or range(1,n):C[w[k-1]][w[k]]+=1
    return tuple(tuple(r) for r in C)

def grouped(T,U,p,q,A,B,vector=None):
    need(A and B and not set(A)&set(B),'disjoint nonempty groups')
    maskA=sum(1<<i for i in A);maskB=sum(1<<i for i in B);value=0;calls=0
    for bits in range(16):
        x=((1<<q) if bits&1 else 0)|(maskA if bits&2 else 0);c=((1<<p) if bits&4 else 0)|(maskB if bits&8 else 0)
        ft,fu=hp(extend(T,x,c)),hp(extend(U,x,c));calls+=2;value+=(-1 if bits.bit_count()%2 else 1)*(ft-fu)
    keep=[i for i in range(len(T)) if i not in (p,q)];index={v:i for i,v in enumerate(keep)};R=induced(T,keep);C=covers_by_words(R)
    want=sum(C[index[j]][index[i]]-C[index[i]][index[j]] for i in A for j in B)
    need(value==want,'literal grouped coefficient and oriented cover skew')
    if vector is not None:
        physical=sum(vector[((1<<q)|(1<<i))|(((1<<p)|(1<<j))<<len(T))] for i in A for j in B)
        need(physical==value,'aliasing sums exactly the physical quartics')
    return value,calls

def main():
    exhaustive=0;calls=0;physical_quartics=0
    for k in (3,4):
        n=k+1;bank=json.loads((ROOT/'simultaneous_response'/f'full_response_n{n}.json').read_text());basis=[a|(c<<n) for a,c in bank['basis']];vectors=[dict(zip(basis,r['vector'])) for r in bank['records']]
        for source in range(1<<(k*(k-1)//2)):
            Q=rows(k,source)
            for z in range(k):
                blocks=[(0,) for _ in range(k)];blocks[z]=(2,0);T,groups=blow(Q,blocks);U,_=blow(opposite(Q),blocks);p,q=groups[z]
                W=list(T);W[p]^=1<<q;W[q]^=1<<p;W=tuple(W)
                vt,vu,vw=(vectors[codeof(X)] for X in (T,U,W));d={m:vt[m]-vu[m] for m in basis};derivative={m:vt[m]-vw[m] for m in basis}
                other=[i for i in range(k) if i!=z]
                for colors in product(range(3),repeat=len(other)):
                    I=[other[i] for i,c in enumerate(colors) if c==1];J=[other[i] for i,c in enumerate(colors) if c==2]
                    if not I or not J:continue
                    A=sum((groups[i] for i in I),[]);B=sum((groups[i] for i in J),[])
                    coef,newcalls=grouped(T,U,p,q,A,B,d);calls+=newcalls;exhaustive+=1;physical_quartics+=len(A)*len(B)
                    twinder=sum(derivative[((1<<q)|(1<<i))|(((1<<p)|(1<<j))<<n)] for i in A for j in B)
                    need(twinder==coef,'whole converse quartic sign reduction')
    # Independently replay every retained larger grouped battery, including nontrivial child restoration.
    producer=json.loads((ROOT/'marked_pair_quotient'/'receipt.json').read_text());larger=[]
    for item in producer['general_grouped_coefficients']+producer['strong_prime_controls']:
        Q=tuple(item['quotient']);blocks=list(map(tuple,item['children']));z=item['pair_index'];T,groups=blow(Q,blocks);U,_=blow(opposite(Q),blocks);p,q=groups[z];I,J=item['marked_child_sets'];A=sum((groups[i] for i in I),[]);B=sum((groups[i] for i in J),[])
        need(list(T)==item['rows_T'] and list(U)==item['rows_U'],'physical source and target rows')
        coef,newcalls=grouped(T,U,p,q,A,B);calls+=newcalls;need(coef==item['quartic_coefficient']==item['cover_skew'],'retained coefficient')
        W=list(T);W[p]^=1<<q;W[q]^=1<<p;W=tuple(W);Wop=opposite(W);derivative=0
        for rec in item['literal_battery']:
            x,c=rec['a'],rec['c'];ft,fu,fw,fo=(hp(extend(R,x,c)) for R in (T,U,W,Wop));calls+=4
            need((ft,fu)==(rec['T'],rec['U']),'retained literal counts');need(fu==fo,'uniform child restoration after W converse');derivative+=rec['sign']*(ft-fw)
        need(derivative==coef,'quartic unchanged under whole converse after grouping')
        larger.append(dict(order=len(T),coefficient=coef))
    small=[rows(n,c) for n in range(1,4) for c in range(1<<(n*(n-1)//2))];cutcases=0
    for A,B in product(small,repeat=2):
        T,groups=blow((2,0),(A,B));C=covers_by_words(T);ab=sum(C[i][j] for i in groups[0] for j in groups[1]);ba=sum(C[j][i] for i in groups[0] for j in groups[1]);ha,hb=hp(A),hp(B);qa,qb=sum(map(sum,covers_by_words(A))),sum(map(sum,covers_by_words(B)))
        need((ab,ba)==(ha*hb,(ha+qa)*(hb+qb)),'four ordered cut cases');need((ba>ab)==(len(A)+len(B)>2),'sharp singleton boundary');cutcases+=1
    need(cutcases==producer['ordered_cut_checks'],'cut census')
    proof=ROOT/'marked_pair_quotient'/'PROOF.md';sourcepins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [proof,ROOT/'marked_pair_quotient'/'verify.py',ROOT/'response'/'PROOF.md']}
    need(sourcepins==producer['source_sha256'],'frozen producer sources')
    need((ROOT/'marked_pair_quotient'/'receipt.json').read_bytes()==(ROOT/'marked_pair_quotient'/'receipt_optimized.json').read_bytes(),'producer normal and optimized receipts')
    return dict(status='PASS',exhaustive_small_grouped_batteries=exhaustive,underlying_named_quartics=physical_quartics,literal_integer_count_comparisons=calls,larger_retained_batteries=larger,ordered_cut_controls=cutcases,source_sha256=sourcepins)
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
