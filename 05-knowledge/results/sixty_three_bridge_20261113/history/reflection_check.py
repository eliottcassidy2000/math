#!/usr/bin/env python3
"""Independent all-word count formula plus literal rooted reflection controls."""
from itertools import product, combinations
from pathlib import Path
import json
from probe import hp, inv, CORES, flip

def require(x, message):
    if not x: raise AssertionError(message)

def apex(w):
    n=len(w)
    rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]+[0]
    for i,b in enumerate(w):
        if b: rows[n]|=1<<i
        else: rows[i]|=1<<n
    return tuple(rows)

def reflection(w): return tuple(1-b for b in reversed(w))

def weighted_count(w):
    return 1+2*sum(1<<max(j-i-2,0) for i,a in enumerate(w)
                   for j,b in enumerate(w) if i<j and a==1 and b==0)

def birth(rows,root,retain):
    n=len(rows); z=list(rows)
    z[root]|=1<<n
    z.append(((1<<n)-1)^(1<<root))
    return tuple(z),root if retain else n

def fvs(rows, vertices=None):
    n=len(rows)
    if vertices is None: vertices=(1<<n)-1
    triangles=[]
    for a,b,c in combinations(range(n),3):
        if (rows[a]>>b&1 and rows[b]>>c&1 and rows[c]>>a&1) or (rows[b]>>a&1 and rows[c]>>b&1 and rows[a]>>c&1):
            triangles.append((1<<a)|(1<<b)|(1<<c))
    best=0; mask=vertices
    while True:
        if mask.bit_count()>best and all(mask&t!=t for t in triangles): best=mask.bit_count()
        if mask==0: break
        mask=(mask-1)&vertices
    return vertices.bit_count()-best

def main():
    allword=0; rearm=0; family=[]
    for n in range(2,9):
        for middle in product((0,1),repeat=n-2):
            w=(1,)+middle+(0,); z=reflection(w)
            T,U=apex(w),apex(z)
            h,S=inv(T); hh,Q=inv(U)
            require(h==hh==weighted_count(w)==weighted_count(z),'count')
            require(T[0]==U[0] and all((T[v]&1)==(U[v]&1) for v in range(1,n+1)), 'root star')
            require(all((T[i]&((1<<n)-1))==(U[i]&((1<<n)-1)) for i in range(n)), 'carrier')
            require(S[0]==[h,weighted_count(w[1:])-w[1],weighted_count(w[1:]),h+w[1]],'before state')
            d=1-w[-2]
            require(Q[0]==[h,weighted_count(w[:-1])-d,weighted_count(w[:-1]),h+d],'after state')
            allword+=1
            if w[1]==0 and w[-2]==0: rearm+=1
    for b in range(1,5):
        B=1<<b; w=(1,)+(0,)*b+(1,)*b+(0,0); z=reflection(w)
        require([i for i in range(len(w))if w[i]!=z[i]]==[1,b+1,2*b+1],'three changed edges')
        h=3*B*B+4*B-1
        before=inv(apex(w))[1][0]; after=inv(apex(z))[1][0]
        require(before==[h,3*B-1,3*B-1,h], 'family before')
        require(after==[h,B*B+2*B,(B+1)**2,h+1],'family after')
        family.append({'b':b,'word':''.join(map(str,w)),'before':before,'after':after})
    hostile=(1,0,1,0)
    require(reflection(hostile)==hostile and inv(apex(hostile))[1][0]==[9,3,3,9], 'hostile mirror')
    T=CORES[0]; E=((2,4),(3,4),(4,7)); U=flip(T,E)
    require(inv(T)[1][6]==[63,11,11,63], 'literal before')
    require(inv(U)[1][6]==[63,24,25,64], 'literal after')
    N,r=birth(U,6,False); NR,r=birth(N,r,True)
    require(hp(NR)==191,'exact NR child')
    rank_controls=0
    for n in range(1,5):
        pairs=list(combinations(range(n),2))
        for bits in range(1<<len(pairs)):
            T=[0]*n
            for i,(a,b) in enumerate(pairs):
                if bits>>i&1: T[a]|=1<<b
                else: T[b]|=1<<a
            for root in range(n):
                V,y=birth(T,root,False); V,y=birth(V,y,True)
                require(fvs(V)==fvs(T,((1<<n)-1)^(1<<root))+1,'NR feedback law')
                rank_controls+=1
    require(fvs(NR)==2,'literal exit rank')
    out={'status':'PASS','universe':'all endpoint-normalized binary words of lengths2..8; familyb1..4; literal63NR',
         'allword_controls':allword,'defect_zero_to_one_controls':rearm,
         'feedback_rank_controls_all_labelled_orders1_to4':rank_controls,
         'family':family,'hostile':{'word':'1010','state':[9,3,3,9]},
         'literal_NR':{'rows':NR,'root':r,'state':inv(NR)[1][r]}}
    Path(__file__).with_name('reflection_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
