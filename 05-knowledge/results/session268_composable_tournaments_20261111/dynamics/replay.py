#!/usr/bin/env python3
"""Exact controls for odd-regular lexicographic GF(2) dynamics.
No census or Chan--Pak completeness claim. Standard library only.
"""
from collections import defaultdict
from functools import lru_cache
from itertools import product
import json

def check(ok, message):
    if not ok: raise RuntimeError(message)

def c3(): return [2,4,1]
def p7(): return [sum(1<<((i+s)%7) for s in (1,2,4)) for i in range(7)]
def lex(a,b):
    m=len(b)
    return [sum(((1<<m)-1)<<(j*m) for j in range(len(a)) if (a[i]>>j)&1) | (b[k]<<(i*m)) for i in range(len(a)) for k in range(m)]
def mul(a,b):
    ans=[]
    for row in a:
        x=0
        while row:
            bit=row & -row; row-=bit; x ^= b[bit.bit_length()-1]
        ans.append(x)
    return ans

def rank(a):
    piv={}
    for row in a:
        while row:
            k=row.bit_length()-1
            if k in piv: row ^= piv[k]
            else: piv[k]=row; break
    return len(piv)

def profile(a):
    n=len(a); h=[0]*(1<<n)
    dp=[[0]*n for _ in h]
    for i in range(n): dp[1<<i][i]=1
    for s in range(1,1<<n):
        for j in range(n):
            if (s>>j)&1 and s!=(1<<j):
                prev=s^(1<<j)
                dp[s][j]=sum(dp[prev][i] for i in range(n) if (prev>>i)&1 and (a[i]>>j)&1)
        h[s]=sum(dp[s])
    @lru_cache(None)
    def q(s,r):
        if not s: return int(r==0)
        if r==0 or s.bit_count()<r: return 0
        out=0; t=s
        while t:
            out+=h[t]*q(s^t,r-1); t=(t-1)&s
        return out
    return [q((1<<n)-1,r) for r in range(n+1)]

def substituted_h(a,qs):
    bounds=[len(q)-1 for q in qs]
    dp={}
    ans=0
    for r in product(*(range(b+1) for b in bounds)):
        if not any(r): continue
        v=[0]*len(a)
        if sum(r)==1: v[r.index(1)]=1
        else:
            for j in range(len(a)):
                if r[j]:
                    p=list(r);p[j]-=1
                    old=dp[tuple(p)]
                    v[j]=sum(old[i] for i in range(len(a)) if (a[i]>>j)&1)
        dp[r]=v
        if min(r)>0:
            factor=1
            for i,k in enumerate(r): factor*=qs[i][k]
            ans+=sum(v)*factor
    return ans

def tower(word):
    a=[0]
    for symbol in reversed(word): a=lex(c3() if symbol=='3' else p7(),a)
    return a

def invariant(word):
    # Multiplicities of t, t+1, and t^2+t+1.
    pref=1; osc=proj=0
    for symbol in word:
        if symbol=='3': osc+=pref;pref*=3
        else: proj+=3*pref;pref*=7
    return (proj+int(len(word)%2==0),proj+int(len(word)%2==1),osc)

def main():
    rows=[]
    for d in range(1,6):
        a=tower('3'*d);n=len(a);i=[1<<j for j in range(n)];j=[(1<<n)-1]*n
        a3=mul(mul(a,a),a)
        check(a3==(i if d%2 else [x^y for x,y in zip(i,j)]),'C3 tower cube')
        check(rank(a)==n-int(d%2==0),'C3 tower rank')
        rows.append({'depth':d,'vertices':n,'rank':rank(a),'fixed_space_dimension':n-rank([x^y for x,y in zip(a,i)])})
    mixed=[]
    for w in ('3','7','33','37','73','77','337','373','733'):
        a=tower(w);n=len(a);i=[1<<j for j in range(n)]
        a2=mul(a,a);a3=mul(a2,a)
        check(mul(a3,a)==a,'t^4=t')
        m0,m1,m2=invariant(w)
        check(n-rank(a)==m0,'zero multiplicity')
        check(n-rank([x^y for x,y in zip(a,i)])==m1,'one multiplicity')
        check(n-rank([x^y^z for x,y,z in zip(a2,a,i)])==2*m2,'quadratic multiplicity')
        mixed.append({'word':w,'vertices':n,'zero_eigen_dimension':m0,'one_eigen_dimension':m1,'quadratic_blocks':m2})
    qs3,qs7=profile(c3()),profile(p7())
    values={'37':substituted_h(c3(),[qs7]*3),'73':substituted_h(p7(),[qs3]*7)}
    check(values['37']!=values['73'],'H context distinguishes order')
    collisions=[]
    for length in range(1,15):
        seen={}
        for tup in product('37',repeat=length):
            w=''.join(tup); key=(w.count('3'),invariant(w))
            if key in seen:
                collisions.append({'length':length,'first':seen[key],'second':w,'multiplicities':list(key[1])});break
            seen[key]=w
        if collisions: break
    print(json.dumps({'scope':'Exact formula controls, not all tournaments or all spectra','C3_towers':rows,'mixed_towers':mixed,'seed_profiles':{'3':qs3,'7':qs7},'H21':values,'first_binary_word_similarity_collision_through_length14':collisions},indent=2,sort_keys=True))
if __name__=='__main__': main()
