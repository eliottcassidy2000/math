#!/usr/bin/env python3
"""Independent finite controls for the source-matching note; works under -O."""
from itertools import combinations, permutations

def require(ok, message):
    if not ok: raise RuntimeError(message)

def graph(n, arcs):
    out=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            x,y=(j,i) if (i,j) in arcs else (i,j)
            out[x]|=1<<y
    return out

def two_apex(w,u,plus):
    m=len(w); a=m; b=m+1; out=[0]*(m+2)
    for i in range(m):
        for j in range(i+1,m): out[i]|=1<<j
        if w[i]=='1': out[a]|=1<<i
        else: out[i]|=1<<a
        if u[i]=='1': out[b]|=1<<i
        else: out[i]|=1<<b
    if plus: out[a]|=1<<b
    else: out[b]|=1<<a
    return out

def paths_by_subset(out):
    n=len(out); dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,c in enumerate(dp[mask]):
            if not c: continue
            avail=out[i]&~mask
            while avail:
                bit=avail&-avail; j=bit.bit_length()-1
                dp[mask|bit][j]+=c; avail-=bit
    h=[sum(row) for row in dp]; h[0]=1
    return h

def h_f2(out):
    h=paths_by_subset(out); n=len(out); allv=(1<<n)-1
    return h[-1],sum(h[s]*h[allv^s] for s in range(1<<n))

def odd_cycles(out):
    n=len(out); result=[]
    for k in range(3,n+1,2):
        for comb in combinations(range(n),k):
            for rest in permutations(comb[1:]):
                cycle=(comb[0],)+rest
                if all((out[x]>>y)&1 for x,y in zip(cycle,cycle[1:]+cycle[:1])):
                    result.append((sum(1<<x for x in cycle),k))
    return result

def profile_coeffs_from_cycles(out):
    # h ≤ 9: at most 3 vertex-disjoint directed odd cycles.
    n=len(out); cycles=odd_cycles(out); coeff={n:1}
    def add(start,used,vertices,num):
        for i in range(start,len(cycles)):
            mask,k=cycles[i]
            if mask&used: continue
            power=n-(vertices+k)+(num+1)
            coeff[power]=coeff.get(power,0)+2**(num+1)
            add(i+1,used|mask,vertices+k,num+1)
    add(0,0,0,0)
    return coeff

def ordinal_pad(out,k):
    n=len(out); result=[0]*(k+n)
    for i in range(k):
        for j in range(i+1,k+n): result[i]|=1<<j
    for i in range(n): result[i+k]=out[i]<<k
    return result

def wrap(out):
    n=len(out); a=n; b=n+1; ans=out[:]+[0,0]
    for i in range(n): ans[i]|=1<<a; ans[b]|=1<<i
    ans[a]|=1<<b
    return ans

def eval_poly(coeff,r): return sum(c*r**p for p,c in coeff.items())

low=graph(4,{(0,2)}); high=graph(4,{(0,2),(0,3)})
require(profile_coeffs_from_cycles(low)=={4:1,2:2},'positive low profile')
require(profile_coeffs_from_cycles(high)=={4:1,2:4},'positive high profile')
require(h_f2(low)==(3,24) and h_f2(high)==(5,32),'positive H,F2')

minus=two_apex('0000111','1000010',False)
plus=two_apex('0000111','1000010',True)
expected_minus={9:1,7:16,5:46,3:32,1:2}
expected_plus={9:1,7:14,5:50,3:32,1:2}
require(profile_coeffs_from_cycles(minus)==expected_minus,'negative minus profile')
require(profile_coeffs_from_cycles(plus)==expected_plus,'negative plus profile')
require(h_f2(minus)==(97,4292),'negative minus H,F2')
require(h_f2(plus)==(99,4164),'negative plus H,F2')
for r in range(8):
    require(eval_poly(expected_plus,r)-eval_poly(expected_minus,r)==2*r**5*(2-r*r),'full difference')
# Literal graph counting controls include the wrappers, not just the scalar formulas.
for k in range(2):
    for (x,h,f) in ((low,3,24),(high,5,32),(minus,97,4292),(plus,99,4164)):
        obs=paths_by_subset(wrap(ordinal_pad(x,k)))[-1]
        require(obs==h+(2**k)*f,'wrapper count')
require(paths_by_subset(wrap(ordinal_pad(low,4)))[-1]==387,'positive source +130')
require(paths_by_subset(wrap(ordinal_pad(high,4)))[-1]==517,'positive target +130')
print('PASS: exact odd-cycle profiles, independent H and F2 subset DP, and literal wrapped graphs')
# Incoming P=6,R=8 pair, checked by a separate literal graph constructor.
def balanced(P,R,plus):
    p=bin(P)[2:]; s=''.join('1' if c=='0' else '0' for c in bin(R)[2:])[::-1]
    return two_apex(p+'01'+s,p+'10'+s,plus)
oldminus=balanced(6,8,False); oldplus=balanced(6,8,True)
require(h_f2(oldminus)==(2625,131004),'incoming negative minus')
require(h_f2(oldplus)==(2627,130876),'incoming negative plus')
require(paths_by_subset(wrap(oldminus))[-1]==133629,'incoming negative wrapped minus')
require(paths_by_subset(wrap(oldplus))[-1]==133503,'incoming negative wrapped plus')
print('PASS: incoming balanced P=6,R=8 literal source and wrapper')
