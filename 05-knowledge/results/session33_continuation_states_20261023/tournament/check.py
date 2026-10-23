#!/usr/bin/env python3
import itertools, math, json, pathlib, hashlib
from functools import lru_cache

def check(ok,msg):
    if not ok: raise RuntimeError(msg)
def tour(n,mask):
    a=[[0]*n for _ in range(n)]
    for b,(i,j) in enumerate(itertools.combinations(range(n),2)):
        a[i][j]=(mask>>b)&1; a[j][i]=1-a[i][j]
    return a

def profile(a):
    n=len(a); f=[0]*(n+1)
    for p in itertools.permutations(range(n)):
        bad=sum(not a[x][y] for x,y in zip(p,p[1:])); good=n-1-bad
        for d in range(bad+1,n+1): f[d]+=math.comb(good,d-1-bad)
    return f

def hp(a):
    n=len(a); dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for s in range(1,1<<n):
        for v in range(n):
            if not dp[s][v]: continue
            for w in range(n):
                if not (s>>w)&1 and a[v][w]: dp[s|1<<w][w]+=dp[s][v]
    return sum(dp[-1])
def substitute(q,fs):
    labels=[(i,v) for i,a in enumerate(fs) for v in range(len(a))]
    return [[(fs[i][u][v] if i==j else q[i][j]) for j,v in labels] for i,u in labels]

def context_coeff(q):
    # Label zero is a module placeholder; other labels are distinct singletons.
    m=len(q)-1; coeff=[0]*(m+2)
    for c in range(1,m+2):
        # Choose placeholder slots, then order exterior vertices independently.
        for slots in itertools.combinations(range(m+c),c):
            if any(y==x+1 for x,y in zip(slots,slots[1:])): continue
            st=set(slots)
            for p in itertools.permutations(range(1,m+1)):
                it=iter(p); word=[0 if j in st else next(it) for j in range(m+c)]
                if all(q[x][y] for x,y in zip(word,word[1:])): coeff[c]+=1
    return coeff

def stirling(n):
    p=[1]+[0]*n
    for i in range(1,n+1): p=[0]+[k*p[k]+p[k-1] for k in range(1,n+1)]
    return [math.factorial(k)*p[k] for k in range(n+1)]
def triangle_row(j):
    g=stirling(j); row=[0]*(j+2)
    for a in range(1,j+2):
        for b in range(1,j+1):
            for c in range(1,j+1):
                w=3 if a==b==c else (1 if max(a,b,c)-min(a,b,c)==1 else 0)
                row[a]+=w*g[b]*g[c]
    return row

def main():
    counts={'context_tournaments':0,'direct_substitution_checks':0,'profile_census':0,'permanent_top_checks':0}
    small=[tour(n,b) for n in range(1,4) for b in range(1<<(n*(n-1)//2))]
    smallf=[profile(a) for a in small]
    for m in range(5):
        for b in range(1<<(m*(m+1)//2)):
            q=tour(m+1,b); coeff=context_coeff(q); counts['context_tournaments']+=1
            check(not any(coeff[m//2+2:]),'run-depth bound')
            outs=[i for i in range(1,m+1) if q[0][i]]
            ins=[i for i in range(1,m+1) if q[i][0]]
            check(not any(coeff[min(len(outs),len(ins))+2:]),'two-side depth bound')
            if m and m%2==0:
                r=m//2; top=0
                if len(outs)==len(ins):
                    permanent=sum(all(q[o][i] for o,i in zip(outs,p)) for p in itertools.permutations(ins))
                    top=math.factorial(r)*permanent
                    check((top==math.factorial(r)**2)==all(q[o][i] for o in outs for i in ins),'extremal observer classification')
                check(coeff[r+1]==top,'top coefficient permanent')
                counts['permanent_top_checks']+=1
            for a,f in zip(small,smallf):
                pred=sum(x*y for x,y in zip(coeff,f))
                check(pred==hp(substitute(q,[a]+[[[0]]]*m)),'independent substitution DP')
                counts['direct_substitution_checks']+=1
    rows=[triangle_row(j) for j in range(1,11)]
    for j,row in enumerate(rows,1): check(row[-1]==math.factorial(j)**2,'triangular diagonal')
    unique={}; pairs={}; fixed_size_observers={}
    for n in range(1,6):
        for b in range(1<<(n*(n-1)//2)):
            f=profile(tour(n,b)); counts['profile_census']+=1
            h=(n-1)//2
            observer=(n,tuple(f[1:h+1]))
            check(observer not in fixed_size_observers or fixed_size_observers[observer]==f,'fixed-size parity observer sufficiency')
            fixed_size_observers[observer]=f
            response=[f[1]]
            for j in range(1,n): response.append(sum(x*y for x,y in zip(triangle_row(j),f)))
            recovered=[0,response[0]]
            for j in range(1,n):
                row=triangle_row(j); rem=response[j]-sum(row[k]*recovered[k] for k in range(1,j+1))
                check(rem%row[j+1]==0,'integrality'); recovered.append(rem//row[j+1])
            check(recovered==f,'profile reconstruction')
            key=(n,f[1]); old=pairs.get(key)
            if old and old[1]!=f and 'same_H' not in unique:
                unique['same_H']={'n':n,'masks':[old[0],b],'profiles':[old[1],f],'triangle_responses':[3*old[1][1]+old[1][2],3*f[1]+f[2]]}
            pairs[key]=(b,f)
    # Reverse orientation invariant holds for every factor, but nonuniform pinned contexts need not.
    q=tour(3,2) # arbitrary explicit triangle checked below
    c3=[[0,1,0],[0,0,1],[1,0,0]]
    a=substitute([[0,1],[0,0]],[[[0]],c3]); ar=[list(r) for r in zip(*a)]
    check(profile(a)==profile(ar),'converse full profile')
    port=None
    for bits in range(16):
        def ext(t):
            r=[row+[int((bits>>i)&1)] for i,row in enumerate(t)]
            r.append([1-int((bits>>i)&1) for i in range(4)]+[0]); return r
        if hp(ext(a))!=hp(ext(ar)):
            port={'factor_a':a,'factor_b':ar,'vertex_to_new_bits':bits,'counts':[hp(ext(a)),hp(ext(ar))]}; break
    check(port is not None,'nonuniform context hostile')
    result={'status':'PASS','counts':counts,'triangle_rows_1_to_10':rows,'same_H_hostile':unique['same_H'],'nonuniform_port_hostile':port,'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__': main()
