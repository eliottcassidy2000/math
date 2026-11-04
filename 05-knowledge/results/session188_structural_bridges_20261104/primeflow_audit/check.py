"""Independent bounded audit of weighted prime-factor closure, without author imports."""
import itertools,json,math,pathlib,hashlib

def require(ok,detail):
    if not ok: raise RuntimeError(detail)

def prime(n):
    return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))

def ordp(p,r):
    return next(d for d in range(1,r) if pow(p,d,r)==1)

def val(n,r):
    require(n>0,'nonpositive valuation input'); k=0
    while n%r==0: n//=r; k+=1
    return k

def profiles(p,targets,N):
    out=[]; s=0
    for n in range(1,N+1):
        s=p*s+1
        out.append(tuple(val(s,r) for r in targets))
    return out

def rank(p,r,k):
    if not k: return 1
    d=ordp(p,r)
    return r**k if d==1 else d*r**max(0,k-val(p**d-1,r))

def formula(p,R,n):
    ans=[]
    for r in R:
        d=ordp(p,r)
        ans.append(val(n,r) if d==1 else (val(p**d-1,r)+val(n,r) if n%d==0 else 0))
    return tuple(ans)

def predict(p,R,euler,u,v):
    ranks=[rank(p,r,k) for r,k in zip(R,u)]
    if any(q%(4 if euler else 2)==0 for q in ranks): return None
    L=math.lcm(2 if euler else 1,*ranks)
    if any(x>cap for x,cap in zip(formula(p,R,L),v)): return None
    if L>1: return L
    for ell in range(3,1000,2):
        if prime(ell) and all(x<=cap for x,cap in zip(formula(p,R,ell),v)): return ell
    raise RuntimeError('prime witness exceeded audit bound')

primes=[3,5,7,11,13,17,19,23,29,31,37,41,43]
counts={'rank_equivalences':0,'source_role_target_pairs':0,'boxes':0,'least_completions':0,'fresh_extensions':0}
for p in primes:
    Rpool=[r for r in primes[:9] if r!=p]
    for R in itertools.combinations(Rpool,2):
        direct=profiles(p,R,240)
        for n,values in enumerate(direct,1):
            require(values==formula(p,R,n),('valuation law',p,R,n))
            for i,r in enumerate(R):
                for k in range(1,4):
                    require((values[i]>=k)==(n%rank(p,r,k)==0),('rank',p,r,n,k))
                    counts['rank_equivalences']+=1
        for euler in (False,True):
            if euler and p%4!=1:continue
            physical=[(n,row) for n,row in enumerate(direct,1) if (n%4==2 if euler else n>=3 and n%2==1)]
            counts['source_role_target_pairs']+=1
            for specs in itertools.product([(a,b) for a in range(3) for b in range(a,3)],repeat=2):
                u=tuple(s[0] for s in specs); v=tuple(s[1] for s in specs)
                predicted=predict(p,R,euler,u,v)
                actual=next((n for n,row in physical if all(a<=x<=b for a,x,b in zip(u,row,v))),None)
                require(actual==(predicted if predicted is not None and predicted<=240 else None),('box',p,R,euler,u,v,predicted,actual))
                counts['boxes']+=1
            for u in itertools.product(range(3),repeat=2):
                ranks=[rank(p,r,k) for r,k in zip(R,u)]
                if any(q%(4 if euler else 2)==0 for q in ranks):continue
                L=math.lcm(2 if euler else 1,*ranks); C=formula(p,R,L)
                require(all(a<=x for a,x in zip(u,C)),('extensive',p,R,u))
                for n,row in physical:
                    if all(a<=x for a,x in zip(u,row)):
                        require(n%L==0 and all(a<=x for a,x in zip(C,row)),('least completion',p,R,u,n))
                counts['least_completions']+=1
                for ell in (47,53):
                    require(formula(p,R,L*ell)==C,('fresh prime',p,R,L,ell))
                    counts['fresh_extensions']+=1

R=(4561,368089,71)
require(all(prime(r) for r in R),'M3 primality')
require(tuple(ordp(3,r) for r in R)==(15,21,35),'M3 orders')
M={n:tuple(val((3**n-1)//2,r) for r in R) for n in (3,15,21,35,105)}
require(M=={3:(0,0,0),15:(1,0,0),21:(0,1,0),35:(0,0,1),105:(1,1,1)},'M3 exact profiles')
require(formula(3,(11,23),55)==(3,1),'weighted join')
require(formula(5,(3,19),18)==(3,1),'Euler overflow')
require(formula(3,(13,),3)==(1,) and formula(13,(3,),18)==(2,) and formula(13,(3,),2)==(0,),'local/global hostile')
root=pathlib.Path(__file__).resolve().parent.parent
source={f:hashlib.sha256((root/'primeflow'/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py']}
print(json.dumps({'status':'PASS','universe':'sources 3..43, two-target subsets of odd primes <=29 excluding source; all admissible roles; direct lengths <=240; bounds 0..2','counts':counts,'M3':M,'source_sha256':source},indent=2))
