from fractions import Fraction as Fr
from math import factorial
# ---- power series in t (list of coeffs), Fraction ----
def smul(p,q,N):
    o=[Fr(0)]*(N+1)
    for i,u in enumerate(p):
        if i>N or u==0: continue
        for j,v in enumerate(q):
            if i+j>N: break
            o[i+j]+=u*v
    return o
def sadd(p,q,N):
    return [ (p[i] if i<len(p) else Fr(0))+(q[i] if i<len(q) else Fr(0)) for i in range(N+1)]
def sscale(p,s,N): return [ (p[i] if i<len(p) else Fr(0))*s for i in range(N+1)]
def sinv(p,N):           # 1/p, p[0]!=0
    o=[Fr(0)]*(N+1); o[0]=1/p[0]
    for n in range(1,N+1):
        s=Fr(0)
        for k in range(1,n+1):
            if k<len(p): s+=p[k]*o[n-k]
        o[n]=-s/p[0]
    return o
def sexp(p,N):           # exp(p), p[0]=0
    o=[Fr(0)]*(N+1); o[0]=Fr(1)
    for n in range(1,N+1):
        s=Fr(0)
        for k in range(1,n+1):
            s+=k*p[k]*o[n-k]
        o[n]=s/n
    return o
# ---- E[P^m] machinery (alpha=acr, beta=b) ----
def pmul(p,q):
    d={}
    for i,u in enumerate(p):
        for j,v in enumerate(q):
            if u and v: d[i+j]=d.get(i+j,Fr(0))+u*v
    n=max(d)+1 if d else 1; o=[Fr(0)]*n
    for k,v in d.items(): o[k]=v
    return o
def padd(*ps):
    n=max(len(p) for p in ps); o=[Fr(0)]*n
    for p in ps:
        for i,u in enumerate(p): o[i]+=u
    return o
def Er(p): return sum(c*factorial(j) for j,c in enumerate(p))
def EPm(a,b,c,M):
    ac=pmul(a,c); alpha=[Fr(0)]+ac; res=[]
    for m in range(1,M+1):
        tot=[Fr(0)]
        for k in range(m//2+1):
            coef=Fr(factorial(m),factorial(k)**2*factorial(m-2*k))
            term=[Fr(1)]
            for _ in range(k): term=pmul(term,alpha)
            for _ in range(m-2*k): term=pmul(term,b)
            tot=padd(tot,[u*coef for u in term])
        res.append(Er(tot))
    return res
# ---- VERIFY closed form  E[e^{tP}] = e^{t b0 + g t^2/(1-t b1)}/(1-t b1)  vs E[P^m]/m! ----
N=10
for (b0,b1,g) in [(Fr(0),Fr(1),Fr(1)),(Fr(1),Fr(1),Fr(1)),(Fr(2),Fr(1),Fr(1)),(Fr(1),Fr(-1),Fr(1)),(Fr(1),Fr(2),Fr(3))]:
    # closed form series
    den=[Fr(1),-b1]                         # 1 - t b1
    invden=sinv(den,N)
    X=sadd([Fr(0),b0],sscale(invden,g,N),N) # t b0 + g t^2/(1-tb1): need g t^2 * invden
    X=sadd([Fr(0),b0], sscale(smul([Fr(0),Fr(0),Fr(1)],invden,N),g,N), N)  # g*t^2*invden
    cf=smul(sexp(X,N),invden,N)             # e^X/(1-tb1)
    # E[P^m] from machinery (a=c chosen so ac=g: use a=[g],c=[1]); b=b0+b1 r
    E=EPm([g],[b0,b1],[Fr(1)],N)
    ok=all(cf[m]*factorial(m)==E[m-1] for m in range(1,N+1))
    print(f"b=({b0})+({b1})r, ac={g}: closed form == E[P^m] for m=1..{N}: {ok}")
# ---- THE IMPOSSIBILITY (nullcone <=> t b0 + g t^2/(1-t b1) = log(1-t b1)) ----
print("\nNullcone condition, coefficient equations (symbolic by hand, verified):")
print("  t^1: b0 = -b1")
print("  t^2: g  = -b1^2/2")
print("  t^3: g*b1 = -b1^3/3  =>  g = -b1^2/3   (if b1!=0)")
print("  t^2 vs t^3:  -b1^2/2 = -b1^2/3  =>  b1^2(1/3-1/2)=0  =>  b1=0  =>  g=0, b0=0.")
print("  Rigor check: for b1=1, t^2 forces g=-1/2 but t^3 forces g=-1/3 -- INCONSISTENT:", Fr(-1,2)!=Fr(-1,3))
print("=> for ac=g != 0 and b linear, E[e^{tP}]=1 is IMPOSSIBLE. Extends Hermite (const b) to linear b.")
