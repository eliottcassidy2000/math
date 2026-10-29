"""Exact independent probes for marked joins and quadratic/symplectic shadows.

Standard library only. Explicit gates remain active with python -O.
"""
from fractions import Fraction as F
from itertools import product, combinations
from math import lcm
from pathlib import Path
import json, random

CHECKS = 0
def need(ok, msg):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(msg)

def eye(n): return [[int(i == j) for j in range(n)] for i in range(n)]
def mul(a,b): return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def inverse(a):
    n=len(a); b=[[F(x) for x in row]+list(map(F,e)) for row,e in zip(a,eye(n))]
    for j in range(n):
        k=next(i for i in range(j,n) if b[i][j]); b[j],b[k]=b[k],b[j]
        t=b[j][j]; b[j]=[x/t for x in b[j]]
        for i in range(n):
            if i!=j:
                t=b[i][j]; b[i]=[x-t*y for x,y in zip(b[i],b[j])]
    return [r[n:] for r in b]

def det(a):
    b=[r[:] for r in a]; old=sign=1; n=len(b)
    for k in range(n-1):
        rows=[i for i in range(k,n) if b[i][k]]
        if not rows: return 0
        i=rows[0]
        if i!=k: b[i],b[k]=b[k],b[i]; sign=-sign
        pivot=b[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                v=b[i][j]*pivot-b[i][k]*b[k][j]
                if v % old: raise RuntimeError('Bareiss division')
                b[i][j]=v//old
            b[i][k]=0
        old=pivot
    return sign*b[-1][-1] if n else 1

def pf(a):
    n=len(a)
    if not n: return 1
    return sum((-1)**(j+1)*a[0][j]*pf([[a[s][t] for t in range(1,n) if t!=j]
                 for s in range(1,n) if s!=j]) for j in range(1,n))

def rank(a,p):
    b=[[x%p for x in r] for r in a]; d=0
    for j in range(len(b[0])):
        inds=[i for i in range(d,len(b)) if b[i][j]]
        if not inds: continue
        i=inds[0]; b[i],b[d]=b[d],b[i]
        t=pow(b[d][j],-1,p); b[d]=[x*t%p for x in b[d]]
        for i in range(d+1,len(b)):
            t=b[i][j]; b[i]=[(x-t*y)%p for x,y in zip(b[i],b[d])]
        d+=1
        if d==len(b): break
    return d

def denominator(a): return lcm(*(x.denominator for r in a for x in r))
def port(inv): return [sum(r) for r in inv]
def join(a,r):
    n=len(a)
    return [[a[i%n][j%n] if i//n==j//n else (1 if i//n<j//n else -1)
             for j in range(n*r)] for i in range(n*r)]
def predicted_inverse(ai,r):
    n=len(ai); q=port(ai)
    return [[ai[i%n][j%n] if i//n==j//n else
             (1 if i//n<j//n else -1)*q[i%n]*q[j%n]
             for j in range(n*r)] for i in range(n*r)]

def tournament(n,mask):
    a=[[0]*n for _ in range(n)]
    for b,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j]=1 if mask>>b&1 else -1; a[j][i]=-a[i][j]
    return a

def seed(p,n):
    # A reproducible bounded search, not a census/minimality assertion.
    rng=random.Random(511+n)
    for trial in range(20000):
        mask=rng.getrandbits(n*(n-1)//2); a=tournament(n,mask)
        if abs(pf(a))!=p: continue
        ai=inverse(a); h=lcm(*(x.denominator for x in port(ai)))
        if h==p:
            need(det(a)==p*p and n-rank(a,p)==2, 'seed group')
            return {'p':p,'n':n,'mask':mask,'trial':trial,'matrix':a}
    raise RuntimeError('seed search exhausted')

def quadratic(v,p):
    n=len(v)//2; return sum(v[i]*v[n+i] for i in range(n))%p
def symplectic(v,w,p):
    n=len(v)//2
    return sum(v[i]*w[n+i]-w[i]*v[n+i] for i in range(n))%p
def polar(v,w,p):
    n=len(v)//2
    return sum(v[i]*w[n+i]+w[i]*v[n+i] for i in range(n))%p

def planes(p):
    # Every 2-subspace of F_p^4 appears once by its RREF basis.
    for pivots in combinations(range(4),2):
        free=[(i,j) for i in range(2) for j in range(pivots[i]+1,4) if j not in pivots]
        for vals in product(range(p),repeat=len(free)):
            b=[[0]*4 for _ in range(2)]
            for i,j in enumerate(pivots): b[i][j]=1
            for (i,j),v in zip(free,vals): b[i][j]=v
            yield b

def run():
    seeds=[seed(3,4),seed(5,6),seed(11,8)]
    rows=[]
    for s in seeds:
        a=s['matrix']; ai=inverse(a); p=s['p']
        need(denominator(ai)==p,'base exponent')
        for r in range(1,9):
            k=join(a,r); ki=predicted_inverse(ai,r); n=len(k)
            need(mul(k,ki)==eye(n),'literal inverse product')
            need(det(k)==p**(2*r),'literal determinant')
            need(denominator(ki)==(p if r==1 else p*p),'exponent saturation')
            need(lcm(*(x.denominator for x in port(ki)))==p,'port order fixed')
            nullity=n-rank(k,p)
            need(nullity==r+(r%2),'p rank')
            if r<=3: need(inverse(k)==ki,'independent inverse elimination')
            rows.append({'p':p,'copies':r,'vertices':n,'order':p**(2*r),
                         'exponent':denominator(ki),'p_rank':nullity,
                         'C_p_factors':2*(r%2),'C_p2_factors':r-r%2})
    # Arbitrary skew lattices and ports, including nonprimitive and zero ports.
    # Congruence-free direct models give mixed exponents and h < d.
    model_rows=[]
    for d in [2,3,4,5,9,11,25]:
        a=[[0,d],[-d,0]]; ai=inverse(a)
        for v in [[0,0],[1,0],[2,0],[3,5]]:
            q=[sum(ai[i][j]*v[j] for j in range(2)) for i in range(2)]
            h=lcm(*(x.denominator for x in q))
            for r in [2,3,4]:
                k=[[a[i%2][j%2] if i//2==j//2 else
                    (1 if i//2<j//2 else -1)*v[i%2]*v[j%2]
                    for j in range(2*r)] for i in range(2*r)]
                ki=inverse(k)
                need(denominator(ki)==lcm(d,h*h),'general port exponent')
                vv=v*r
                newq=[sum(ki[i][j]*vv[j] for j in range(2*r)) for i in range(2*r)]
                need(newq==q*r,'marked repeated port')
            model_rows.append({'d':d,'port':v,'h':h,'join_exponent':lcm(d,h*h)})
    # Two-channel attachments escape the one-channel prime-support trap.
    base=[[0,1,1,1],[-1,0,-1,1],[-1,1,0,-1],[-1,-1,1,0]]
    ai=inverse(base); e=[1]*4; w=[-2,0,-2,-2]
    def contract(x,a,y): return sum(x[i]*a[i][j]*y[j] for i in range(len(x)) for j in range(len(y)))
    alpha=contract(e,ai,w); attachments=[]
    for cm,prime in [(6,5),(7,11)]:
        v=[int(cm>>j&1) for j in range(4)]
        b=[[1+w[i]*v[j] for j in range(4)] for i in range(4)]
        k=[base[i]+b[i] for i in range(4)]+[[-b[j][i] for j in range(4)]+base[i] for i in range(4)]
        beta=contract(e,ai,v)
        need(all(abs(x)==1 for row in b for x in row),'genuine tournament cut')
        need(rank(b,101)==2,'attachment rank lower bound')
        need(pf(k)==pf(base)**2*(1-alpha*beta)==prime,'two-channel Pfaffian')
        need(det(k)==prime**2 and denominator(inverse(k))==prime,'new prime group')
        attachments.append({'base_matrix':base,'column_mask':cm,'cut':b,
                            'alpha':str(alpha),'beta':str(beta),'pfaffian':prime,
                            'group_cyclic_factors':[prime,prime]})
    geometry=[]
    for p in [2,3,5,11]:
        total=sp=oq=0
        for u,v in planes(p):
            total+=1
            sp+=symplectic(u,v,p)==0
            oq+=quadratic(u,p)==quadratic(v,p)==polar(u,v,p)==0
        need(total==(p*p+1)*(p*p+p+1),'Grassmannian enumeration')
        need(sp==(p+1)*(p*p+1),'symplectic Lagrangians')
        need(oq==2*(p+1),'quadratic maximal isotropics')
        graphs_sp=graphs_oq=0
        for a,b,c,d in product(range(p),repeat=4):
            u=[1,0,a,c]; v=[0,1,b,d]
            si=symplectic(u,v,p)==0
            qi=quadratic(u,p)==quadratic(v,p)==polar(u,v,p)==0
            need(si==(b==c),'symmetric graph criterion')
            need(qi==(a==d==0 and (b+c)%p==0),'alternating graph criterion')
            graphs_sp+=si; graphs_oq+=qi
        need(graphs_sp==p**3 and graphs_oq==p,'graph count')
        geometry.append({'q':p,'all_planes':total,'symplectic':sp,'quadratic':oq,
                         'symplectic_graphs':graphs_sp,'quadratic_graphs':graphs_oq})
    # Golden period controls: ramification at 5, split roots at 11.
    golden=[]
    f=[[0,1],[1,1]]
    for p in [5,11]:
        power=eye(2)
        for n in range(1,121):
            power=[[x%p for x in row] for row in mul(power,f)]
            if power==eye(2): break
        roots=[x for x in range(p) if (x*x-x-1)%p==0]
        need(n=={5:20,11:10}[p],'golden clock')
        golden.append({'p':p,'roots':roots,'matrix_period':n})
    # Arithmetic calibration: the cited Monsky interface is an input, not
    # re-proved by these small linear-algebra checks.
    arithmetic=[]
    for n in [5,11]:
        d2=int(pow(2,(n-1)//2,n)==n-1)
        dm2=int(pow(n-2,(n-1)//2,n)==n-1)
        monsky=[[d2,d2],[d2,dm2]]
        nul=2-rank(monsky,2)
        kernel=[list(v) for v in product(range(2),repeat=2)
                if all(sum(x*y for x,y in zip(row,v))%2==0 for row in monsky)]
        need(len(kernel)==2**nul,'independent Monsky kernel')
        need(nul=={5:1,11:0}[n],'Monsky 5/11')
        arithmetic.append({'n':n,'Monsky':monsky,'kernel':kernel,'pure_dimension':nul})
    x,y=F(25,4),F(75,8)
    need(y*y==x*x*x-25*x and x.denominator!=1,'E5 nontorsion point input to Nagell-Lutz')
    # Incoming clocked cokernel: hyperbolization is a classical paired repair.
    def digit_sign(n):
        value=1
        while n:
            if n%3==2:value=-value
            n//=3
        return value
    word=[digit_sign(i) for i in range(9)]
    c=[[word[(j-i)%9] for j in range(9)] for i in range(9)]
    hyper=[[0]*9+row for row in c]+[[-c[j][i] for j in range(9)]+[0]*9 for i in range(9)]
    need(det(c)==4864 and det(hyper)==4864**2,'clocked hyperbolization order')
    hi=inverse(hyper)
    need(all(hi[i][j]==-hi[j][i] for i in range(18) for j in range(18)),'hyperbolic alternating pairing')
    clock=[[int(j==(i-1)%9) for j in range(9)] for i in range(9)]
    lift=[clock[i]+[0]*9 for i in range(9)]+[[0]*9+clock[i] for i in range(9)]
    need(mul(list(map(list,zip(*lift))),mul(hyper,lift))==hyper,'paired cyclic clock')
    return {'status':'PASS','checks':CHECKS,'seeds':seeds,'joins':rows,
            'general_ports':model_rows,'two_channel_attachments':attachments,
            'geometry':geometry,'golden':golden,'arithmetic':arithmetic,
            'clock_hyperbolization':{'word':word,'original_order':4864,
                'paired_order':4864**2,'original_odd_19_valuation':1,
                'paired_odd_19_valuation':2}}

if __name__=='__main__':
    print(json.dumps(run(),indent=2,sort_keys=True))
