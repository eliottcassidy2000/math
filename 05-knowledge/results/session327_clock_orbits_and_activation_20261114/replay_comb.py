"""S327 exact comb/star geometry certificate. Copied elementary helpers credited to S325/S326."""
"""Literal certificates for 99*2^m + [13,15,...,35] and the independent overlaps."""
from itertools import product

def need(q):
    if not q:
        raise RuntimeError('incorrect corridor certificate')

def tt(n): return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
def ordinal(A,B):return tuple(x|(((1<<len(B))-1)<<len(A)) for x in A)+tuple(x<<len(A) for x in B)
def add(A,mask): return tuple(x|((1<<len(A)) if not mask>>i&1 else 0) for i,x in enumerate(A))+(mask,)
def star(I,U,W,ic,wc):return add(ordinal(ordinal(I,U),W),ic|(wc<<(len(I)+len(U))))
def reach(T,i):
    seen={i}; todo=[i]
    while todo:
        a=todo.pop()
        for j in range(len(T)):
            if T[a]>>j&1 and j not in seen:
                seen.add(j); todo.append(j)
    return seen
def components(T):
    R=[reach(T,i) for i in range(len(T))]
    return {tuple(j for j in sorted(ri) if i in R[j]) for i,ri in enumerate(R)}
def induce(T,vs):return tuple(sum(1<<k for k,j in enumerate(vs) if T[i]>>j&1) for i in vs)
def validate(T):
    need(all(not (r>>i&1) and r>>len(T)==0 for i,r in enumerate(T)))
    need(all((T[i]>>j&1)+(T[j]>>i&1)==1 for i in range(len(T)) for j in range(i+1,len(T))))
def endtable(T):
    n=len(T); M=(1<<n)-1; E=[[0]*n for _ in range(M+1)]
    for i in range(n): E[1<<i][i]=1
    for S in range(1,M+1):
        for i,a in enumerate(E[S]):
            if a:
                r=T[i]&~S
                while r:
                    b=r&-r; r-=b; E[S|b][b.bit_length()-1]+=a
    return E
def profile(T,mask):
    M=(1<<len(T))-1
    e=endtable(T); p=endtable(tuple(M^r^(1<<i) for i,r in enumerate(T)))
    h=[1]+[sum(r) for r in e[1:]]
    W=[sum(k for i,k in enumerate(r) if mask>>i&1) for r in p]
    F=[sum(k for i,k in enumerate(r) if mask>>i&1) for r in e]
    E=[sum(k for i,k in enumerate(r) if not(mask>>i&1)) for r in e]
    return (h[M],
            sum(W[s]*h[M^s] for s in range(1,M+1)),
            sum(F[s]*W[M^s] for s in range(1,M)),
            h[M]+sum(E[s]*h[M^s] for s in range(1,M+1)),
            sum(endtable(add(T,mask))[-1]))
def count(T): return sum(endtable(T)[-1])

from itertools import permutations
from collections import defaultdict,Counter
from fractions import Fraction
from pathlib import Path
import json
logical=need
def extend(A,cut):
    n=len(A)
    return tuple(r|(((cut>>i&1)^1)<<n) for i,r in enumerate(A))+(cut,)

def canon(A):
    groups=defaultdict(list)
    for i,r in enumerate(A): groups[r.bit_count()].append(i)
    best=None
    for runs in product(*(tuple(permutations(groups[d])) for d in sorted(groups))):
        p=tuple(x for run in runs for x in run)
        row=tuple(sum((A[i]>>j&1)<<k for k,j in enumerate(p)) for i in p)
        if best is None or row<best: best=row
    return best

def strong(A):
    for i in range(len(A)):
        rb=1<<i
        while True:
            b=rb
            for j,r in enumerate(A):
                if rb>>j&1: b|=r
            if b==rb: break
            rb=b
        if rb!=(1<<len(A))-1: return False
    return True

def hp_all(A):
    n=len(A); N=1<<n
    E=[[0]*n for _ in range(N)]; H=[1]+[0]*(N-1)
    for i in range(n): E[1<<i][i]=1
    for s in range(1,N):
        H[s]=sum(E[s])
        for i,e in enumerate(E[s]):
            if not e: continue
            W=A[i]&~s
            while W:
                b=W&-W; W-=b
                E[s|b][b.bit_length()-1]+=e
    return H

def pair_profile(A):
    H=hp_all(A)
    return H[-1],sum(a*b for a,b in zip(H,H[::-1]))

def catalog(qmax):
    logical(qmax<=319) # orders 7 and greater are excluded by triangle lower bounds
    T={(0,)}; levels=[1]; prof={(1,2):(0,)}
    for n in range(2,7):
        T={canon(extend(A,c)) for A in T for c in range(1<<(n-1))}
        levels.append(len(T))
        for A in T:
            if strong(A):
                h,q=pair_profile(A)
                if q<=qmax: prof.setdefault((h,q),A)
    logical(levels==[1,1,2,4,12,56])
    reached={(1,1):()}
    for b in range(1,qmax+1):
        for (h,q),L in tuple(reached.items()):
            if q!=b: continue
            for (a,c),X in prof.items():
                if q*c<=qmax: reached.setdefault((h*a,q*c),L+(X,))
    return reached


def main():
    P=catalog(304);need(len(P)==36)
    D={}
    for h,q in P:
        if q>1:D.setdefault(q,set()).add(h)
    expected={2:{1},4:{1},8:{1},12:{3},16:{1},24:{3},32:{1,5},48:{3},64:{1,5},84:{9},96:{3},100:{11},104:{13},108:{15},120:{15},128:{1,5},144:{9},168:{9},192:{3},200:{11},208:{13},216:{15},224:{17},240:{15},256:{1,5,19},272:{23},288:{9,27},296:{23},304:{25}}
    need(D==expected)
    B=json.loads(Path(__file__).with_name('star_bank.json').read_text());phase={int(c):z['shift']for c,z in B.items()}
    need(set(phase)==set(range(1,198,2)))
    cert=[];tight=0
    for c,z in B.items():
        c=int(c);I=tuple(z['I']);W=tuple(z['W']);ic=z['ic'];wc=z['wc']
        h,alpha,lam,_,_=profile(I,ic);k,_,_,beta,d=profile(W,wc)
        s=ic.bit_count();y=len(W)-wc.bit_count();A=alpha*beta
        need(s>=1 and alpha>=s*h and beta>=d and lam<=alpha and beta>=(y+1)*k)
        need(A==99*2**z['shift'] and h*d-k*lam==c)
        need(A-s*c==beta*(alpha-s*h)+s*h*(beta-d)+s*k*lam>=0)
        need(A+(y+1)*c==alpha*(beta-(y+1)*k)+(y+1)*k*(alpha-lam)+(y+1)*h*d>0)
        r=0
        while 99*2**r<s*c:r+=1
        tight+=r==z['shift']
        cert.append((c,s,y,A,c))
    need(tight==67)
    lo=29327;hi=30119;raw={};actual={}
    for N in range(lo,hi+1,2):
        for h,q in P:
            if q<=1 or (N-99*q)%h:continue
            c=(N-99*q)//h
            if c not in range(1,198,2):continue
            raw.setdefault(N,[]).append((h,q,c));r=phase[c]
            if q%2**r==0 and (h,q//2**r) in P and q//2**r>1:
                actual.setdefault(N,[]).append((h,q,c,r))
    comb={N for N in range(lo,hi+1,2)if N%18==9 or N%46 in (25,41)}
    upper={N for N in range(lo,hi+1,2)if N%18==9 or N%46==25}
    need(set(raw)==comb and len(raw)==75 and len(upper)==60 and len(actual)==43)
    need({(h,q)for zs in raw.values()for h,q,c in zs}=={(23,272),(9,288),(27,288),(23,296)})
    need(all(105<=c<=137 and c>99 for zs in raw.values()for h,q,c in zs if(h,q)==(23,272)))
    need(all((23,272//2**r)not in P for r in range(1,9) if 272%2**r==0))
    need(set(actual)<=upper)
    need(len([n for n in range(1,415,2)if n%18==9 or n%46 in(25,41)])==39)
    need(len([n for n in range(1,415,2)if n%18==9 or n%46==25])==31)
    need(not any({n,n+2,n+4}<=comb for n in range(lo,hi-3,2)))
    seq=sorted(actual);need(max(b-a for a,b in zip(seq,seq[1:]))==44 and {29833,29835}<=actual.keys())
    Q=sorted(D)
    ratios=[(Fraction(197*max(D[q])-min(D[r]),99*(r-q)),q,r)for q,r in zip(Q,Q[1:])]
    need(max(ratios)==(Fraction(662,99),288,296))
    gaps=[(99*(r-q)*8+min(D[r])-197*max(D[q]),q,r)for q,r in zip(Q,Q[1:])]
    need(min(gaps)==(622,104,108))
    direct=[]
    for t,h,q,c,N,X in ((2,13,104,123,42783,(2,12,9,17,7)),(2,15,108,1,42783,(16,5,9,3,14)),
                       (3,13,104,197,84929,(2,12,9,17,7)),(3,15,108,1,85551,(16,5,9,3,14))):
        need(pair_profile(X)==(h,q));z=B[str(c)];r=z['shift'];need(t>=r)
        U=ordinal(X,tt(t-r));T=star(tuple(z['I']),U,tuple(z['W']),z['ic'],z['wc']);validate(T)
        C=[S for S in components(T)if len(S)>1];need(len(C)==1 and all(len(z['I'])+j in C[0]for j in range(len(U))))
        core=induce(T,C[0]);need(strong(core) and count(core)==N and N==99*q*2**t+c*h)
        direct.append(dict(t=t,h=h,q=q,c=c,phase=r,N=N,order=len(core),rows=core))
    for t in range(3,30):
        L=10296*2**t+2561;R=10692*2**t+15
        need((R-L)//2-1==198*2**t-1274)
        for q,hs in D.items():
            need((99*q*2**t+197*max(hs)<=L)or(99*q*2**t+min(hs)>=R))
    print(json.dumps(dict(status='PASS',physical_pairs=35,star_slack_certificates=99,phase_lower_bound_tight=67,
        relaxed_local_targets=len(raw),universal_star_upper_targets=len(upper),published_local_targets=len(actual),
        published_targets={str(n):v for n,v in actual.items()},max_adjacent_clock_ratio='662/99',
        sharp_direct_witnesses=direct,minimum_clock_gap_at_t3=622,first_omitted_odds=310),sort_keys=True,indent=2))
if __name__=='__main__':main()
