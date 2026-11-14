"""S328 optimal complete phase companions. Literal helpers adapted from S325."""
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
from pathlib import Path
from math import gcd,isqrt
from itertools import combinations
import json,hashlib

def pair_profile(T):
    H=[1]+[sum(e)for e in endtable(T)[1:]]
    return H[-1],sum(a*b for a,b in zip(H,H[::-1]))
def strong(T):return len(components(T))==1

def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

def main():
    root=Path(__file__).resolve().parent;lines=(root/'catalog4096.txt').read_text().splitlines();P={};literal=0
    for line in lines:
        q,h,n,*T=map(int,line.split());T=tuple(T)
        need(len(T)==n and q<=4096 and q>=2**n);validate(T)
        need(pair_profile(T)==(h,q));P.setdefault(q,{})[h]=T;literal+=1
    need(literal==919 and sorted(P[4096])==[1,5,19,25,217,223])
    additive=[];sign=[]
    for q,hs in P.items():
        for a,b in combinations(sorted(hs),2):
            if gcd(a,99)!=1 or gcd(b,99)!=1:continue
            if (b-a)%396==198:additive.append((q,a,b))
            if (a+b)%396==0:sign.append((q,a,b))
    need(additive==[(4096,19,217),(4096,25,223)] and not sign)
    # Exhaustive finite controls for the all-odd-slope stabilizer proof.
    for a in (3,5,9,11,33,99,101):
        S=set(range(1,2*a,2));O=set(range(1,4*a,2));stable=[];flip=[]
        for u in range(1,4*a,2):
            if gcd(u,a)>1:continue
            T={u*c%(4*a)for c in S}
            if T==S:stable.append(u)
            if T==O-S:flip.append(u)
        need(set(stable)=={1,2*a-1} and set(flip)=={2*a+1,4*a-1})
    C5=(2,4,9,3);C19=(2,20,9,3,45,15)
    S={19:ordinal(C19,tt(4)),25:ordinal(ordinal(C5,C5),tt(2)),217:(2,20,41,19,101,203,143,31),223:(2,36,49,7,75,217,143,31)}
    need(strong(C5) and strong(C19) and pair_profile(C5)==(5,32) and pair_profile(C19)==(19,256))
    for h,T in S.items():validate(T);need(pair_profile(T)==(h,4096))
    need(strong(S[217]) and strong(S[223]))
    bank=json.loads((root/'selected_templates.json').read_text());phase={int(c):r for c,r in json.loads((root/'phase_table.json').read_text()).items()}
    out=[]
    for h,c in ((19,1),(217,1),(25,1),(223,1),(25,5),(217,7)):
        z=bank[str(c)];need(z['shift']==0)
        I=tuple(z['I']);W=tuple(z['W']);hi,alpha,lam,_,_=profile(I,z['ic']);hw,_,_,beta,d=profile(W,z['wc'])
        need((alpha*beta,hi*d-hw*lam)==(99,c))
        T=star(I,S[h],W,z['ic'],z['wc']);validate(T);C=[p for p in components(T)if len(p)>1]
        need(len(C)==1 and all(len(I)+j in C[0]for j in range(len(S[h]))))
        U=induce(T,C[0]);N=405504+c*h
        need(strong(U) and count(U)==N)
        other={19:217,217:19,25:223,223:25}[h]
        digit=N*pow(other,-1,198)%198;need(digit==c)
        need((N-c*other)%99==0);qbad=(N-c*other)//99
        need(qbad%4==2)
        out.append(dict(H=h,c=c,N=N,order=len(U),other_H=other,impossible_Q=qbad,rows=U))
    need(prime(405529) and prime(407023) and 405629==7*57947 and prime(57947))
    flips=0
    for a,b in ((19,217),(25,223)):
        for c,r in phase.items():
            for k in (r,r+1,r+4):
                Q=4096*2**k
                for h,other in ((a,b),(b,a)):
                    N=99*Q+c*h;need(N*pow(other,-1,198)%198==c)
                    bad=(N-c*other)//99;need((N-c*other)%99==0 and bad%4==2);flips+=1
        need(gcd(a,b)==1 and b>197 and 99*4096>197*max(a,b))
    text=(root/'catalog4096.txt').read_bytes()
    need(text==(root/'catalog4096.incoming.txt').read_bytes())
    need((root/'catalog4096.log').read_bytes()==(root/'catalog4096.incoming.log').read_bytes())
    expected=[1,1,2,4,12,56,456,2160,716,151,18,1]
    got=[int(line.split()[2])for line in (root/'catalog4096.log').read_text().splitlines()if line.startswith('FINAL ')]
    need(got==expected and sum(got)==3578)
    print(json.dumps(dict(status='PASS',budget=4096,complete_class_counts=got,total_classes=sum(got),literal_profile_checks=literal,
        additive_companions=additive,sign_companions=sign,full_digit_phase_rejections=flips,direct_strong_outputs=out,
        catalog_sha256=hashlib.sha256(text).hexdigest()),sort_keys=True,indent=2))
if __name__=='__main__':main()
