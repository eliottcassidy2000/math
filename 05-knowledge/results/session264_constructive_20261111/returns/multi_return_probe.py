#!/usr/bin/env python3
"""Literal multi-return constructors and exact square-zero seam algebra.
Standard library only. All tests concern the newly constructed seams.
"""
from collections import Counter
from fractions import Fraction
from itertools import permutations
from pathlib import Path
import hashlib, json, random

# A square-zero polynomial is {vertex_mask: rational coefficient}.
def clean(p): return {m:c for m,c in p.items() if c}
def add(*ps):
    q=Counter()
    for p in ps:
        for m,c in p.items(): q[m]+=c
    return clean(q)
def scale(p,c): return clean({m:v*c for m,v in p.items()})
def mul(p,q):
    r=Counter()
    for a,x in p.items():
        for b,y in q.items():
            if not a&b: r[a|b]+=x*y
    return clean(r)
def odd(p): return {m:c for m,c in p.items() if m.bit_count()%2}
def negate_vars(p): return {m:c*(-1)**m.bit_count() for m,c in p.items()}
def exp_nil(p,n):
    ans={0:1}; power={0:1}
    for j in range(1,n+1):
        power=scale(mul(power,p),Fraction(1,j))
        if not power: break
        ans=add(ans,power)
    return ans
def inverse(p,n):
    assert p.get(0)==1
    q=scale(add(p,{0:-1}),-1); power={0:1}; ans={0:1}
    for _ in range(n):
        power=mul(power,q)
        if not power: break
        ans=add(ans,power)
    return ans
def mmul(A,B):
    k=len(A)
    return [[add(*(mul(A[i][a],B[a][j]) for a in range(k))) for j in range(k)] for i in range(k)]
def determinant(A):
    k=len(A); ans={}
    for p in permutations(range(k)):
        term={0:(-1)**sum(p[i]>p[j] for i in range(k) for j in range(i+1,k))}
        for i in range(k): term=mul(term,A[i][p[i]])
        ans=add(ans,term)
    return ans

def tournament(n,bits):
    A=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if bits>>k&1: A[i]|=1<<j
            else: A[j]|=1<<i
            k+=1
    return tuple(A)
def transitive(n): return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def join(A,B,ports):
    a=len(A); b=len(B)
    T=[x|(((1<<b)-1)<<a) for x in A]+[x<<a for x in B]
    assert len(set(i for i,j in ports))==len(ports)==len(set(j for i,j in ports))
    for i,j in ports: T[i]&=~(1<<(a+j)); T[a+j]|=1<<i
    return tuple(T)
def paths(A,start=None,end=None):
    ans=[]
    def walk(p,mask):
        if end is None or p[-1]==end: ans.append((p,mask))
        if end is not None and p[-1]==end: return
        for v in range(len(A)):
            if not mask>>v&1 and A[p[-1]]>>v&1: walk(p+(v,),mask|1<<v)
    for v in (range(len(A)) if start is None else [start]): walk((v,),1<<v)
    return ans
def cycles(A):
    ans=[]
    for s in range(len(A)):
        def walk(p,mask):
            if len(p)>=3 and len(p)%2 and A[p[-1]]>>s&1: ans.append((p,mask))
            for v in range(s+1,len(A)):
                if not mask>>v&1 and A[p[-1]]>>v&1: walk(p+(v,),mask|1<<v)
        walk((s,),1<<s)
    return ans
def cycle_poly(C):
    q=Counter()
    for _,m in C: q[m]+=1
    return dict(q)
def packing(C):
    P={(0,0):1}
    for _,m in C:
        Q=Counter(P)
        for (s,k),c in P.items():
            if not s&m: Q[s|m,k+1]+=c
        P=dict(Q)
    return P
def packed_at(P,t):
    q=Counter()
    for (m,k),c in P.items(): q[m]+=c*t**k
    return clean(q)
def pack_mul(P,Q):
    ans=Counter()
    for (a,k),x in P.items():
        for (b,l),y in Q.items():
            if not a&b: ans[a|b,k+l]+=x*y
    return dict(ans)
def packing_exponential(C,n):
    P={(0,0):1}; p={0:1}
    for k in range(1,n//3+1):
        p=scale(mul(p,C),Fraction(1,k))
        for m,c in p.items(): P[m,k]=c
    return P

def seam_matrix(A,B,ports,delete_forbidden=True):
    a=len(A); forbidden=set(ports); k=len(ports)
    starts=[paths(A,start=i) for i,j in ports]
    ends=[paths(B,end=j) for i,j in ports]
    M=[[Counter() for _ in range(k)] for _ in range(k)]
    for i in range(k):
        for j in range(k):
            for p,s in starts[i]:
                for q,t in ends[j]:
                    if delete_forbidden and (p[-1],q[0]) in forbidden: continue
                    M[i][j][s|(t<<a)]+=1
    return [[dict(p) for p in row] for row in M]
def literal_cross_cycles(A,B,ports):
    a=len(A); k=len(ports); forbidden=set(ports)
    legs=[[[] for _ in range(k)] for _ in range(k)]
    for i,(alpha,_) in enumerate(ports):
        for j,(_,beta) in enumerate(ports):
            for p,pm in paths(A,start=alpha):
                for q,qm in paths(B,end=beta):
                    if (p[-1],q[0]) not in forbidden:
                        legs[i][j].append((p+tuple(v+a for v in q),pm|(qm<<a)))
    result=[]
    for s in range(k):
        def walk(i,seen,word,used):
            for j in range(s,k):
                if j!=s and j in seen: continue
                for leg,mask in legs[i][j]:
                    if used&mask: continue
                    w=word+leg
                    if j==s:
                        if len(w)%2:
                            r=w.index(min(w)); result.append(w[r:]+w[:r])
                    else: walk(j,seen|{j},w,used|mask)
        walk(s,{s},(),0)
    assert len(set(result))==len(result)
    return set(result)

def cross_inventory(M):
    k=len(M); P=M; C={}; by_returns=[]
    for ell in range(1,k+1):
        q=odd(scale(add(*(P[i][i] for i in range(k))),Fraction(1,ell)))
        assert all(c.denominator==1 if isinstance(c,Fraction) else True for c in q.values())
        by_returns.append(q); C=add(C,q); P=mmul(P,M)
    return C,by_returns

def path_counts(A):
    n=len(A); h=[0]*(1<<n); h[0]=1; dp={}
    for mask in range(1,1<<n):
        row=[0]*n
        for v in range(n):
            if mask>>v&1:
                prev=mask^(1<<v)
                row[v]=1 if not prev else sum(dp[prev][u] for u in range(n) if prev>>u&1 and A[u]>>v&1)
        dp[mask]=row; h[mask]=sum(row)
    return h
def q_direct(A,r):
    H=path_counts(A); f=[0]*len(H); f[0]=1
    for _ in range(r):
        g=[0]*len(H)
        for mask in range(len(H)):
            sub=mask
            while True:
                g[mask]+=H[sub]*f[mask^sub]
                if not sub: break
                sub=(sub-1)&mask
        f=g
    return f[-1]
def q_from_pack(P,n,r): return sum(c*2**k*r**(n-m.bit_count()+k) for (m,k),c in P.items())
def ints(p): return {str(m):int(c) for m,c in sorted(p.items())}
def check(A,B,ports,full=False):
    n=len(A)+len(B); T=join(A,B,ports); C=cycles(T); M=seam_matrix(A,B,ports)
    cross,by=cross_inventory(M)
    oldA=cycles(A); oldB=[(tuple(v+len(A) for v in p),m<<len(A)) for p,m in cycles(B)]
    direct_cross=[(p,m) for p,m in C if m&((1<<len(A))-1) and m>>len(A)]
    assert cross==cycle_poly(direct_cross)
    assert literal_cross_cycles(A,B,ports)=={p for p,m in direct_cross}
    # Independently label each actual cycle by its actual number of reversed seams.
    for ell,pred in enumerate(by,1):
        direct=[(p,m) for p,m in direct_cross if sum(p[i]>=len(A) and p[(i+1)%len(p)]<len(A) for i in range(len(p)))==ell]
        assert pred==cycle_poly(direct)
    old=pack_mul(packing(oldA),packing(oldB))
    predicted=pack_mul(old,packing_exponential(cross,n))
    actual=packing(C)
    assert predicted==actual
    assert q_from_pack(predicted,n,1)==path_counts(T)[-1]
    if len(ports)==2:
        M11,M12=M[0]; M21,M22=M[1]
        singles=add(odd(M11),odd(M22),odd(mul(M12,M21)))
        pair=mul(odd(M11),odd(M22))
        direct={(0,0):1}
        direct.update({(m,1):c for m,c in singles.items()})
        direct.update({(m,2):c for m,c in pair.items()})
        assert direct==packing_exponential(cross,n)
    if full:
        k=len(M)
        minus=[[add({0:int(i==j)},scale(M[i][j],-1)) for j in range(k)] for i in range(k)]
        plus=[[add({0:int(i==j)},scale(negate_vars(M[i][j]),-1)) for j in range(k)] for i in range(k)]
        ratio=mul(determinant(plus),inverse(determinant(minus),n))
        assert ratio==exp_nil(scale(cross,2),n)
        for r in (1,2,3):
            rp={0:1}
            for _ in range(r): rp=mul(rp,ratio)
            allpack=mul(packed_at(old,2*r),rp)
            extraction=sum(c*r**(n-m.bit_count()) for m,c in allpack.items())
            assert extraction==q_direct(T,r)==q_from_pack(actual,n,r)
    return {'A':A,'B':B,'ports':ports,'rows':T,'odd_cycles':len(C),'cross_cycles_by_return_count':[int(sum(p.values())) for p in by],'packing_size_coefficients':dict(sorted(Counter({k:sum(c for (m,l),c in actual.items() if l==k) for k in range(n//3+1)}).items())),'H':path_counts(T)[-1],'Q2':q_from_pack(actual,n,2),'Q3':q_from_pack(actual,n,3),'M_support_coefficients':[[ints(p) for p in row] for row in M]}

def main():
    totals=Counter(); examples=[]
    # All labelled shores of orders 1..3, prefix ports, every pairing.
    for a in range(1,4):
        for b in range(1,4):
            for ia in range(1<<(a*(a-1)//2)):
                for ib in range(1<<(b*(b-1)//2)):
                    A=tournament(a,ia); B=tournament(b,ib)
                    for k in range(1,min(a,b)+1):
                        for bp in permutations(range(k)):
                            check(A,B,list(zip(range(k),bp)))
                            totals[str(k)]+=1
    rng=random.Random(264)
    for trial in range(24):
        a=4; b=3 if trial<12 else 4; k=2+trial%2
        A=tournament(a,rng.randrange(1<<(a*(a-1)//2))); B=tournament(b,rng.randrange(1<<(b*(b-1)//2)))
        ports=list(zip(rng.sample(range(a),k),rng.sample(range(b),k)))
        r=check(A,B,ports,full=True); totals['random_7_8']+=1
        if trial<4: examples.append(r)
    tri=(2,4,1)
    for ports in ([(0,0)],[(0,0),(1,1)],[(0,1),(1,0)],[(0,0),(1,1),(2,2)]):
        examples.append(check(tri,tri,ports,full=True))
    # A hostile to the tempting product-of-endpoint-paths seam formula.
    hostA=transitive(3); hostB=transitive(3); ports=[(0,2),(2,0)]
    good=check(hostA,hostB,ports,full=True)
    naive,_=cross_inventory(seam_matrix(hostA,hostB,ports,delete_forbidden=False))
    correct,_=cross_inventory(seam_matrix(hostA,hostB,ports))
    assert naive!=correct
    hostile={'actual':good,'naive_cross_count':int(sum(naive.values())),'correct_cross_count':int(sum(correct.values())),'spurious_supports':ints(add(naive,scale(correct,-1)))}
    # The two explicit six-vertex gates are a physical rooted dipole.
    trade=[]
    for ports in [[(0,0),(1,1)],[(0,1),(1,0)]]:
        T=join(tri,tri,ports); H=path_counts(T); C=cycles(T); P=packing(C); full=(1<<6)-1; deleted=full^1
        U=sum(H[full^m] for p,m in paths(T,start=0))
        V=sum(H[full^m] for p,m in paths(T,end=0))
        d=sum((-1)**(len(p)+1)*H[full^m] for p,m in paths(T,start=0))
        e=sum((-1)**(len(p)+1)*H[full^m] for p,m in paths(T,end=0))
        qdel=sum(H[s]*H[deleted^s] for s in range(64) if not s&~deleted)
        kdel=sum((-1)**s.bit_count()*H[s]*H[deleted^s] for s in range(64) if not s&~deleted)
        Z=Counter()
        for (m,k),c in P.items(): Z[k,m.bit_count()]+=c
        deletion_profiles=[]
        for z in range(6):
            qz=Counter()
            for (m,j),c in packing([(p,m) for p,m in C if not m>>z&1]).items():
                qz[5-m.bit_count()+j]+=c*2**j
            deletion_profiles.append(dict(sorted(qz.items())))
        qp=Counter()
        for (m,k),c in P.items(): qp[6-m.bit_count()+k]+=c*2**k
        trade.append({'rows':T,'cycles':[p for p,m in C],'compatible_address_pairs':[(i+1,j+1) for i in range(len(C)) for j in range(i+1,len(C)) if not C[i][1]&C[j][1]],'Z':{f'{k},{n}':v for (k,n),v in sorted(Z.items())},'Q_coefficients':dict(sorted(qp.items())),'vertex_deleted_Q_coefficients':deletion_profiles,'root0_state':{'h':H[deleted],'Q':qdel,'k':kdel,'H':H[-1],'U':U,'V':V,'d':d,'e':e},'return_singleton_H':path_counts(join(T,(0,),[(0,0)]))[-1]})
    assert trade[0]['Q_coefficients']==trade[1]['Q_coefficients']=={2:16,4:12,6:1}
    assert trade[0]['vertex_deleted_Q_coefficients'][0]==trade[1]['vertex_deleted_Q_coefficients'][0]=={1:2,3:6,5:1}
    assert trade[0]['root0_state']=={'h':9,'Q':84,'k':0,'H':29,'U':52,'V':41,'d':4,'e':5}
    assert trade[1]['root0_state']=={'h':9,'Q':84,'k':0,'H':29,'U':50,'V':43,'d':4,'e':5}
    assert [t['return_singleton_H'] for t in trade]==[77,75]
    result={'status':'FINITE-EXACT','rooted_dipole_trade':trade,'scope':'All labelled shores of orders 1..3, prefix marks and every pairing; 24 seeded orders 7/8; explicit controls. Not an all-tournament census.','counts':dict(totals),'examples':examples,'forbidden_seam_hostile':hostile,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name('multi_return_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'counts':dict(totals),'examples':[{k:v for k,v in x.items() if k!='M_support_coefficients'} for x in examples],'hostile':{k:v for k,v in hostile.items() if k!='actual'}},indent=2))
if __name__=='__main__': main()
