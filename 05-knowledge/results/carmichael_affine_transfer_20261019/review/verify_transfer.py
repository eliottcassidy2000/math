#!/usr/bin/env python3
"""Exact enumerations and separately labelled numerical Fourier checks."""
import cmath
import hashlib
import itertools
import json
import math
import random
from collections import Counter
from fractions import Fraction
from pathlib import Path

rng = random.Random(105720261019)
checks = Counter()
tol = 2e-9

def require(cond, msg):
    if not cond:
        raise RuntimeError(msg)

def cyclic_dp(seq, m):
    out = [[0]*m for _ in range(len(seq)+1)]
    out[0][0] = 1
    for i,g in enumerate(seq):
        for k in range(i+1,0,-1):
            for h in range(m):
                out[k][(h+g)%m] += out[k-1][h]
    return out

def cyclic_brute(seq,m):
    out = [[0]*m for _ in range(len(seq)+1)]
    for mask in range(1<<len(seq)):
        out[mask.bit_count()][sum(g for i,g in enumerate(seq) if mask>>i&1)%m]+=1
    return out

# Exact count identities and numerical Fourier inequalities on a frozen sample.
min_numerical_slack = 1.0
for trial in range(600):
    m=rng.randrange(2,24)
    n=rng.randrange(2,11)
    seq=[rng.randrange(m) for _ in range(n)]
    dp=cyclic_dp(seq,m)
    require(dp==cyclic_brute(seq,m), "cyclic DP/brute disagreement")
    checks["exact_cyclic_list_comparisons"]+=1
    for k in range(1,n):
        rho=k/n
        b=math.comb(n,k)*rho**k*(1-rho)**(n-k)
        require(b+tol>=1/(n+1), "binomial mode bound")
        for ell in range(1,m):
            roots=[cmath.exp(2j*math.pi*ell*g/m) for g in seq]
            S=sum(roots)
            actual=abs(sum(dp[k][g]*cmath.exp(2j*math.pi*ell*g/m) for g in range(m))/math.comb(n,k))
            bound=math.exp(-rho*(1-rho)*(n-abs(S)))/b
            slack=bound-actual
            min_numerical_slack=min(min_numerical_slack,slack)
            require(slack>=-tol, "fixed-size dispersion bound")
            pair=sum(abs(z-w)**2 for i,z in enumerate(roots) for w in roots[i+1:])
            require(abs(pair-(n*n-abs(S)**2))<tol, "pair dispersion identity")
            checks["numerical_fixed_size_character_checks"]+=1

def prune_cyclic(seq,m):
    N=len(seq)
    D=max(1,(m-1).bit_length())
    T=list(range(N))
    step,a=1,0
    history=[]
    while True:
        h=m//step
        found=None
        for ell in range(1,h):
            fibers={}
            for i in T:
                require((seq[i]-a)%step==0,"affine containment")
                co=((seq[i]-a)//step)%h
                fibers.setdefault((ell*co)%h,[]).append(i)
            keep=max(fibers.values(),key=len)
            if (len(T)-len(keep))*2*D<N:
                found=(ell,keep)
                break
        if found is None:
            break
        ell,keep=found
        oldstep=step
        step*=h//math.gcd(ell,h)
        a=seq[keep[0]]
        history.append((oldstep,step,len(T)-len(keep)))
        T=keep
    require(2*len(T)>=N,"pruning retention")
    require(len(history)<=D,"pruning depth")
    for ell in range(1,m//step):
        c=Counter((ell*((seq[i]-a)//step))%(m//step) for i in T)
        require((len(T)-max(c.values()))*2*D>=N,"fiber pruning terminal condition")
    return a,step,T

for m in range(2,42):
    ambient=list(range(m))*3
    for trial in range(12):
        chosen=rng.sample(range(len(ambient)),rng.randrange(1,len(ambient)+1))
        seq=[ambient[i] for i in chosen]
        a,d,T=prune_cyclic(seq,m)
        require(len(T)*d<=len(ambient),"ambient index bound A0=0")
        D=max(1,(m-1).bit_length())
        R=math.ceil(64*len(ambient)/len(seq))
        lam=len(seq)/R**2*min(2/D,3/(4*math.pi**2))
        for ell in range(1,m//d):
            S=sum(cmath.exp(2j*math.pi*ell*((seq[i]-a)//d)/(m//d)) for i in T)
            require(len(T)-abs(S)+tol>=lam,"affine energy bound")
            checks["numerical_affine_energy_checks"]+=1
        checks["exact_affine_pruning_cases"]+=1

# Translation of Fejer kernels: ambient roots are exactly uniform.
for m in range(2,35):
    for R in range(2,min(17,m)+1):
        phi=math.pi*rng.random()
        vals=[abs(sum(cmath.exp(1j*j*(2*math.pi*g/m-phi)) for j in range(R)))**2/R for g in range(m)]
        require(abs(sum(vals)-m)<tol,"translated Fejer total")
        arc=sum(abs(cmath.phase(cmath.exp(1j*(2*math.pi*g/m-phi))))<=1/R for g in range(m))
        require(arc<=8*m/R+tol,"translated arc bound")
        checks["numerical_translated_fejer_cases"]+=1

def gf2rank(vectors):
    piv={}
    for x in vectors:
        while x:
            j=x.bit_length()-1
            if j in piv: x^=piv[j]
            else:
                piv[j]=x
                break
    return len(piv)

def control_columns(F,t):
    b=1<<(F-1)
    out=[]
    for j in range(t):
        out.append(b)
        b^=b>>1
    return out

def xor_dp(cols):
    F=max((x.bit_length() for x in cols),default=0)
    out=[[0]*(1<<F) for _ in range(len(cols)+1)]
    out[0][0]=1
    for i,c in enumerate(cols):
        for k in range(i+1,0,-1):
            for s in range(1<<F):
                out[k][s^c]+=out[k-1][s]
    return out

def enum_xor(cols,F):
    counts=[[0]*(1<<F) for _ in range(len(cols)+1)]
    for mask in range(1<<len(cols)):
        s=0
        for i,c in enumerate(cols):
            if mask>>i&1:s^=c
        counts[mask.bit_count()][s]+=1
    return counts

def kraw(k,w,n):
    return sum((-1)**j*math.comb(w,j)*math.comb(n-w,k-j)
               for j in range(max(0,k-(n-w)),min(w,k)+1))

def code_weight(F,P):
    if F==0:return Counter({0:1})
    cols=control_columns(F,P)
    return Counter(sum((v&c).bit_count()%2 for c in cols) for v in range(1<<F))

def Wrec(F,P):
    if F==0:return Counter({0:1})
    if P==1:return Counter({0:1,1:1})
    L=P//2
    if F<=L:return Counter({2*w:c for w,c in Wrec(F,L).items()})
    out=Counter()
    for w,c in Wrec(F-L,L).items():
        for j in range(L-w+1):
            out[w+2*j]+=c*(2**w)*math.comb(L-w,j)
    return out

automata_examples={}
for F in range(1,13):
    P=1<<(F-1).bit_length()
    W=code_weight(F,P)
    require(W==Wrec(F,P),"Pascal code recursion")
    require(min(w for w in W if w)>=(1 if F==P else 2),"code distance lower")
    require(min(w for w in W if w)==(1 if F==P else 2),"code distance exact")
    require(gf2rank(control_columns(F,F))==F,"control rank")
    b=control_columns(F,F)
    require(gf2rank([x^b[0] for x in b])==F-1,"control affine rank")
    for j,c in enumerate(control_columns(F,P)):
        expected=sum((math.comb(j,r)%2)<<(F-1-r) for r in range(F))
        require(c==expected,"Pascal versus matrix columns")
    if F in [3,5,7]:automata_examples[str(F)]={str(k):v for k,v in sorted(W.items())}
    checks["exact_pascal_codes"]+=1

for F in range(1,7):
    for t in range(F,min(F+5,12)):
        cols=control_columns(F,t)
        dp=xor_dp(cols)
        require(dp==enum_xor(cols,F),"control subset DP/brute")
        for k in range(t+1):
            require(all(c==0 or ((s>>(F-1))&1)==(k%2) for s,c in enumerate(dp[k])),"control weight parity")
            for v in range(1<<F):
                lhs=sum(c*((-1)**((v&s).bit_count()%2)) for s,c in enumerate(dp[k]))
                w=sum((v&c).bit_count()%2 for c in cols)
                require(lhs==kraw(k,w,t),"control Krawtchouk")
                checks["exact_control_krawtchouk_checks"]+=1
        # q=1/3, integer probabilities with denominator 3^t.
        probs=[sum(dp[k][s]*2**(t-k) for k in range(t+1)) for s in range(1<<F)]
        for v in range(1<<F):
            w=sum((v&c).bit_count()%2 for c in cols)
            transform=sum(c*((-1)**((v&s).bit_count()%2)) for s,c in enumerate(probs))
            require(transform==3**(t-w),"biased-control Fourier")
        checks["exact_control_laws"]+=1

tournament_examples={}
for v in range(2,7):
    edges=list(itertools.combinations(range(v),2))
    cols=[(1<<i)|(1<<j) for i,j in edges]
    m=len(cols)
    dp=xor_dp(cols)
    require(dp==enum_xor(cols,v),"tournament subset DP/brute")
    total=[sum(row[s] for row in dp) for s in range(1<<v)]
    feasible=1<<(m-v+1)
    require(all(c==(feasible if s.bit_count()%2==0 else 0) for s,c in enumerate(total)),"fair parity fibers")
    for S in range(1<<(v-1)):
        cut=S.bit_count()*(v-S.bit_count())
        for k in range(m+1):
            transform=sum(c*((-1)**((S&s).bit_count()%2)) for s,c in enumerate(dp[k]))
            require(transform==kraw(k,cut,m),"tournament Krawtchouk")
            checks["exact_tournament_krawtchouk_checks"]+=1
    tournament_examples[str(v)]={"edges":m,"parity_fiber":feasible,"parity_vectors":1<<(v-1)}
    checks["exact_complete_tournament_universes"]+=1

# Every connected labelled simple graph on 2..5 vertices.
for v in range(2,6):
    all_edges=list(itertools.combinations(range(v),2))
    for mask in range(1,1<<len(all_edges)):
        edges=[e for j,e in enumerate(all_edges) if mask>>j&1]
        adj=[[] for _ in range(v)]
        for i,j in edges:adj[i].append(j);adj[j].append(i)
        color={0:0}
        queue=[0]
        bip=True
        for i in queue:
            for j in adj[i]:
                if j not in color:color[j]=1-color[i];queue.append(j)
                elif color[j]==color[i]:bip=False
        if len(color)<v:continue
        cols=[(1<<i)|(1<<j) for i,j in edges]
        rank=gf2rank([c^cols[0] for c in cols])
        require(rank==v-2+int(not bip),"graph affine rank/bipartiteness")
        checks["exact_connected_graph_affine_ranks"]+=1

# Mod-p^2 arithmetic lifts: brute congruence versus subset-syndrome counts.
for p in [2,3,5]:
    completed=0
    while completed<100:
        rows=rng.randrange(1,4);a=rng.randrange(1,8)
        M=[[rng.randrange(-4,5) for _ in range(a)] for _ in range(rows)]
        v=[rng.randrange(-p,p*p) for _ in range(a)]
        Mv=[sum(c*x for c,x in zip(row,v)) for row in M]
        if any(z%p for z in Mv):continue
        d=[z//p%p for z in Mv]
        direct=Counter();syndrome=Counter()
        for mask in range(1<<a):
            w=[(mask>>i)&1 for i in range(a)]
            lifted=[v[i]+p*w[i] for i in range(a)]
            if all(sum(c*x for c,x in zip(row,lifted))%(p*p)==0 for row in M):direct[mask.bit_count()]+=1
            if all((sum(c*x for c,x in zip(row,w))+dd)%p==0 for row,dd in zip(M,d)):syndrome[mask.bit_count()]+=1
        require(direct==syndrome,"binary lift versus syndrome")
        completed+=1;checks["exact_binary_lifting_cases"]+=1

hostile_seq=[1]*8
hostile=cyclic_dp(hostile_seq,7)
require([i for i,c in enumerate(hostile[3]) if c]==[3],"constant-label affine hostile")
require((3+3*2)%9==0 and all((3+3*w)%9 for w in [0,1]),"odd-prime digit hostile")

# Complementary Pascal duals and their exact MacWilliams transforms.
for P in [1,2,4,8,16]:
    for F in range(P+1):
        C=code_weight(F,P)
        dual=code_weight(P-F,P)
        for h in range(F):
            for k in range(P-F):
                dot=sum((math.comb(j,h)%2)*(math.comb(j,k)%2) for j in range(P))%2
                require(dot==0,"Pascal complementary dual orthogonality")
        for w in range(P+1):
            numerator=sum(c*kraw(w,j,P) for j,c in dual.items())
            require(numerator==C[w]*(1<<(P-F)),"Pascal MacWilliams")
        checks["exact_pascal_duality_codes"]+=1

# Fourier positivity: independent exponent-two choices and inverse pairs.
for trial in range(100):
    F=rng.randrange(1,6)
    n=rng.randrange(1,10)
    cols=[rng.randrange(1<<F) for _ in range(n)]
    counts=enum_xor(cols,F)
    weighted=[sum(counts[k][x]*2**(n-k) for k in range(n+1)) for x in range(1<<F)]
    rank=gf2rank(cols)
    require(weighted[0]*(1<<rank)>=3**n,"binary identity-positive bound")
    require(weighted[0]==max(weighted),"binary modal identity")
    checks["exact_identity_positive_binary_cases"]+=1
for m in range(2,18):
    for trial in range(8):
        seed=[rng.randrange(m) for _ in range(4)]
        seq=[g for a in seed for g in (a,(-a)%m)]
        counts=cyclic_dp(seq,m)
        weighted=[sum(counts[k][x] for k in range(9)) for x in range(m)]
        require(weighted[0]*m>=2**8 and weighted[0]==max(weighted),"inverse-paired modal identity")
        checks["exact_inverse_paired_cases"]+=1
require(sum(row[0] for row in cyclic_dp([1,1],3))*3<4,"C3 positivity hostile")

pruning_hostile=cyclic_dp([0]*9+[1],2)
require(pruning_hostile[8][0]==9 and 9<Fraction(math.comb(10,8),4),"pruned reservoir entropy hostile")
checks["exact_pruned_entropy_hostile"]+=1

# Literal finite Korselt compiler: prime labels retained, no analytic hypotheses asserted.
def prime(n):
    return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))

carmichael_banks={}
for L in [80,720,1440]:
    bank=[d+1 for d in range(1,L+1) if L%d==0 and prime(d+1) and math.gcd(d+1,L)==1]
    direct=Counter()
    least=None
    dp=[[0]*L for _ in range(len(bank)+1)]
    dp[0][1]=1
    for i,p in enumerate(bank):
        for k in range(i+1,0,-1):
            for r in range(L):
                dp[k][r*p%L]+=dp[k-1][r]
    for mask in range(1<<len(bank)):
        if mask.bit_count()<2:continue
        ps=[p for i,p in enumerate(bank) if mask>>i&1]
        n=math.prod(ps)
        if n%L!=1:continue
        require(all((n-1)%(p-1)==0 for p in ps),"literal Korselt compiler")
        direct[len(ps)]+=1
        least=n if least is None else min(least,n)
    require(all(dp[k][1]==direct[k] for k in range(2,len(bank)+1)),"Korselt DP/brute")
    carmichael_banks[str(L)]={"prime_bank":bank,"counts_by_exact_size":dict(sorted(direct.items())),"smallest_constructed":least}
    checks["exact_finite_korselt_banks"]+=1

report={
 "status":"PASS",
 "exact_checks":{k:v for k,v in checks.items() if k.startswith("exact")},
 "numerical_checks":{"status":"VERIFIED floating checks, not exact proofs","tolerance":tol,
                     "counts":{k:v for k,v in checks.items() if k.startswith("numerical")},
                     "minimum_observed_F1_slack":min_numerical_slack},
 "automata_weight_enumerators":automata_examples,
 "tournament_counts":tournament_examples,
 "hostiles":{"constant_labels_C7":{"n":8,"k":3,"only_target":3,"count":math.comb(8,3)},
             "odd_prime_binary_digit":{"p":3,"matrix":[[1]],"lift":3,"unrestricted_correction":2,"binary_corrections":[]}},
 "literal_carmichael_banks":carmichael_banks,
 "seed":105720261019,
 "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
}
print(json.dumps(report,indent=2,sort_keys=True))
