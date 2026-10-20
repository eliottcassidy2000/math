#!/usr/bin/env python3
"""Independent finite controls; no imports from the reviewed verification code.
Integer/Fraction assertions are exact. Complex/Decimal checks are numerical only.
"""
from collections import Counter, defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import product, combinations
from math import comb, gcd, lcm, log, pi, ceil
import cmath, hashlib, json, random
from pathlib import Path

BASE=Path(__file__).resolve().parent
counts=Counter(); numerical=Counter(); tol=1e-9

def check(ok,label):
    if not ok: raise RuntimeError(label)

def elements(mods): return list(product(*(range(m) for m in mods)))
def add(x,y,mods): return tuple((a+b)%m for a,b,m in zip(x,y,mods))
def sub(x,y,mods): return tuple((a-b)%m for a,b,m in zip(x,y,mods))
def phase(c,x,mods):
    e=lcm(*mods)
    return sum(a*b*(e//m) for a,b,m in zip(c,x,mods))%e

def subset_counts(seq,mods):
    zero=tuple(0 for _ in mods); ans=Counter()
    for mask in range(1<<len(seq)):
        s=zero
        for j,x in enumerate(seq):
            if mask>>j&1:s=add(s,x,mods)
        ans[(mask.bit_count(),s)]+=1
    return ans

# Exhaustive labelled lists, including noncyclic groups, independently enumerated.
for mods in [(2,),(3,),(2,2)]:
    G=elements(mods); exponent=lcm(*mods)
    for n in range(2,6):
        for seq in product(G,repeat=n):
            dist=subset_counts(seq,mods)
            for k in range(1,n):
                b=Fraction(comb(n,k)*k**k*(n-k)**(n-k),n**n)
                check(b>=Fraction(1,n+1),'exact binomial mode')
                check(sum(dist[k,g] for g in G)==comb(n,k),'exact coefficient total')
                for c in G[1:]:
                    z=[cmath.exp(2j*pi*phase(c,x,mods)/exponent) for x in seq]
                    S=sum(z); delta=n-abs(S)
                    coefficient=sum(dist[k,g]*cmath.exp(2j*pi*phase(c,g,mods)/exponent) for g in G)/comb(n,k)
                    bound=cmath.exp(-k*(n-k)/n**2*delta).real/float(b)
                    check(abs(coefficient)<=bound+tol,'F1 exhaustive numerical')
                    pair=sum(abs(z[i]-z[j])**2 for i in range(n) for j in range(i+1,n))
                    check(abs(pair-(n*n-abs(S)**2))<=tol,'F5 numerical')
                    numerical['exhaustive_product_group_F1_and_F5']+=1
                counts['binomial_mode_and_total']+=1
            counts['exhaustive_labelled_lists']+=1

# Genuinely nonvacuous F1 cases: exact C2 coefficients, 80-digit Decimal bound.
nonvacuous=0; strongest=None
with localcontext() as ctx:
    ctx.prec=80
    for n in [24,40,80,160,240]:
        for w in [n//4,n//2,3*n//4]:
            for k in [n//4,n//2,3*n//4]:
                signed=sum((-1)**j*comb(w,j)*comb(n-w,k-j)
                    for j in range(max(0,k-n+w),min(k,w)+1))
                total=comb(n,k)
                delta=n-abs(n-2*w)
                rho=Decimal(k)/Decimal(n)
                b=Decimal(total)*rho**k*(1-rho)**(n-k)
                rhs=(-rho*(1-rho)*Decimal(delta)).exp()/b
                lhs=Decimal(abs(signed))/Decimal(total)
                check(lhs<=rhs,'high-precision nonvacuous F1')
                if rhs<1:
                    nonvacuous+=1
                    if strongest is None or rhs<Decimal(strongest['bound']):
                        strongest={'n':n,'ones':w,'k':k,'actual':str(lhs),'bound':str(rhs)}
                if rhs<=Decimal('0.5'):
                    check((total+signed)*2>=total and (total-signed)*2>=total,'exact F3 multiplicity')
                    counts['effective_F3_multiplicity']+=1
                numerical['high_precision_binary_F1']+=1
check(nonvacuous>0,'nonvacuous tests exist')

def prune(seq,mods):
    G=elements(mods); H=set(G); a=G[0]; T=list(range(len(seq)))
    N=len(seq); D=max(1,(len(G)-1).bit_length()); history=[]
    while True:
        choice=None
        for c in G:
            if all(phase(c,h,mods)==0 for h in H):continue
            fibers=defaultdict(list)
            for j in T:fibers[phase(c,sub(seq[j],a,mods),mods)].append(j)
            keep=max(fibers.values(),key=len)
            if 2*D*(len(T)-len(keep))<N:
                choice=(c,keep);break
        if choice is None:break
        c,keep=choice
        old=len(H); H={h for h in H if phase(c,h,mods)==0}
        check(old>=2*len(H),'strict subgroup reduction')
        history.append(len(T)-len(keep));T=keep;a=seq[T[0]]
    check(2*len(T)>=N and len(history)<=D,'pruning depth and retention')
    check(all(sub(seq[j],a,mods) in H for j in T),'retained coset')
    for c in G:
        if all(phase(c,h,mods)==0 for h in H):continue
        fs=Counter(phase(c,sub(seq[j],a,mods),mods) for j in T)
        check(2*D*(len(T)-max(fs.values()))>=N,'terminal fibers')
    period=next(k for k in range(1,len(G)+1) if tuple(k*x%m for x,m in zip(a,mods)) in H)
    check((len(G)//len(H))%period==0,'period divides index')
    return a,H,T,period

rng=random.Random(820261019)
for mods in [(6,),(8,),(2,3),(2,2,2),(4,2)]:
    G=elements(mods)
    for trial in range(100):
        extras=[rng.choice(G) for _ in range(trial%3)]
        ambient=G+G+extras; A0=len(extras); Q=len(ambient)
        chosen=rng.sample(range(Q),rng.randint(1,min(Q,12)))
        seq=[ambient[j] for j in chosen]
        a,H,T,h=prune(seq,mods);n=len(T);d=len(G)//len(H)
        if n>A0:
            check(d*(n-A0)<=Q-A0,'exact F7 with nonzero Fourier upper bound')
            counts['nonuniform_ambient_index']+=1
        dist=subset_counts([seq[j] for j in T],mods)
        for (k,target),multiplicity in dist.items():
            coset_shift=tuple(k*x%m for x,m in zip(a,mods))
            check(sub(target,coset_shift,mods) in H,'every size has predicted affine coset')
            if target==tuple(0 for _ in mods):check(k%h==0,'identity size period necessity')
        counts['product_group_affine_pruning']+=1

# Exercise high-order Fejer branch, unlike tests confined below order 64.
for m in [127,257,509]:
    for R in [64,96]:
        for phi in [0.0,0.123,1.71,pi]:
            def fejer(theta):
                return abs(sum(cmath.exp(1j*j*theta) for j in range(R)))**2/R
            vals=[fejer(2*pi*x/m-phi) for x in range(m)]
            total=sum(vals)+fejer(-phi) # uniform cyclic ambient plus one extra zero
            expansion=m+fejer(-phi)
            check(abs(total-expansion)<tol,'translated Fejer exact-model total numerically')
            check(total<=m+R+tol,'nonzero ambient A0=1 bound')
            close=sum(abs(cmath.phase(cmath.exp(1j*(2*pi*x/m-phi))))<=1/R for x in range(m))
            close+=int(abs(cmath.phase(cmath.exp(-1j*phi)))<=1/R)
            check(close<=8*(m+1)/R,'translated arc count')
            for theta in [0,1/(2*R),1/R,-1/R]:
                check(fejer(theta)+tol>=R/4,'central Fejer floor')
            numerical['genuine_high_order_translated_Fejer']+=1

# Pascal codes constructed directly via Lucas parity, without matrix iteration.
def columns(F,t):
    return [sum((int((h&j)==h))<<h for h in range(F)) for j in range(t)]
def weights(F,P):
    cs=columns(F,P)
    return Counter(sum((c&v).bit_count()%2 for c in cs) for v in range(1<<F))
def kraw(n,w,k):
    return sum((-1)**j*comb(w,j)*comb(n-w,k-j) for j in range(max(0,k-n+w),min(w,k)+1))
def xor_counts(cs,F):
    out=[[0]*(1<<F) for _ in range(len(cs)+1)];out[0][0]=1
    for step,c in enumerate(cs):
        newer=[row[:] for row in out]
        for k in range(step+1):
            for s,n in enumerate(out[k]):newer[k+1][s^c]+=n
        out=newer
    return out

examples={}
for P in [1,2,4,8,16]:
    for F in range(P+1):
        W=weights(F,P)
        if P>1:
            L=P//2
            if F<=L:expected=Counter({2*w:n for w,n in weights(F,L).items()})
            else:
                expected=Counter()
                for w,n in weights(F-L,L).items():
                    for j in range(L-w+1):expected[w+2*j]+=n*2**w*comb(L-w,j)
            check(W==expected,'Pascal recursion')
        if P//2<F<=P:
            check(min(w for w in W if w)==(1 if F==P else 2),'exact Pascal code distance')
        if (F,P) in [(3,4),(5,8),(7,8)]:examples[str(F)]=dict(sorted(W.items()))
        counts['Pascal_weight_enumerators']+=1

for F in range(1,7):
    P=1<<(F-1).bit_length()
    for periods in [1,2]:
        t=P*periods;cs=columns(F,t);dist=xor_counts(cs,F)
        for k,row in enumerate(dist):
            check(all(not n or (s&1)==k%2 for s,n in enumerate(row)),'fixed-weight terminal coordinate')
            for v in range(1<<F):
                signed=sum(n*(-1)**((s&v).bit_count()%2) for s,n in enumerate(row))
                w=sum((c&v).bit_count()%2 for c in cs)
                check(signed==kraw(t,w,k),'fixed-weight automaton transform')
                counts['automaton_Krawtchouk']+=1
        for q in [Fraction(2,5),Fraction(2,3)]:
            a=1-2*q
            prob=[sum(Fraction(dist[k][s])*q**k*(1-q)**(t-k) for k in range(t+1)) for s in range(1<<F)]
            W=weights(F,P)
            value=sum(n*(a**periods)**w for w,n in W.items())/2**F
            check(prob[0]==value,'exact complete-period return law')
            tv=sum(abs(p-Fraction(1,2**F)) for p in prob)/2
            rhs=(sum(n*(a**(2*periods))**w for w,n in W.items())-1)/4
            check(tv*tv<=rhs,'exact squared TV inequality')
            counts['rational_biased_automaton_return_and_TV']+=1

# Tournament orientation counts calculated from outdegrees themselves.
for v in range(2,6):
    edges=list(combinations(range(v),2));e=len(edges);by_size=Counter()
    base=sum(((v-1-i)&1)<<i for i in range(v))
    for bits in product([0,1],repeat=e):
        degree=[0]*v
        for (i,j),bit in zip(edges,bits):degree[j if bit else i]+=1
        par=sum((d&1)<<i for i,d in enumerate(degree))
        by_size[sum(bits),par^base]+=1
    for z in range(1<<v):
        count=sum(by_size[k,z] for k in range(e+1))
        check(count==(0 if z.bit_count()%2 else 2**(e-v+1)),'tournament parity fibers')
    for s in range(1<<(v-1)):
        cut=s.bit_count()*(v-s.bit_count())
        for k in range(e+1):
            value=sum(by_size[k,z]*(-1)**((z&s).bit_count()%2) for z in range(1<<v))
            check(value==kraw(e,cut,k),'tournament cut transform')
            counts['tournament_cut_Krawtchouk']+=1

def rank(vec):
    piv={}
    for x in vec:
        while x:
            k=x.bit_length()
            if k not in piv:piv[k]=x;break
            x^=piv[k]
    return len(piv)
for v in range(2,6):
    edges=list(combinations(range(v),2))
    for mask in range(1,1<<len(edges)):
        es=[e for j,e in enumerate(edges) if mask>>j&1]
        cs=[(1<<i)|(1<<j) for i,j in es]
        if rank(cs)!=v-1:continue # incidence rank gives connectedness
        bip=any(all(((s>>i)^(s>>j))&1 for i,j in es) for s in range(1<<v))
        check(rank([c^cs[0] for c in cs])==v-2+int(not bip),'connected host affine rank')
        counts['connected_graph_affine_rank']+=1

# Binary digit lifts, including the counterexample to a pruning-based a-count.
for p in [2,3,5]:
    for _ in range(70):
        M=[[rng.randrange(p) for j in range(4)] for i in range(2)]
        lift=[rng.randrange(2*p) for j in range(4)]
        raw=[sum(x*y for x,y in zip(row,lift)) for row in M]
        if any(x%p for x in raw):continue
        target=tuple(-(x//p)%p for x in raw)
        syndromes=Counter();actual=Counter()
        for w in product([0,1],repeat=4):
            syndrome=tuple(sum(x*y for x,y in zip(row,w))%p for row in M)
            syndromes[sum(w),syndrome]+=1
            corrected=[x+p*y for x,y in zip(lift,w)]
            if all(sum(x*y for x,y in zip(row,corrected))%(p*p)==0 for row in M):actual[sum(w)]+=1
        check(all(actual[k]==syndromes[k,target] for k in range(5)),'binary correction compiler')
        counts['binary_correction_instances']+=1
check([w for w in [0,1] if (3+3*w)%9==0]==[],'odd-prime lift hostile')
pruning_hostile={'p':2,'columns':[0]*9+[1],'defect':0,'k':8,'actual':comb(9,8),'incorrect_unpruned_lower_bound':str(Fraction(comb(10,8),4))}
check(pruning_hostile['actual']<Fraction(comb(10,8),4),'retained-size hostile')

# Envelope: independent direct lcm of every smooth integer, using exact powers.
def primes(n):
    return [q for q in range(2,n+1) if all(q%d for d in range(2,int(q**0.5)+1))]
for X in range(2,101):
    ps=primes(X)
    for Y in range(2,X+1):
        allowed=[q for q in ps if q<=Y]
        E=1;layered=1
        for q in allowed:
            power=q
            while power*q<=X:power*=q
            E*=power
        for j in range(1,X.bit_length()+1):
            for q in allowed:
                if q**j<=X:layered*=q
        direct=1
        for n in range(1,X+1):
            residue=n
            for q in allowed:
                while residue%q==0:residue//=q
            if residue==1:direct=lcm(direct,n)
        check(E==direct==layered,'prime-power envelope and theta-layer identity')
        counts['exact_prime_power_envelopes']+=1

# Complementary Pascal duality and coefficientwise MacWilliams transform.
for P in [1,2,4,8,16]:
    rows=[sum(int((h&j)==h)<<j for j in range(P)) for h in range(P)]
    for F in range(P+1):
        check(all((rows[h]&rows[k]).bit_count()%2==0 for h in range(F) for k in range(P-F)), 'Pascal dual orthogonality')
        check(rank(rows[:F])==F and rank(rows[:P-F])==P-F, 'Pascal complementary dimensions')
        W=weights(F,P);dual=weights(P-F,P)
        for w in range(P+1):
            check(W[w]*2**(P-F)==sum(n*kraw(P,j,w) for j,n in dual.items()), 'coefficientwise MacWilliams')
        counts['Pascal_duality_and_MacWilliams']+=1

# Exact rational identity positivity, using direct probability propagation.
def bernoulli_distribution(seq, probabilities, mods):
    law={tuple(0 for _ in mods):Fraction(1)}
    for g,p in zip(seq,probabilities):
        new=defaultdict(Fraction)
        for z,mass in law.items():
            new[z]+=mass*(1-p);new[add(z,g,mods)]+=mass*p
        law=new
    return law
for mods in [(2,),(2,2),(2,2,2),(2,2,2,2)]:
    G=elements(mods)
    for trial in range(30):
        seq=[rng.choice(G) for _ in range(8)]
        probs=[rng.choice([Fraction(0),Fraction(1,5),Fraction(1,2)]) for _ in seq]
        law=bernoulli_distribution(seq,probs,mods);identity=law[G[0]]
        check(identity>=Fraction(1,len(G)) and all(identity>=mass for mass in law.values()), 'exponent-two identity positivity')
        counts['exact_exponent_two_positivity']+=1
for mods in [(3,),(5,),(8,),(2,3)]:
    G=elements(mods)
    for trial in range(30):
        seq=[];probs=[]
        for _ in range(4):
            g=rng.choice(G);p=rng.choice([Fraction(1,3),Fraction(2,3),Fraction(1)])
            seq.extend([g,tuple(-x%m for x,m in zip(g,mods))]);probs.extend([p,p])
        law=bernoulli_distribution(seq,probs,mods);identity=law[G[0]]
        check(identity>=Fraction(1,len(G)) and all(identity>=mass for mass in law.values()), 'inverse-paired identity positivity')
        counts['exact_inverse_paired_positivity']+=1
check(bernoulli_distribution([(1,),(1,)],[Fraction(1,2)]*2,(3,))[(0,)]==Fraction(1,4), 'nonpositive C3 hostile')

beta=Fraction(95,96);delta=sigma=Fraction(1,1024);eta=Fraction(1,128)
check(beta+2*delta+2*sigma<1 and 2*delta+2*sigma<eta<1-sigma,'conditional parameter sample')
gap=eta-2*delta-2*sigma
scales=[{'log_y':u,'log_of_ratio':float(gap)*u-2*log(u)} for u in [100,1000,10000,100000]]
root=BASE.parent
metadata=json.loads((root/'SOURCE_METADATA.json').read_text())
actual_source_hashes={}
for relative,expected in metadata['included_source_sha256'].items():
    actual=hashlib.sha256((root/relative).read_bytes()).hexdigest()
    check(actual==expected,'included source hash mismatch: '+relative)
    actual_source_hashes[relative]=actual
report={'status':'PASS','independence':'Standalone implementation; no execution/import of reviewed verification code.','exact_integer_or_fraction_checks':dict(counts),'numerical_checks':dict(numerical),'complex_absolute_tolerance':tol,'decimal_precision':80,'nonvacuous_F1_cases':nonvacuous,'strongest_F1_control':strongest,'Pascal_examples':examples,'pruning_wording_hostile':pruning_hostile,'log_budget_example':{'beta':str(beta),'delta':str(delta),'sigma':str(sigma),'eta':str(eta),'power_gap':str(gap),'exponent':str(1-sigma-eta),'scales':scales},'input_sha256':actual_source_hashes,'pinned_external_sources':metadata['external_sources'],'source_metadata_sha256':hashlib.sha256((root/'SOURCE_METADATA.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Finite controls do not prove analytic (S),(A), asymptotic prime supply, Lean correctness, or novelty.'}
(BASE/'controls.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
