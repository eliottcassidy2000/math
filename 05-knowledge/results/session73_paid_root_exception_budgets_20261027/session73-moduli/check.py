#!/usr/bin/env python3
"""Exact finite controls for REPORT.md; explicit checks survive python -O."""
from fractions import Fraction as F
from itertools import product, combinations, combinations_with_replacement
from math import gcd, lcm, isqrt
from functools import reduce
import hashlib
import json
import random


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def dot(a, b):
    return sum((x*y for x,y in zip(a,b)), F(0))


def weighted(a,b,mu):
    return sum((x*y*w for x,y,w in zip(a,b,mu)),F(0))


def rank(rows):
    if not rows:
        return 0
    a=[[F(x) for x in row] for row in rows]
    r=0
    for c in range(len(a[0])):
        p=next((i for i in range(r,len(a)) if a[i][c]),None)
        if p is None:
            continue
        a[r],a[p]=a[p],a[r]
        v=a[r][c]
        a[r]=[x/v for x in a[r]]
        for i in range(r+1,len(a)):
            if a[i][c]:
                v=a[i][c]
                a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a):
            break
    return r


def basis(rows):
    # Incremental echelon reduction; keep original independent rows.
    echelon=[]
    out=[]
    for raw in rows:
        v=list(map(F,raw))
        for c,e in echelon:
            t=v[c]
            if t:
                v=[x-t*y for x,y in zip(v,e)]
        c=next((j for j,x in enumerate(v) if x),None)
        if c is not None:
            d=v[c]
            v=[x/d for x in v]
            echelon.append((c,v))
            out.append(list(map(F,raw)))
    return out


def solve(a,b):
    n=len(a)
    need(n==len(b),"solve shape")
    if not n:
        return []
    x=[list(map(F,row))+[F(v)] for row,v in zip(a,b)]
    for c in range(n):
        p=next((j for j in range(c,n) if x[j][c]),None)
        need(p is not None,"singular solve")
        x[c],x[p]=x[p],x[c]
        t=x[c][c]
        x[c]=[z/t for z in x[c]]
        for i in range(n):
            if i!=c and x[i][c]:
                t=x[i][c]
                x[i]=[u-t*v for u,v in zip(x[i],x[c])]
    return [row[-1] for row in x]


def project(f,gs,mu):
    C=[[weighted(g,h,mu) for h in gs] for g in gs]
    y=[weighted(g,f,mu) for g in gs]
    alpha=solve(C,y)
    out=[sum((a*g[j] for a,g in zip(alpha,gs)),F(0)) for j in range(len(mu))]
    return out,weighted(out,out,mu)


def affine_period(word):
    A,D,Q=1,1,1
    for a,b,d,h,m in word:
        Q=lcm(Q,D*h//gcd(A,D*h))
        A*=a
        D*=d
        Q=lcm(Q,D*m//gcd(A,D*m))
    return Q


def affine_run(x,word,guard_residues):
    obs=[]
    for (a,b,d,h,m),r in zip(word,guard_residues):
        if x%h!=r:
            return None
        z=a*x+b
        if z%d:
            return None
        x=z//d
        obs.append(x%m)
    return tuple(obs)


def check_affine():
    rng=random.Random(731042)
    tested=0
    shifts=0
    reset_cases=0
    for _ in range(400):
        x0=rng.randrange(-12,13)
        x=x0
        word=[]
        guards=[]
        for _ in range(rng.randrange(1,5)):
            a=rng.choice([-3,-1,0,1,2,3])
            d=rng.randrange(1,4)
            h=rng.randrange(1,5)
            m=rng.randrange(1,6)
            nxt=rng.randrange(-12,13)
            b=d*nxt-a*x
            word.append((a,b,d,h,m))
            guards.append(x%h)
            x=nxt
        Q=affine_period(word)
        baseline=affine_run(x0,word,guards)
        need(baseline is not None,"constructed legal source")
        # All shifts in two periods, including negative representatives.
        for delta in range(-Q, Q+1):
            need((affine_run(x0+delta,word,guards)==baseline)==(delta%Q==0),
                 "least guarded affine period")
            shifts+=1
        tested+=1
        reset_cases+=int(any(a==0 for a,b,d,h,m in word))
    hostile=[(1,0,2,2,1),(2,0,1,1,3)]
    need(affine_period(hostile)==6,"cancelled guard hostile period")
    need(affine_run(2,hostile,[0,0])==(0,2),"hostile positive source")
    need(affine_run(5,hostile,[0,0]) is None,"hostile rejected source")
    return {"words":tested,"source_shifts":shifts,"with_zero_slope":reset_cases,
            "cancelled_endpoint_period":3,"guarded_period":6}


def cylinder(ks):
    h=len(ks)
    S=B=0
    for k in ks:
        B=3*B+2**S
        S+=k
    modulus=2**(S+1)
    c=((2**S-B)*pow(3**h,-1,modulus))%modulus
    d=(3**h*c+B)//2**S
    need(c%2==d%2==1,"odd cylinder endpoints")
    return (c-1)//2,(d-1)//2,S,h


def native(n,h):
    out=[]
    clock=0
    for _ in range(h):
        z=3*n+1
        k=0
        while z%2==0:
            z//=2
            k+=1
        out.append(k)
        clock+=1+k
        n=z
    return tuple(out),n,clock


def word_precision(ks,q):
    return 2**sum(ks)*q//gcd(3**len(ks),q)


def check_collatz():
    words=[()]
    for h in range(1,4):
        words.extend(product(range(1,5),repeat=h))
    cnt=transitions=images=compositions=0
    examples=[]
    for ks in words:
        c,d,S,h=cylinder(ks)
        for q in range(1,19):
            g=gcd(3**h,q)
            per=q//g
            need(word_precision(ks,q)==2**S*per,"source-index precision")
            observed=[]
            for t in range(2*q):
                z=c+2**S*t
                got,n,clock=native(2*z+1,h)
                need(got==tuple(ks),"native valuation word")
                need((n-1)//2==d+3**h*t,"literal endpoint address")
                need(clock==h+S,"physical integer clock")
                observed.append(((n-1)//2)%q)
                transitions+=h
            for r in range(q):
                good=[t for t in range(2*q) if observed[t]==r]
                need(bool(good)==((r-d)%g==0),"endpoint image guard")
                if good:
                    t0=0 if per==1 else ((r-d)//g*pow(3**h//g,-1,per))%per
                    need(good==[t for t in range(2*q) if t%per==t0],"exact endpoint fibre")
                images+=1
            # Direct finite minimality in the normalized branch parameter.
            need(all((3**h*j)%q for j in range(1,per)),"least endpoint period")
            for cut in range(h+1):
                u,v=ks[:cut],ks[cut:]
                need(word_precision(u,word_precision(v,q))==word_precision(ks,q),
                     "precision composition")
                cu,du,Su,hu=cylinder(u)
                cv,dv,Sv,hv=cylinder(v)
                # Compose the actual word at a legal common source, not just Q sizes.
                for t in [0,1]:
                    z=c+2**S*t
                    ku,mid,clocku=native(2*z+1,hu)
                    kv,end,clockv=native(mid,hv)
                    need(ku==tuple(u) and kv==tuple(v),"composed guards")
                    need((mid-1)//2%2**Sv==cv,"suffix source cylinder")
                    need(end==2*(d+3**h*t)+1,"composed affine endpoint")
                    need(clocku+clockv==h+S,"clock composition")
                compositions+=1
            cnt+=1
        if ks in [(),(1,),(1,2),(2,1,3)]:
            examples.append({"word":ks,"source_base":c,"endpoint_base":d,
                             "valuation_sum":S,"clock":h+S,
                             "precision_q9":word_precision(ks,9),
                             "precision_q27":word_precision(ks,27)})
    return {"words":len(words),"word_modulus_pairs":cnt,"native_transitions":transitions,
            "endpoint_image_tests":images,"guarded_split_compositions":compositions,
            "examples":examples}


def crt(a,q,b,m):
    g=gcd(q,m)
    if (b-a)%g:
        return None
    k=m//g
    t=0 if k==1 else ((b-a)//g*pow(q//g,-1,k))%k
    return (a+q*t)%lcm(q,m)


def rows_for(q,m):
    L=lcm(q,m)
    return [[int(x%q==r) for x in range(L)] for r in range(q)]+[
            [int(x%m==r) for x in range(L)] for r in range(m)]


def check_crt_rank():
    cells=frames=0
    ranks=[]
    for q in range(1,13):
        for m in range(1,13):
            g,L=gcd(q,m),lcm(q,m)
            for a in range(q):
                for b in range(m):
                    got=crt(a,q,b,m)
                    actual=[x for x in range(L) if x%q==a and x%m==b]
                    need(actual==([] if got is None else [got]),"CRT exact intersection")
                    cells+=1
            gs=rows_for(q,m)
            rr=rank(gs)
            need(rr==q+m-g,"CRT incidence rank")
            need(L-rr==g*(q//g-1)*(m//g-1),"CRT cycle nullity")
            mu=[F(1,L)]*L
            h=[F(L)]+[F(0)]*(L-1)
            ph=[F(q*int(x%q==0)+m*int(x%m==0)-g*int(x%g==0)) for x in range(L)]
            u=[x-y for x,y in zip(h,ph)]
            need(all(weighted(row,u,mu)==0 for row in gs),"explicit residual zero marginals")
            kappa=F(L-q-m+g)
            need(weighted(u,u,mu)==kappa and u[0]==kappa,"uniform unresolved energy")
            if kappa:
                zero=[1-x/kappa for x in u]
                need(zero[0]==0 and min(zero)>=0,"nonnegative zero-root hostile")
                need(weighted(zero,zero,mu)==1+1/kappa,"sharp zero-root payment")
                need(all(weighted(row,zero,mu)==weighted(row,[F(1)]*L,mu) for row in gs),
                     "zero-root positive marginals")
            frames+=1
            if (q,m) in [(2,3),(4,6),(5,7),(8,12),(12,12)]:
                ranks.append({"q":q,"m":m,"L":L,"gcd":g,"rank":rr,"kernel":L-rr})
    return {"frames":frames,"cell_pairs":cells,"examples":ranks}


def resistance(q,m,mu,e):
    # Independent grounded Laplacian calculation, rather than a Gram inverse.
    N=q+m
    mat=[[F(0) for _ in range(N)] for _ in range(N)]
    for x,w in enumerate(mu):
        u=x%q
        v=q+x%m
        mat[u][u]+=w
        mat[v][v]+=w
        mat[u][v]-=w
        mat[v][u]-=w
    g=gcd(q,m)
    # One right-shore ground per component, labeled by residues 0..g-1.
    keep=[x for x in range(N) if x not in range(q,q+g)]
    b=[F(0)]*N
    b[e%q]=1
    b[q+e%m]=-1
    pot=solve([[mat[i][j] for j in keep] for i in keep],[b[i] for i in keep])
    return dot(pot,[b[i] for i in keep])


def check_weighted_gram():
    frames=targets=intervals=resistances=0
    examples=[]
    for q in range(1,7):
        for m in range(1,7):
            L=lcm(q,m)
            mu=[F(x+1,L*(L+1)//2) for x in range(L)]
            rows=rows_for(q,m)
            gs=basis(rows)
            need(len(gs)==q+m-gcd(q,m),"independent observation basis")
            for e in sorted({0,L//2,L-1}):
                h=[F(int(x==e),mu[e]) for x in range(L)]
                ph,pe=project(h,gs,mu)
                u=[x-y for x,y in zip(h,ph)]
                kappa=1/mu[e]-pe
                need(kappa>=0,"positive unresolved energy")
                need(weighted(u,u,mu)==kappa and u[e]==kappa,"source residual energy")
                need(all(weighted(row,u,mu)==0 for row in rows),"all redundant observations retained")
                R=resistance(q,m,mu,e)
                need(kappa==1/mu[e]-R,"weighted electrical identity")
                resistances+=1
                for seed in range(3):
                    f=[F(((x+1)*(seed+3))%11-5,seed+1) for x in range(L)]
                    f0,E0=project(f,gs,mu)
                    G=weighted(f,f,mu)
                    need((f[e]-f0[e])**2<=kappa*(G-E0),"sharp feasible interval")
                    need(all(weighted(row,f0,mu)==weighted(row,f,mu) for row in rows),"least-norm data fit")
                    for t in [F(-2),F(-1,2),F(0),F(1,3),F(2)]:
                        witness=[x+t*y for x,y in zip(f0,u)]
                        need(all(weighted(row,witness,mu)==weighted(row,f,mu) for row in rows),
                             "extremal witness preserves all data")
                        need(weighted(witness,witness,mu)==E0+t*t*kappa,"extremal exact energy")
                        need(witness[e]==f0[e]+t*kappa,"extremal exact evaluation")
                        intervals+=1
                if (q,m,e)==(2,3,0):
                    examples.append({"q":q,"m":m,"source":e,"mu":str(mu[e]),
                                     "resistance":str(R),"kappa":str(kappa)})
                targets+=1
            frames+=1
    # Cell-boundary algebra with independently supplied positive source/alias mass.
    onecell=0
    for mu,E,S,G in product([F(1,7),F(2,5),F(1)], [F(1,9),F(1,2),F(2)],
                             [F(-2),F(0),F(1,3),F(2)], [F(1,3),F(1),F(5),F(17)]):
        M=mu+E
        E0=S*S/M
        if G<E0:
            continue
        a=S/M
        kappa=E/(mu*M)
        need((a>0 and a*a>kappa*(G-E0))==(S>0 and S*S>E*G),"one-cell sharp threshold")
        onecell+=1
    return {"nonuniform_frames":frames,"named_targets":targets,
            "exact_extremal_checks":intervals,"independent_laplacians":resistances,
            "onecell_threshold_checks":onecell,"examples":examples}


def check_original_weights():
    ks=(1,2)
    c,d,S,h=cylinder(ks)
    # Explicit finite synthetic original prior. Never asserted to be native W or lambda.
    N=256
    weights=[F(z+1,N*(N+1)//2) for z in range(N)]
    direct={}
    affine={}
    direct_clock={}
    for z,weight in enumerate(weights):
        word,end,clock=native(2*z+1,h)
        if word!=ks:
            continue
        w=(end-1)//2
        key=(w%6,w%10)
        direct[key]=direct.get(key,F(0))+weight
        direct_clock.setdefault(key,set()).add(clock)
    for z,weight in enumerate(weights):
        if (z-c)%2**S:
            continue
        t=(z-c)//2**S
        w=d+3**h*t
        a,b=w%6,w%10
        both=crt(a,6,b,10)
        need(both==w%30,"shared endpoint CRT address")
        affine[(a,b)]=affine.get((a,b),F(0))+weight
    need(direct==affine,"original weighted endpoint pushforward")
    need(all(v=={h+S} for v in direct_clock.values()),"weighted bins retain physical clock")
    masses=list(direct.values())
    need(len(set(masses))>1,"nonuniform posterior hostile really nonuniform")
    need(sum(masses)==sum(w for z,w in enumerate(weights) if (z-c)%2**S==0),"mass conservation")
    # Dropping a product guard admits an incompatible pair.
    need(crt(0,6,1,10) is None,"incompatible shared endpoint rejected")
    return {"original_sources":N,"word":ks,"joint_cells":len(direct),
            "physical_clock":h+S,"accepted_source_mass":str(sum(masses)),
            "joint_weight_digest":hashlib.sha256(json.dumps(
                [[list(k),str(v)] for k,v in sorted(direct.items())],separators=(',',':')).encode()).hexdigest()}




def bank_terms(bank):
    out=[]
    for k in range(1,len(bank)+1):
        for I in combinations(bank,k):
            out.append(((-1)**(k+1),reduce(gcd,I)))
    return out


def check_modulus_banks():
    banks=[]
    for k in range(1,4):
        banks.extend(combinations_with_replacement(range(1,9),k))
    banks.extend([(2,3,5,7),(2,4,6,8),(2,3,6),(2,2,3)])
    frames=0
    examples=[]
    for bank in banks:
        L=lcm(*bank)
        terms=bank_terms(bank)
        R=sum(s*q for s,q in terms)
        for N in [L,2*L]:
            rows=[[int(x%q==r) for x in range(N)] for q in bank for r in range(q)]
            need(rank(rows)==R,"finite bank inclusion-exclusion rank")
            mu=[F(1,N)]*N
            h=[F(N)]+[F(0)]*(N-1)
            ph=[F(sum(sign*q*int(x%q==0) for sign,q in terms)) for x in range(N)]
            u=[a-b for a,b in zip(h,ph)]
            need(all(weighted(g,u,mu)==0 for g in rows),"bank projection has zero marginals")
            need(weighted(ph,ph,mu)==R,"bank projection energy")
            need(weighted(u,u,mu)==N-R,"bank residual energy")
            joint=[F(L*int(x%L==0)) for x in range(N)]
            alias=[a-b for a,b in zip(h,joint)]
            alignment=[a-b for a,b in zip(joint,ph)]
            need(weighted(alias,alias,mu)==N-L,"bank alias debt")
            need(weighted(alignment,alignment,mu)==L-R,"bank alignment debt")
            need(weighted(alias,alignment,mu)==0,"orthogonal two-debt split")
            if N>R:
                z=[1-a/(N-R) for a in u]
                need(z[0]==0 and min(z)>=0,"bank nonnegative zero-source witness")
                need(weighted(z,z,mu)==1+F(1,N-R),"bank sharp threshold")
            frames+=1
            if bank in [(2,3,5),(2,3,6),(2,2,3),(4,6),(2,3,5,7)] and N==L:
                examples.append({"bank":bank,"N":N,"rank":R,"kappa":N-R})
    weighted_frames=0
    for q in range(1,6):
        for m in range(1,6):
            L=lcm(q,m)
            N=2*L
            mu=[F(x+1,N*(N+1)//2) for x in range(N)]
            rows=[[int(x%a==r) for x in range(N)] for a in [q,m] for r in range(a)]
            gs=basis(rows)
            e=L//2
            h=[F(int(x==e),mu[e]) for x in range(N)]
            ph,pe=project(h,gs,mu)
            cell=[F(int(x%L==e%L)) for x in range(N)]
            M=dot(cell,mu)
            joint=[x/M for x in cell]
            alias=[a-b for a,b in zip(h,joint)]
            alignment=[a-b for a,b in zip(joint,ph)]
            need(weighted(alias,alias,mu)==1/mu[e]-1/M,"original-weight alias debt")
            need(weighted(alignment,alignment,mu)==1/M-pe,"original-weight incidence debt")
            need(weighted(alias,alignment,mu)==0,"original-weight debt orthogonality")
            repaired,_=project(h,basis(rows+[cell]),mu)
            need(repaired==joint,"single named joint-cell repair")
            weighted_frames+=1
    return {"nonempty_banks":len(banks),"uniform_fixed_N_frames":frames,
            "nonuniform_two_debt_frames":weighted_frames,"examples":examples,
            "empty_bank":"not in theorem domain"}

def isprime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))


def check_incoming_patterns():
    cases=odd_cases=aliases=0
    prime_examples=[]
    for p in range(3,200,2):
        if not isprime(p):
            continue
        for d in range(1,61):
            S={pow(x,d,p) for x in range(p)}-{0}
            neg={(-x)%p for x in S}
            tournament=(not(S&neg) and len(S)==(p-1)//2)
            need(tournament==(gcd(d,p-1)==2 and p%4==3),"prime monomial criterion")
            cases+=1
        if p<121 and p%4==3:
            squares={x*x%p for x in range(p)}
            for d in [3,4,5,13]:
                other={(x*x+p*x**d)%p for x in range(p)}
                need(other==squares,"arbitrary polynomial degree aliases")
                aliases+=1
        if p in [3,7,11,19]:
            prime_examples.append({"p":p,"fourth_power_nonzero":sorted({pow(x,4,p) for x in range(p)}-{0})})
    for m in range(3,121):
        for d in range(1,60,2):
            S={pow(x,d,m) for x in range(m)}-{0}
            need(1 in S and m-1 in S,"odd-degree opposite pair")
            odd_cases+=1
    S={x*x%15 for x in range(15)}-{0}
    U={x*x%15 for x in range(15) if gcd(x,15)==1}
    need(14 not in U and 6 in S and 9 in S,"composite nonunit square hostile")
    need({pow(x,4,7) for x in range(7)}-{0}=={1,2,4},"HYP9419 counterexample")
    return {"status":"HYP-9419(4) polynomial exclusion REFUTED",
            "incoming_commit":"64a7985e7054821c5958dc85f6dd489c844111f2",
            "prime_monomial_cases":cases,"odd_degree_opposite_controls":odd_cases,
            "polynomial_degree_aliases":aliases,"examples":prime_examples,
            "mod15_nonzero_squares":sorted(S),"mod15_unit_squares":sorted(U)}

def main():
    out={"status":"PASS","arithmetic":"exact integer/Fraction; no removable assertions",
         "scope":"finite controls supplement the scoped proofs; no native ROOT positivity claim",
         "affine":check_affine(),"collatz":check_collatz(),"crt_rank":check_crt_rank(),
         "weighted_gram":check_weighted_gram(),"original_weights":check_original_weights(),
         "incoming_patterns":check_incoming_patterns(),
         "modulus_banks":check_modulus_banks()}
    print(json.dumps(out,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
