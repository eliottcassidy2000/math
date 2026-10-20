#!/usr/bin/env python3
"""Self-contained exact controls; does not read an existing atlas or theorem output."""
from itertools import combinations, permutations, product
from collections import Counter, defaultdict
from fractions import Fraction
from math import factorial, gcd, comb
from pathlib import Path
import json
import sympy as sp

HERE=Path(__file__).resolve().parent

def partitions(n,lo=1):
    if n==0: yield ()
    for a in range(lo,n+1):
        for tail in partitions(n-a,a):yield (a,)+tail

def burnside(n):
    ordinary=Fraction(0);twisted=Fraction(0)
    for lam in partitions(n):
        c=sum(a//2 for a in lam)+sum(gcd(a,b) for a,b in combinations(lam,2))
        z=1
        for a,m in Counter(lam).items():z*=a**m*factorial(m)
        if all(a%2 for a in lam):ordinary+=Fraction(1<<c,z)
        if lam.count(1)<=1 and all(a==1 or a%4==2 for a in lam):twisted+=Fraction(1<<c,z)
    assert ordinary.denominator==twisted.denominator==1
    merged=(ordinary+twisted)/2;assert merged.denominator==1
    return {'ordinary':int(ordinary),'self_converse':int(twisted),'merged':int(merged)}

def pairs(n):return tuple(combinations(range(n),2))
def freepairs(n):return tuple((i,j) for i in range(n) for j in range(i+2,n))
def out_from_code(n,c):
    out=[0]*n
    for k,(i,j) in enumerate(pairs(n)):
        a,b=(i,j) if c>>k&1 else (j,i);out[a]|=1<<b
    return out

def fixed_code(n,x):
    f={p:k for k,p in enumerate(freepairs(n))}
    return sum((1 if j==i+1 else x>>f[(i,j)]&1)<<k for k,(i,j) in enumerate(pairs(n)))

def relabel_code(n,out,p):
    return sum((out[p[i]]>>p[j]&1)<<k for k,(i,j) in enumerate(pairs(n)))

def reflected(n,x):
    f=freepairs(n);index={p:k for k,p in enumerate(f)}
    return sum((x>>k&1)<<index[(n-1-j,n-1-i)] for k,(i,j) in enumerate(f))

def tiling_of_path(n,out,p):
    assert all(out[p[i]]>>p[i+1]&1 for i in range(n-1))
    return sum((out[p[i]]>>p[j]&1)<<k for k,(i,j) in enumerate(freepairs(n)))

def canonical_certificate(n):
    """Generate positive controls, retaining explicit maps to representatives."""
    m=len(freepairs(n));N=1<<m;mask=(1<<len(pairs(n)))-1
    cache={};data=[];perms=list(permutations(range(n)))
    for x in range(N):
        code=fixed_code(n,x)
        if code not in cache:
            out=out_from_code(n,code)
            orbit=[(relabel_code(n,out,p),p) for p in perms]
            key,pbest,epsbest=min((c if e==0 else c^mask,p,e) for c,p in orbit for e in (0,1))
            inv=[0]*n
            for i,v in enumerate(pbest):inv[v]=i
            for c,p in orbit:
                q=[inv[v] for v in p]
                cache[c]=(key,q,epsbest)
                cache[c^mask]=(key,q,epsbest^1)
        data.append(cache[code])
    reps=sorted({t[0] for t in data});index={c:u for u,c in enumerate(reps)}
    halfedges=[{'node':index[key],'coordinate':x,'partner':x^(N-1),'blue':reflected(n,x)==x,'permutation':p,'converse':eps} for x,(key,p,eps) in enumerate(data)]
    return {'n':n,'representatives':reps,'halfedges':halfedges}

def verify_certificate(cert):
    """Verifier does not canonicalize and does not enumerate permutations."""
    n=cert['n'];m=len(freepairs(n));N=1<<m;reps=cert['representatives'];hs=cert['halfedges'];k=len(reps)
    assert k==burnside(n)['merged'], 'wrong complete node count'
    assert len(hs)==N and {h['coordinate'] for h in hs}==set(range(N)), 'coordinates not a bijection'
    assert {h['node'] for h in hs}==set(range(k)), 'empty or invalid node'
    assert all(0<=c<1<<len(pairs(n)) for c in reps)
    rout=[out_from_code(n,c) for c in reps]
    arc_checks=0
    for i,h in enumerate(hs):
        j=h['partner'];x=h['coordinate'];p=h['permutation'];e=h['converse'];u=h['node']
        assert 0<=j<N and j!=i and hs[j]['partner']==i, 'invalid line pairing'
        assert hs[j]['coordinate']==x^(N-1), 'incoherent complement'
        assert h['blue']==hs[j]['blue']==(reflected(n,x)==x), 'incorrect line color'
        assert sorted(p)==list(range(n)) and e in (0,1), 'invalid transporter'
        out=out_from_code(n,fixed_code(n,x))
        for a,b in pairs(n):
            assert (out[a]>>b&1)==((rout[u][p[a]]>>p[b]&1)^e), f'arc mismatch at halfedge {i}, pair {a,b}'
            arc_checks+=1
    return {'n':n,'halfedges':N,'nodes':k,'lines':N//2,'arc_equalities_checked':arc_checks,'status':'PASS'}

def assignment_matrices(n,fi,k):
    N=len(fi);B=[[0]*k for _ in range(k)];K=[[0]*k for _ in range(k)]
    for x,u in enumerate(fi):
        v=fi[x^(N-1)];A=B if reflected(n,x)==x else K;A[u][v]+=1
    return B,K

def score(n,code):return tuple(sorted(a.bit_count() for a in out_from_code(n,code)))
def folded(n,s):return min(tuple(sorted(s)),tuple(sorted(n-1-a for a in s)))
def triangles(n,s):return comb(n,3)-sum(comb(a,2) for a in s)

def endpoint_predictions(n,reps):
    """Independent representative HP/Aut endpoint-inventory calculation."""
    allmask=(1<<len(pairs(n)))-1;res=[]
    for code in reps:
        out=out_from_code(n,code);s=[a.bit_count() for a in out];aut=0;sc=False;counts=[Counter(),Counter()]
        for p in permutations(range(n)):
            q=relabel_code(n,out,p)
            aut+=q==code;sc|=q==(code^allmask)
            if not all(out[p[i]]>>p[i+1]&1 for i in range(n-1)):continue
            t=tiling_of_path(n,out,p);blue=reflected(n,t)==t
            target=[n-1-a for a in s];target[p[0]]+=1;target[p[-1]]-=1
            counts[not blue][(s[p[0]]-s[p[-1]]-1,folded(n,target))]+=1
        kap=1 if sc else 2
        scaled=[]
        for C in counts:
            assert all(v*kap%aut==0 for v in C.values())
            scaled.append(Counter({q:v*kap//aut for q,v in C.items()}))
        res.append(scaled)
    return res

def matrix_endpoint_rows(n,reps,B,K):
    s=[score(n,c) for c in reps];c3=[triangles(n,a) for a in s];ans=[]
    for u in range(len(reps)):
        cols=[]
        for A in (B,K):
            C=Counter()
            for v in range(len(reps)):
                if A[u][v]:C[(c3[v]-c3[u],folded(n,s[v]))]+=A[u][v]
            cols.append(C)
        ans.append(cols)
    return ans

def spectral_audit(n,fi,k):
    m=len(freepairs(n));N=len(fi);mass=[fi.count(u) for u in range(k)]
    W=[sp.zeros(k) for _ in range(m+1)]
    for x in range(N):
        for y in range(N):W[(x^y).bit_count()][fi[x],fi[y]]+=1
    Z=sp.zeros(k,N)
    for x,u in enumerate(fi):
        for s in range(N):Z[u,s]+=(-1)**((x&s).bit_count())
    for s in range(N):assert sum(Z[:,s])==(N if s==0 else 0)
    # Boolean idempotence, a full phase-compatible condition stronger than PSD.
    for u in range(k):
        for t in range(N):assert sum(Z[u,s]*Z[u,s^t] for s in range(N))==N*Z[u,t]
    Q=[];ranks=[]
    for j in range(m+1):
        q=sp.zeros(k)
        for d in range(m+1):
            Kjd=sum((-1)**a*comb(d,a)*comb(m-d,j-a) for a in range(max(0,j-m+d),min(j,d)+1))
            q+=Kjd*W[d]
        cols=[s for s in range(N) if s.bit_count()==j];z=Z[:,cols]
        assert q==z*z.T
        ranks.append(q.rank());Q.append(q)
    assert all(fi[reflected(n,x)]==fi[x] for x in range(N))
    assert all(Z[u,s]==Z[u,reflected(n,s)] for u in range(k) for s in range(N))
    f=(n-1)//2;q=(m-f)//2
    caps=[]
    for j in range(m+1):
        Fj=sum(comb(f,a)*comb(q,(j-a)//2) for a in range(f+1) if 0<=j-a<=2*q and (j-a)%2==0)
        cap=(comb(m,j)+Fj)//2;caps.append(cap);assert ranks[j]<=cap
        for u in range(k):
            target=comb(m,j) if mass[u]%2 else 4*Fj
            assert Q[j][u,u]%8==target%8
    D=sp.diag(*mass);leak=m*D+2*W[2]-W[1]*D.inv()*W[1] if m>=2 else sp.zeros(k)
    local=sp.zeros(N,k)
    for x in range(N):
        for b in range(m):local[x,fi[x^(1<<b)]]+=1
    direct=sp.zeros(k)
    for u in range(k):
        rows=[x for x in range(N) if fi[x]==u];A=local[rows,:];v=sp.ones(1,len(rows))*A
        direct+=A.T*A-v.T*v/mass[u]
    assert direct==leak
    return {'mass':mass,'top_charge':list(map(int,Z[:,N-1])),'ranks':ranks,'reflection_caps':caps,'leak_rank':leak.rank(),'boolean_convolution':'PASS','all_distance_gram':'PASS','reflection':'PASS','diagonal_congruences':'PASS','covariance_identity':'PASS'},W,Z

def compatible_node_permutations(mass):
    d=defaultdict(list)
    for u,m in enumerate(mass):d[m].append(u)
    groups=list(d.values())
    for images in product(*(list(permutations(g)) for g in groups)):
        p=list(range(len(mass)))
        for g,h in zip(groups,images):
            for a,b in zip(g,h):p[a]=b
        yield p

def hostile_audit():
    n=5;cert=canonical_certificate(n);reps=cert['representatives'];fi=[h['node'] for h in cert['halfedges']];bad=fi.copy();bad[4],bad[26]=bad[26],bad[4];k=len(reps)
    B,K=assignment_matrices(n,fi,k);Bb,Kb=assignment_matrices(n,bad,k)
    good,_,Z=spectral_audit(n,fi,k);evil,_,Zb=spectral_audit(n,bad,k)
    assert good['mass']==evil['mass'] and good['top_charge']==evil['top_charge']
    assert K==Kb and [sum(r) for r in B]==[sum(r) for r in Bb]
    assert [B[u][u] for u in range(k)]==[Bb[u][u] for u in range(k)]
    expect=endpoint_predictions(n,reps)
    assert matrix_endpoint_rows(n,reps,B,K)==expect==matrix_endpoint_rows(n,reps,Bb,Kb)
    changes=[{'codes':[reps[u],reps[v]],'old':B[u][v],'new':Bb[u][v]} for u in range(k) for v in range(u,k) if B[u][v]!=Bb[u][v]]
    isos=[];tried=0
    for p in compatible_node_permutations(good['mass']):
        tried+=1
        if all(B[u][v]==Bb[p[u]][p[v]] and K[u][v]==Kb[p[u]][p[v]] for u in range(k) for v in range(k)):isos.append(p)
    assert not isos
    # An explicit path-renormalization relation severed by the hostile.
    out=out_from_code(n,fixed_code(n,4));broken=None
    for p in permutations(range(n)):
        if all(out[p[i]]>>p[i+1]&1 for i in range(n-1)):
            y=tiling_of_path(n,out,p)
            if fi[y]==fi[4] and bad[y]!=bad[4]:broken={'source':4,'path':p,'normalized':y,'true_code':reps[fi[4]],'hostile_source_code':reps[bad[4]],'hostile_target_code':reps[bad[y]]};break
    assert broken
    # Counterexample certificate fails with unchanged witnessed transporters.
    altered=json.loads(json.dumps(cert));altered['halfedges'][4]['node'],altered['halfedges'][26]['node']=bad[4],bad[26]
    try:verify_certificate(altered);raise RuntimeError('hostile accepted')
    except AssertionError as e:reject=str(e)
    # Exhaustive direct non-isomorphism for the swapped representatives.
    out4=out_from_code(n,fixed_code(n,4));target=reps[bad[4]];mask=(1<<len(pairs(n)))-1
    no_transport=all(relabel_code(n,out4,p) not in (target,target^mask) for p in permutations(range(n)))
    assert no_transport
    return {'n':n,'swapped_tilings':[4,26],'codes':reps,'changes_blue_lines':changes,'good_cube':good,'hostile_cube':evil,'endpoint_and_score_shell_rows':'PASS unchanged','black_matrix':'PASS unchanged','blue_degree_and_loop_rows':'PASS unchanged','abstract_graph_nonisomorphic':True,'mass_preserving_node_permutations_tested':tried,'broken_path_transport':broken,'certificate_rejection':reject,'no_alternative_source_transporter':True,'old_B':B,'new_B':Bb,'K':K}

def brute_burnside_check(n):
    """Check cycle-type formula via direct signed action on unordered arc pairs."""
    ps=pairs(n);idx={e:i for i,e in enumerate(ps)};total=[0,0];compatible=[0,0]
    for p in permutations(range(n)):
        for twist in (0,1):
            seen=set();cycles=0;valid=True
            for e in range(len(ps)):
                if e in seen:continue
                v=e;par=0;cycles+=1
                while v not in seen:
                    seen.add(v);a,b=ps[v];a,b=p[a],p[b];par^=(a>b)^twist;v=idx[tuple(sorted((a,b)))]
                if par:valid=False
            if valid:compatible[twist]+=1;total[twist]+=1<<cycles
    b=burnside(n)
    assert total[0]==factorial(n)*b['ordinary'] and total[1]==factorial(n)*b['self_converse']
    return {'n':n,'compatible_permutations':compatible,'fixed_tournament_sums':total,**b}

def run():
    controls=[]
    for n in (3,4,5,6):
        cert=canonical_certificate(n);r=verify_certificate(cert);r['burnside']=brute_burnside_check(n);controls.append(r)
        if n<=5:(HERE/f'certificate_n{n}.json').write_text(json.dumps(cert,indent=2)+'\n')
        print('certificate',r,flush=True)
    hostile=hostile_audit();print('hostile', {k:v for k,v in hostile.items() if k not in ('old_B','new_B','K','good_cube','hostile_cube')},flush=True)
    report={'positive_controls':controls,'burnside_counts_n1_to12':[{'n':n,**burnside(n)} for n in range(1,13)],'hostile':hostile}
    (HERE/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print('ALL EXACT CHECKS PASS',flush=True)

if __name__=='__main__':run()
