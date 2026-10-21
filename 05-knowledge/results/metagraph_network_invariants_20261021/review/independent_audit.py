#!/usr/bin/env python3
"""Independent exact audit. Standard library only; no producer imports.
Run anywhere: python3 independent_audit.py --repo PATH --producer PATH.
Inputs are content-pinned before use. Assertions are deliberately not used.
"""
import argparse, hashlib, itertools, json, math
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser()
ap.add_argument('--repo',type=Path,default=HERE.parent/'math')
ap.add_argument('--producer',type=Path,default=HERE.parent/'session20-metagraph-networks')
args=ap.parse_args()
SRC=args.repo/'05-knowledge/results/merged_metagraph_black_quotient_bridges_20261020'
pins=json.loads((HERE/'frozen_inputs.json').read_text())
checks=0

def check(ok,why):
    global checks
    checks+=1
    if not ok:raise RuntimeError(why)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for name,h in pins['producer'].items():check(sha(args.producer/name)==h,'Changed producer '+name)
for name,h in pins['atlas'].items():check(sha(SRC/name)==h,'Changed frozen atlas '+name)
for name,h in json.loads((HERE/'inheritance_pins.json').read_text()).items():check(sha(args.repo/name)==h,'Changed inherited proof '+name)


def determinant(a):
    """Integer Bareiss elimination, independent of SymPy determinant."""
    a=[list(x) for x in a];n=len(a)
    if not n:return 1
    sign=1;old=1
    for k in range(n-1):
        p=next((i for i in range(k,n) if a[i][k]),None)
        if p is None:return 0
        if p!=k:a[k],a[p]=a[p],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                num=a[i][j]*pivot-a[i][k]*a[k][j]
                check(num%old==0,'Bareiss division')
                a[i][j]=num//old
            a[i][k]=0
        old=pivot
    return sign*a[-1][-1]


def inverse(a):
    n=len(a);w=[[F(v) for v in row]+[F(i==j) for j in range(n)] for i,row in enumerate(a)]
    for k in range(n):
        p=next((i for i in range(k,n) if w[i][k]),None)
        if p is None:raise RuntimeError('Singular reduced Laplacian')
        w[k],w[p]=w[p],w[k];v=w[k][k];w[k]=[x/v for x in w[k]]
        for i in range(n):
            if i!=k:
                v=w[i][k]
                if v:w[i]=[x-v*y for x,y in zip(w[i],w[k])]
    return [row[n:] for row in w]


def rank2(rows):
    # Column-by-column row elimination rather than producer's bit-basis insertion.
    rows=list(rows);rank=0;bits=max(rows,default=0).bit_length()
    for j in range(bits):
        p=next((k for k in range(rank,len(rows)) if (rows[k]>>j)&1),None)
        if p is None:continue
        rows[rank],rows[p]=rows[p],rows[rank]
        for k in range(len(rows)):
            if k!=rank and (rows[k]>>j)&1:rows[k]^=rows[rank]
        rank+=1
    return rank


def components(n,edges):
    p=list(range(n))
    def root(x):
        while p[x]!=x:x=p[x]
        return x
    for (u,v),w in edges.items():
        if w:p[root(u)]=root(v)
    d={}
    for i in range(n):d.setdefault(root(i),[]).append(i)
    return list(d.values())


def lap(n,edges):
    a=[[0]*n for _ in range(n)]
    for (u,v),w in edges.items():
        check(w>=0,'Negative edge multiplicity')
        if u!=v:
            a[u][u]+=w;a[v][v]+=w;a[u][v]-=w;a[v][u]-=w
    return a


def reduced(a):return [r[1:] for r in a[1:]]
def tau(n,e):return determinant(reduced(lap(n,e)))
def deg(n,e):
    d=[0]*n
    for (u,v),w in e.items():d[u]+=w;d[v]+=w
    return d


def tree_polynomial(n,B,Q):
    """Determinant polynomial by exact finite differences at 0,...,n-1."""
    vals=[tau(n,B+Counter({e:k*w for e,w in Q.items()})) for k in range(n)]
    coef=[F(0)]*n;basis=[F(1)]
    for k in range(n):
        for j,x in enumerate(basis):coef[j]+=vals[0]*x
        vals=[b-a for a,b in zip(vals,vals[1:])]
        nxt=[F(0)]*(len(basis)+1)
        for j,x in enumerate(basis):nxt[j]-=k*x/(k+1);nxt[j+1]+=x/(k+1)
        basis=nxt
    check(all(x.denominator==1 for x in coef),'Integral tree coefficients')
    return {k:int(v) for k,v in enumerate(coef) if v}


def contract(n,e,parts):
    owner={v:i for i,c in enumerate(parts) for v in c};out=Counter()
    for (u,v),w in e.items():
        a,b=sorted((owner[u],owner[v]))
        if a!=b:out[a,b]+=w
    return len(parts),out


def endpoint(n,internal,cross):
    cs=components(n,internal);out=tau(*contract(n,cross,cs))
    for c in cs:
        ix={v:i for i,v in enumerate(c)}
        e=Counter({tuple(sorted((ix[u],ix[v]))):w for (u,v),w in internal.items() if u in ix and v in ix})
        out*=tau(len(c),e)
    return out


def trim(p):
    p=list(p)
    while p and p[-1]==0:p.pop()
    return p


def polyrem(a,b):
    a=trim(a);b=trim(b)
    if not b:raise RuntimeError('Polynomial zero divisor')
    while len(a)>=len(b):
        v=a[-1]/b[-1];k=len(a)-len(b)
        for j,x in enumerate(b):a[k+j]-=v*x
        a=trim(a)
    return a


def all_roots_negative(p):
    # Sturm on the square-free part; repeated roots retain the same location.
    f=[F(x) for x in p];g=[i*f[i] for i in range(1,len(f))]
    if not g:return True
    a,b=f,g
    while b:a,b=b,polyrem(a,b)
    gcd=[x/a[-1] for x in a]
    rem=list(f);quot=[F(0)]*(len(f)-len(gcd)+1)
    while len(rem)>=len(gcd):
        v=rem[-1]/gcd[-1];k=len(rem)-len(gcd);quot[k]=v
        for j,x in enumerate(gcd):rem[k+j]-=v*x
        rem=trim(rem)
    check(not rem,'Exact square-free quotient')
    seq=[trim(quot),[i*quot[i] for i in range(1,len(quot))]]
    while seq[-1]:
        r=polyrem(seq[-2],seq[-1])
        if not r:break
        seq.append([-x for x in r])
    def changes(values):
        signs=[1 if x>0 else -1 for x in values if x]
        return sum(a!=b for a,b in zip(signs,signs[1:]))
    minus_infty=changes([p[-1]*(-1)**(len(p)-1) for p in seq])
    zero=changes([p[0] for p in seq])
    return minus_infty-zero==len(quot)-1


def validate_tree(n,B,Q,expected=None):
    p=tree_polynomial(n,B,Q);a=len(components(n,B))-1;b=n-len(components(n,Q));r=b-a
    check(set(p)==set(range(a,b+1)),'Full coefficient interval')
    check(all(v>0 for v in p.values()),'Positive coefficients')
    check(all_roots_negative([p[k] for k in range(a,b+1)]),'Exact negative-root Sturm certificate')
    check(p[a]==endpoint(n,B,Q),'Low coefficient contraction')
    check(p[b]==endpoint(n,Q,B),'High coefficient contraction')
    check(sum(v*2**k for k,v in p.items())==tau(n,B+Counter({e:2*w for e,w in Q.items()})),'Literal factor two')
    for j in range(1,r):
        check(F(p[a+j],math.comb(r,j))**2 >= F(p[a+j-1],math.comb(r,j-1))*F(p[a+j+1],math.comb(r,j+1)),'Newton inequality')
    if expected is not None:check(p=={int(k):v for k,v in expected.items()},'Stored tree polynomial')
    return p


def canonical_full(n,tiles,t):
    """Classify full adjacency bitstrings under ALL permutations, no HP filtering."""
    pairs=list(itertools.combinations(range(n),2));arcs={(n-i,n-i-1) for i in range(1,n)}
    # Coordinates are the atlas labels 1,...,n, with base path n -> ... -> 1.
    arcs={(i,i-1) for i in range(2,n+1)}
    for k,(x,y) in enumerate(tiles):arcs.add((y,x) if (t>>k)&1 else (x,y))
    def code(p):return sum(int((p[i],p[j]) in arcs)<<k for k,(i,j) in enumerate(pairs))
    codes=[code(p) for p in itertools.permutations(range(1,n+1))]
    return min(codes),min(((1<<len(pairs))-1)^x for x in codes)


records=[]
producer=json.loads((args.producer/'probe.json').read_text())['records']
actual={}
for n in range(3,8):
    d=json.loads((SRC/f'census_n{n}.json').read_text());N=len(d['nodes']);tiles=[tuple(x) for x in d['tiles']];m=len(tiles);allbits=(1<<m)-1
    check(len(d['pi'])==1<<m,'Full tiling cube')
    ix={xy:i for i,xy in enumerate(tiles)};reflection=[ix[n+1-y,n+1-x] for x,y in tiles]
    check(reflection==d['reflection_permutation'],'Reflection coordinate list')
    def reflect(t):return sum(((t>>i)&1)<<j for i,j in enumerate(reflection))
    check(any(i==j for i,j in enumerate(reflection)),'Reflection fixed coordinate')
    lines={min(t,t^allbits) for t in range(1<<m)};B=Counter();K=Counter();Q=Counter();black_orbits=set()
    for t in lines:
        rt=reflect(t);rr=min(rt,rt^allbits);edge=tuple(sorted((d['pi'][t],d['pi'][t^allbits])))
        check(d['pi'][rt]==d['pi'][t],'Owner reflection equivariance')
        if rt==t:B[edge]+=1
        else:
            check(rr!=t,'Black reflection orbit free')
            K[edge]+=1;black_orbits.add(tuple(sorted((t,rr))))
    for t,s in black_orbits:
        edge=tuple(sorted((d['pi'][t],d['pi'][t^allbits])))
        check(edge==tuple(sorted((d['pi'][s],d['pi'][s^allbits]))),'Paired endpoints')
        Q[edge]+=1
    for name,e in [('blue_edges',B),('black_edges',K),('paired_black_edges',Q)]:
        check(e==Counter({(u,v):w for u,v,w in d[name]}),'Rebuilt '+name)
    check(K==Counter({e:2*w for e,w in Q.items()}),'Literal/paired multiplicity')
    G=B+K;check(len(components(N,G))==1,'Finite actual connectivity')
    masses=Counter(d['pi']);bd=deg(N,B);qd=deg(N,Q)
    for u,node in enumerate(d['nodes']):
        check(masses[u]==node['mass']==bd[u]+2*qd[u],'Stub mass, including loops')
        check(bd[u]==node['blue'],'Blue stub count')
        check(qd[u]==node['reflection_pairs'],'Paired degree')
        check(qd[u]%2==node['quotient_black_charge'],'Charge')
        check((bd[u]%2==1) if node['SC'] else (bd[u]==0 and qd[u]%2==1),'SC/NS parity')
    if n==5:
        regular=d['pi'][11];neighbors={v if u==regular else u for (u,v),w in B.items() if w and regular in (u,v)}
        check(masses[regular]==bd[regular]==3 and qd[regular]==0 and len(neighbors)==3,'Actual pure-blue nonleaf hostile')
    if n<=5:
        cls={};merged={};sc={};canon=[]
        for t in range(1<<m):
            c,co=canonical_full(n,tiles,t);canon.append(c)
            cls.setdefault(c,[]).append(t);merged.setdefault(min(c,co),[]).append(t);sc[t]=c==co
        expected_classes={frozenset(t for t,c in enumerate(d['class_by_mask']) if c==j) for j in range(len(d['classes']))}
        expected_merged={frozenset(t for t,u in enumerate(d['pi']) if u==j) for j in range(N)}
        check({frozenset(v) for v in cls.values()}==expected_classes,'Independent full-adjacency ordinary classification')
        check({frozenset(v) for v in merged.values()}==expected_merged,'Independent converse-merged classification')
        check(all(sc[t]==d['nodes'][d['pi'][t]]['SC'] for t in range(1<<m)),'Independent SC labels')
    LB=lap(N,B);LG=lap(N,G);cb=len(components(N,B));lb_rank=rank2([sum((v%2)<<i for i,v in enumerate(row)) for row in LB])
    # Literal blue incidence, with a column for every parallel edge and zero loop columns.
    incident=[0]*N;col=0
    for (u,v),w in B.items():
        for _ in range(w):
            if u!=v:incident[u]^=1<<col;incident[v]^=1<<col
            col+=1
    incidence_rank=rank2(incident);bicycle=incidence_rank-lb_rank
    rho=N-1-rank2([sum((v%2)<<i for i,v in enumerate(row)) for row in reduced(LG)])
    check(incidence_rank==N-cb,'Blue incidence rank')
    check(rho==cb-1+bicycle,'Critical group/bicycle rank')
    saved=producer[n-3]
    for key,value in [('vertices',N),('NS',sum(not x['SC'] for x in d['nodes'])),('blue_components',cb),('blue_bicycle_dimension',bicycle),('critical_group_2_rank',rho)]:check(saved[key]==value,'Stored '+key)
    rec={'n':n,'classification':'independent full adjacency permutations' if n<=5 else 'inherited frozen owner labels','vertices':N,'blue_components':cb,'Q_components':len(components(N,Q)),'bicycle':bicycle,'two_rank':rho,'literal_edges':sum(G.values()),'loops':sum(w for (u,v),w in G.items() if u==v)}
    if n<=5:rec['tree_polynomial']=validate_tree(N,B,Q,saved['tree_polynomial_B_plus_xQ'])
    if n<=6:
        t=tau(N,G);check(t==saved['spanning_trees'],'Stored tree count');rec['trees']=t
        inv=inverse(reduced(LG));exp=math.lcm(*(x.denominator for row in inv for x in row));check(exp==max(saved['smith_invariants']),'Independent group exponent from inverse')
        # The claimed Smith lists are determined by order, exponent, and 2-rank here.
        check(math.prod(saved['smith_invariants'])==t,'Smith order')
        check(all(x in (1,2,exp) for x in saved['smith_invariants']),'Claimed Smith shape')
        check(t//exp==2**max(0,rho-1),'Independent Smith shape: all remaining nonunits are two')
        cols=[[F(0)]*N]+[[F(0)]+[inv[i][u-1] for i in range(N-1)] for u in range(1,N)]
        count=0;oddres=0
        for u in range(N):
            if d['nodes'][u]['SC']:continue
            for v in range(N):
                if u==v:continue
                phi=[cols[u][i]-cols[v][i] for i in range(N)]
                check(all(sum(LG[i][j]*phi[j] for j in range(N))==int(i==u)-int(i==v) for i in range(N)),'Full Kirchhoff equations')
                den=math.lcm(*(x.denominator for x in phi));check(den%2==0,'Even full voltage denominator')
                count+=1;oddres+=((phi[u]-phi[v]).denominator%2==1)
                if n==5 and u==2 and v==0:
                    check(phi[2]==F(6092,5137) and phi[1]==F(11327,10274) and den==10274,'Exact Kron witness')
                    check(phi[1]/phi[2]==F(11327,12184),'Exact harmonic extension')
                    # Unit boundary source voltage: all interior equations vanish; net current 1/R.
                    boundary=[x/phi[2] for x in phi]
                    currents=[sum(LG[i][j]*boundary[j] for j in range(N)) for i in range(N)]
                    check(currents[2]==F(5137,6092) and currents[0]==-F(5137,6092) and all(currents[i]==0 for i in range(N) if i not in (0,2)),'Exact two-terminal response')
        rec.update({'group_exponent':exp,'NS_pairs_checked':count,'odd_resistance_denominator_pairs':oddres})
    actual[n]=(d,B,Q);records.append(rec)

# Degree/colour/loop-preserving abstract switch. NOT a realizability assertion.
d,B,Q=actual[5];N=len(d['nodes']);q=Q.copy()
for e in [(1,4),(2,3)]:check(q[e]>=1,'Switch edge exists');q[e]-=1
for e in [(1,3),(2,4)]:q[e]+=1
check(all(w>=0 for w in q.values()),'Switch nonnegative')
check(deg(N,Q)==deg(N,q),'Switch degrees')
check({e:w for e,w in Q.items() if e[0]==e[1]}=={e:w for e,w in q.items() if e[0]==e[1]},'Switch loops')
check(sum(Q.values())==sum(q.values()),'Switch edge count')
check(tau(N,B+Counter({e:2*w for e,w in q.items()}))==39232,'Switch actual physical tree count')
check(tau(N,B+Counter({e:2*w for e,w in Q.items()}))==41096,'Actual physical tree count')

# Standard finite multigraph hostiles, with disconnected graphs treated COMPONENTWISE.
generic=0
for n in range(1,5):
    pairs=list(itertools.combinations(range(n),2))
    for colors in itertools.product(range(3),repeat=len(pairs)):
        B=Counter({e:1 for e,c in zip(pairs,colors) if c==1});Q=Counter({e:1 for e,c in zip(pairs,colors) if c==2})
        for comp in components(n,B+Q):
            ix={v:i for i,v in enumerate(comp)}
            b=Counter({tuple(sorted((ix[u],ix[v]))):w for (u,v),w in B.items() if u in ix and v in ix})
            q=Counter({tuple(sorted((ix[u],ix[v]))):w for (u,v),w in Q.items() if u in ix and v in ix})
            validate_tree(len(comp),b,q);generic+=1
# Blue parallel edges are a bicycle; paired black alone guarantees two-torsion.
validate_tree(2,Counter({(0,1):2,(0,0):7}),Counter({(0,1):3,(1,1):4}))
check(tau(2,Counter({(0,1):2,(0,0):7,(1,1):4}))==2,'Loops electrically invisible')
check(rank2([0,0])==0,'Double blue edge gives one bicycle and one two-rank')
# One Q bridge lifts to two nonbridge literal black lines.
check(tau(2,Counter({(0,1):2}))==2 and tau(2,Counter({(0,1):1}))==1,'Q bridge versus literal nonbridge')
# Loop holding convention: two stubs per loop changes transition clock, never the DC L.
check(deg(2,Counter({(0,1):1,(0,0):3}))==[7,1],'Loop holding stubs')

result={'status':'PASS','explicit_checks':checks,'generic_component_cases':generic,'atlas_pin_enforced':True,'producer_pin_enforced':True,'independent_methods':['full-adjacency permutation classifier through n=5','C/R line and pair reconstruction through n=7','integer Bareiss determinants','exact finite-difference tree polynomial and Sturm negative-root check','full literal blue incidence rank','Fraction Gauss-Jordan voltage solves through n=6','group exponent from inverse denominators'],'records':records,'scope':'No expanded atlas census, no all-order connectivity, no alternative owner-map realizability, no reserved theorem inputs.'}
print(json.dumps(result,sort_keys=True,indent=2))
(HERE/'independent_results.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
