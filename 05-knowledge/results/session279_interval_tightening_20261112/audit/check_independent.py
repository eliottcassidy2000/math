"""Independent small controls for Session279; no inherited partition generator."""
from pathlib import Path
import importlib.util,sys,json
from fractions import Fraction as F
from math import gcd,factorial,comb
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'enumeration'))
import backbone
import second_jet

def occupancies(n,a=1):
    if a>n:
        if not n:yield {}
        return
    for j in range(n//a+1):
        for d in occupancies(n-j*a,a+2):
            if j:d={a:j,**d}
            yield d

def brute(n,B):
    u=c=c2=F(0)
    for d in occupancies(n):
        exponent=sum(j*(a-1)//2+a*comb(j,2) for a,j in d.items())
        moment=square=F(0)
        ls=sorted(d)
        for i,a in enumerate(ls):
            for b in ls[i+1:]:
                base=gcd(a,b) if a in B or b in B else a
                exponent+=d[a]*d[b]*base
                q=1-F(1,1<<(base-gcd(a,b)))
                moment+=d[a]*d[b]*q
                square+=d[a]*d[b]*q*q
        den=1
        for a,j in d.items():den*=a**j*factorial(j)
        w=F(2**exponent,den)
        u+=w;c+=w*moment;c2+=w*(moment*moment-square)/2
    return u,c,c2

rows=[]
for n in range(25):
    for B in [(),(3,),(5,),(3,7),(5,9),(3,5,7,9),(31,)]:
        u,c,_=backbone.envelope(n,B,True)
        assert (u,c)==brute(n,B)[:2],(n,B)
        rows.append((n,B))
jet_rows=0
for n in range(21):
    for B in [(),(3,5),(5,9)]:
        cs,_=second_jet.second_jet(n,B)
        assert tuple(cs)==brute(n,B),(n,B)
        jet_rows+=1
report={'status':'PASS','independent_first_moment_controls':len(rows),'independent_second_moment_controls':jet_rows,'first_moment_orders':[0,24],'second_moment_orders':[0,20]}
# Marked-graph compiler compared with separate occupancy/edge-polynomial summation.
sys.path.insert(0,str(P/'marked'))
import marked_graph

def direct_coefficient(n,B,k):
    total=F(0)
    for d in occupancies(n):
        e=sum(j*(a-1)//2+a*comb(j,2) for a,j in d.items())
        coeff=[F(1)]+[F(0)]*k
        ls=sorted(d)
        for i,a in enumerate(ls):
            for b in ls[i+1:]:
                base=gcd(a,b) if a in B or b in B else a
                mult=d[a]*d[b];e+=mult*base
                q=1-F(1,1<<(base-gcd(a,b)))
                coeff=[sum((coeff[h]*comb(mult,j-h)*q**(j-h) for h in range(j+1) if j-h<=mult),F(0)) for j in range(k+1)]
        z=1
        for a,j in d.items():z*=a**j*factorial(j)
        total+=F(2**e,z)*coeff[k]
    return total

marked_checks=0
for n,B,k in [(n,B,k) for n in range(16) for B in [(),(3,5)] for k in range(3)]+[(14,(),3),(18,(),3),(22,(3,),3)]:
    actual,_=marked_graph.marked_coefficient(n,k,B)
    assert actual==direct_coefficient(n,B,k),(n,B,k)
    marked_checks+=1
report['independent_marked_graph_controls']=marked_checks
import integer_jet
int_controls=0
for n in range(25):
    for B in [(),(3,5),(5,9)]:
        u,c,_=integer_jet.integer_jet(n,B)
        assert (u,c)==brute(n,B)[:2]
        int_controls+=1
report['independent_integer_jet_controls']=int_controls

sys.path.insert(0,str(P/'cross'))
import ground_moments as gm
from itertools import combinations,permutations,product
# Independent root polynomials over all possible four-contact vertices.
vertices=0
for m in range(5,15):
    ls,us=gm.duals(m)
    for roots in combinations(range(1,m+1),4):
        vals=[]
        for x in range(m+1):
            val=F(1)
            for r in roots:val*=F(r-x,r)
            vals.append(val)
        upper=all(v>=0 for v in vals[1:])
        lower=all(v<=0 for v in vals[1:])
        assert upper==(roots in us)
        assert lower==(roots in ls)
        vertices+=1
report['independent_dual_vertices_classified']=vertices
# Arbitrary positive mixtures, unrelated to tournament partitions.
mixture_controls=0
for m in range(5,15):
    for t in range(1,8):
        dist=[F(((x+1)*(t+3))%11+1,17) for x in range(m+1)]
        cs=[sum((v*comb(x,j) for x,v in enumerate(dist) if x>=j),F(0)) for j in range(5)]
        lo,hi,_,_=gm.moment_interval(m,cs)
        assert lo<=dist[0]<=hi
        mixture_controls+=1
report['arbitrary_positive_moment_mixtures']=mixture_controls

# Exact independent matrix/exterior calculation and direct permutation packets.
def prodmat(a,b):
    return tuple(tuple(sum(x*y for x,y in zip(row,col)) for col in zip(*b)) for row in a)
M=((3,3,-1),(3,4,0),(1,0,0)); pairs=((0,1),(0,2),(1,2))
T=tuple(tuple(M[i][a]*M[j][b]-M[i][b]*M[j][a] for a,b in pairs) for i,j in pairs)
T4=prodmat(prodmat(T,T),prodmat(T,T))
assert max(sum(abs(x) for x in row) for row in T4)==992
rays=[(1,1,0),(1,1,1),(0,1,0)];Us={}
def act(v):return tuple(sum(x*y for x,y in zip(row,v)) for row in M)
for k in range(1,29):
    rays=list(map(act,rays))
    if k==3:
        growth=[F(act(v)[0],v[0]) for v in rays]
        g,G=min(growth),max(growth)
        assert (g,G)==(F(97,15),F(797,123))
    if k>=25:
        E=max(abs(a[i]*b[j]-a[j]*b[i]) for a,b in combinations(rays,2) for i,j in pairs)
        U=F(96*max(v[0] for v in rays)*E,min(v[0] for v in rays)**2)
        assert U<1;Us[k]=str(U)
assert 992*(G/g**2)**4<1
report['independent_uniform_threshold_bounds']=Us

spec=importlib.util.spec_from_file_location('physical_verifier',P/'physical'/'verify.py')
pv=importlib.util.module_from_spec(spec);spec.loader.exec_module(pv)
from functools import lru_cache

def perm_packet(words):
    # Build graph separately with an adjacency relation on (block,within-block).
    vertices=[(b,j) for b,w in enumerate(words) for j in range(len(w)+1)]
    def edge(u,v):
        b,i=u;c,j=v
        if b!=c:
            return b>c if i==j==0 and abs(b-c)==1 else b<c
        if i==0:return words[b][j-1]=='1'
        if j==0:return words[b][i-1]=='0'
        return i<j
    @lru_cache(None)
    def h(sub):
        return sum(all(edge(p[i],p[i+1]) for i in range(len(p)-1)) for p in permutations(sub))
    root=(len(words)-1,0);o=e=0
    others=tuple(v for v in vertices if v!=root)
    for length in range(len(vertices)):
        for tail in permutations(others,length):
            path=(root,)+tail
            if not all(edge(path[i],path[i+1]) for i in range(len(path)-1)):continue
            remainder=tuple(v for v in vertices if v not in path)
            if len(path)%2:o+=h(remainder)
            else:e+=h(remainder)
    return h(tuple(vertices)),o+e,o-e
words=[''.join(w) for l in range(1,3) for w in product('01',repeat=l)]
physical_controls=0
for ws in [(w,) for w in words]+list(product(words,repeat=2))+[('0','10'),('10','10')]:
    assert perm_packet(ws)==pv.packet(ws),(ws,perm_packet(ws),pv.packet(ws))
    physical_controls+=1
report['independent_permutation_physical_packets']=physical_controls


import divisor_jet
for n in range(25):
    for B in [(),(3,5),(5,9)]:
        u,c,_=divisor_jet.divisor_jet(n,B)
        assert (u,c)==brute(n,B)[:2]
report['independent_divisor_jet_controls']=75

# Cross-representation exact integer at n64, not a timing comparison.
mark64=json.loads((P/'marked'/'order64.json').read_text())
cert=gm.certify(json.loads((P/'cross'/'packet_n64_B11.json').read_text()))['ground_with_R']
assert cert['low']==cert['high']==mark64['exact_count']
report['order64_independent_representations_agree_bits']=int(cert['low']).bit_length()

# Characteristic coefficients by minors, not by the producer's scalar formula.
criterion_rows=[]
for size in range(2,13):
    failures=[]
    for middle in product('01',repeat=size-2):
        word='1'+''.join(middle)+'0';A=pv.block(word)
        tr=sum(A[i][i] for i in range(3))
        minors=sum(A[i][i]*A[j][j]-A[i][j]*A[j][i] for i,j in pairs)
        det=sum((1 if perm in [(0,1,2),(1,2,0),(2,0,1)] else -1)*A[0][perm[0]]*A[1][perm[1]]*A[2][perm[2]] for perm in permutations(range(3)))
        assert det==2**size
        f1=1-tr+minors-det;fm1=-1-tr-minors-det
        if not (f1<0 and fm1<0):failures.append(word)
    criterion_rows.append({'length':size,'failures':len(failures)})
    if size<8:assert failures==[]
    if size==8:assert failures==['10101010']
report['independent_endpoint_classification']=criterion_rows


def core_parts(n,lo=13):
    if n==0:yield ()
    for a in range(lo,n+1,2):
        for tail in core_parts(n-a,a):yield (a,)+tail
bad=[]
for size in range(71):
    for core in core_parts(size):
        defects=[min(a,b)-gcd(a,b) for a,b in combinations(core,2) if min(a,b)>gcd(a,b)]
        if defects and sum(d>min(defects) for d in defects)>4:bad.append((size,core))
assert bad and min(x[0] for x in bad)==70
assert (70,(13,15,17,25)) in bad
report['ground_four_moment_first_core_rank5']=bad
print(json.dumps(report,indent=2))
