"""Independent exact referee: complete individual-danger wall partition.

No producer imports. Integer product atlas, complement-mask full profiles,
center-pair geometry, and a larger 548-wall partition independently certify
the actual entry and uniform correlated count. Output is actual raw LF.
"""
from fractions import Fraction as F
from math import gcd,prod,isqrt
from functools import reduce
from itertools import combinations
from pathlib import Path
import json,hashlib,argparse,sys
sys.stdout.reconfigure(newline='\n')
N=0
def need(ok,label):
    global N
    if not ok:raise RuntimeError(label)
    N+=1
def allowed_sums():
    out={1}
    for p in range(2,357):
        if p%3==2 and all(p%d for d in range(2,isqrt(p)+1)):
            for a in tuple(out):
                for b in (p,p*p):
                    if a*b<=356:out.add(a*b)
    return out-{1}
def graph(v,allowed):return [(i,j) for i,j in combinations(range(len(v)),2) if (v[i]+v[j])//gcd(v[i],v[j]) in allowed]
def components(n,ee):
    parent=list(range(n))
    def find(i):
        while parent[i]!=i:i=parent[i]
        return i
    for i,j in ee:parent[find(j)]=find(i)
    return sorted(sorted(i for i in range(n) if find(i)==r) for r in {find(i) for i in range(n)})
def rank(a):
    a=[row[:] for row in a];r=0
    for j in range(len(a[0])):
        k=next((k for k in range(r,len(a)) if a[k][j]),None)
        if k is None:continue
        a[r],a[k]=a[k],a[r]
        for k in range(r+1,len(a)):
            if a[k][j]:
                p,q=a[r][j],a[k][j]
                a[k]=[p*u-q*v for u,v in zip(a[k],a[r])]
                d=reduce(gcd,map(abs,a[k]))
                if d:a[k]=[v//d for v in a[k]]
        r+=1
        if r==len(a):break
    return r
def packet(row,d):
    safe=(1<<d)-1
    for v in row:
        bad=0;r=0
        for j in range(d):
            if 14*min(r,d-r)<d:bad|=1<<j
            r=(r+v)%d
        safe &= ~bad
    return tuple(j for j in range(d) if safe>>j&1)
def clear(row,a,b):return F(min(min(v*a%b,(-v*a)%b) for v in row),b)
def lengths(p,q):
    p,q=sorted((p,q));out=[];den=p*q
    for i in range(p):
        for j in range(q):
            k=(q*i-p*j)%den;delta=min(k,den-k)
            if 14*delta<p+q:out.append(F(min(2*p,p+q-14*delta),14*den))
    return sorted(out)
def ceil(x):return -(-x.numerator//x.denominator)
def literal_count(p,q,t,alpha):
    a,b=alpha.numerator,alpha.denominator;den=b*t
    count=0
    for j in range(t):
        num=a+b*j
        count+=14*min(p*num%den,(-p*num)%den)<den and 14*min(q*num%den,(-q*num)%den)<den
    return count

here=Path(__file__).resolve().parent
default=here/'overnight12_20260906_lrc_decoder_descent_inherited_profiles.json'
if not default.exists():default=Path('C:/w/s0905/04-computation')/default.name
ap=argparse.ArgumentParser();ap.add_argument('--profiles',type=Path,default=default)
raw=ap.parse_args().profiles.read_bytes()
need(hashlib.sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','full raw profile pin')
profiles={int(k):{(d,tuple(w)) for d,w in v['profiles']} for k,v in json.loads(raw)['levels'].items()}
Q=91**6;leaves=(198,215,251,257,263);groups=(107,149,167,173,179,191,197)
P=prod(leaves);H=prod(groups);V=tuple(sorted([P]+[83*P//v for v in leaves]));U=tuple(sorted(H//v for v in groups));t=968
need(P==722214566370 and H==3102310371501629,'literal cofactor products')
for xs in (leaves,groups):
    for a,b in combinations(xs,2):need(gcd(a,b)==1,'pairwise cofactor coprimality')
need(reduce(gcd,V)==reduce(gcd,U)==1 and gcd(prod(V),H)==1,'primitive disjoint component supports')
allowed=allowed_sums();ve,ue=graph(V,allowed),graph(U,allowed)
need(len(ve)==5 and components(6,ve)==[list(range(6))],'complete actual six-tree')
need(len(ue)==6 and components(7,ue)==[list(range(7))],'complete actual seven-tree')
need(sorted({(U[i]+U[j])//gcd(U[i],U[j]) for i,j in ue})==[274,298,340,346],'literal accepted seven-tree sums')
dmin=min(gcd(a,b) for a,b in combinations(V,2))
need(dmin==886860810,'minimum primitive pair gcd')
need(max(Q//dmin+1,388*Q//min(V)+1)==967,'exact sufficient class threshold')
need(967*dmin>Q and 967*min(V)>388*Q,'strict crossing inequalities for every later scale')
need(1042*sum(V)+sum(U)<Q*Q,'entire class box')
need(min(V)>F(3*Q,28) and min(U)>1 and min(V)>1 and max(U)>28,'scope of unit/minimum shortcuts')
need(clear(V,1,4)==F(1,4),'primitive quarter phase')
alllengths={p:lengths(*p) for p in combinations(groups,2)}
for (p,q),ls in alllengths.items():
    need(q>=149 and max(ls)<=F(1,7*q),'complete individual length cap')
    need(all(1042*l<1 and ceil(1042*l)-1==0 for l in ls),'all individual lower credits vanish through class endpoint')
row=tuple(t*v for v in V)+U
need(967<=t<=1042 and gcd(t,H)==1,'main scale class')
need(reduce(gcd,row)==1 and len(set(row))==13 and sum(row)==2075281984219247<Q*Q,'actual primitive distinct boxed row')
ee=graph(row,allowed)
need(len(ee)==11 and components(13,ee)==[list(range(6)),list(range(6,13))],'complete actual physical graph')
matrix=[]
for i,j in ee:
    d=gcd(row[i],row[j]);r=[0]*13;r[i],r[j]=row[j]//d,-row[i]//d
    need(sum(v*c for v,c in zip(row,r))==0 and max(map(abs,r))<=355,'literal bounded weighted edge')
    matrix.append(r)
need(rank(matrix)==11,'independent fraction-free weighted rank')
fwd=[];rev=[]
for triple in combinations(range(13),3):
    left=[i for i in triple if i<6];right=[i for i in triple if i>=6]
    if len(left)==2 and len(right)==1:
        i,j=left;k=right[0];d=gcd(row[i],row[j]);c=d//gcd(d,row[k]);need(c>Q,'forward complete support');fwd.append(c)
    elif len(left)==1 and len(right)==2:
        i,j=right;k=left[0];d=gcd(row[i],row[j]);r=F(row[k]//gcd(d,row[k]),Q*((row[i]+row[j])//d));need(r>1,'reverse complete support');rev.append(r)
need(len(fwd)==105 and len(rev)==126,'all 231 mixed supports')
need(min(fwd)==858481264080 and min(rev)==F(55157421217140,55083317447977),'crossing margins')
total=0
for r in range(1,7):
    need((1,(1,)*r) in profiles[r],'full all-one supplier')
    for outside in combinations(range(13),r):
        mask=sum(1<<i for i in outside)
        d=reduce(gcd,(row[i] for i in range(13) if not(mask>>i&1)))
        word=tuple(sorted(gcd(d,row[i]) for i in outside))
        need(d==1 and (d,word) in profiles[r],'full complement-mask profile')
        total+=1
need(total==4095,'entire retained profile universe')

# A different complete wall partition: every individual danger endpoint,
# including endpoints irrelevant to the intersection. Each is a rational
# residue t*(k +/-1/14)/speed mod1. Between these walls no indicator changes.
p,q=107,167
walls=sorted({(F(t*(14*k+sign),14*speed))%1 for speed in (p,q) for k in range(speed) for sign in (-1,1)})
need(len(walls)==548,'complete individual-danger wall universe')
midpoints=[((a+walls[(i+1)%len(walls)]+(i==len(walls)-1))/2)%1 for i,a in enumerate(walls)]
wallcounts=[literal_count(p,q,t,a) for a in walls]
cellcounts=[literal_count(p,q,t,a) for a in midpoints]
for c in wallcounts+cellcounts:need(15<=c<=23,'literal count on every wall and open cell')
need(min(wallcounts+cellcounts)==15 and max(wallcounts+cellcounts)==23,'exact all-translate minimum and maximum')
active=[i for i,c in enumerate(wallcounts) if c!=cellcounts[i] or c!=cellcounts[i-1]]
need(len(active)==78,'independent recovery of all active intersection walls')
need(literal_count(p,q,t,F(0))==15,'uniform minimum attained at zero')
ls=alllengths[(p,q)]
need(len(ls)==39 and sum(ls)==F(2553,125083) and sum(ceil(t*l)-1 for l in ls)==0,'same selected geometry, zero sum of separate minima')
need(gcd(H//(p*q),t)==1 and all(gcd(u,t)==1 for u in U),'cofactor and marginal clock bijections')
E=7*ceil(F(t,7))-t
need(E==5 and 15-E==10,'positive full-grid consumer')

for d in range(2,7):need(not packet(row,d),'first physical denominator exclusions')
need(packet(row,7)==tuple(range(1,7)) and clear(row,1,7)==F(1,7),'simpler safety explicitly retained')
safe=[]
for j in range(t):
    c=clear(row,4*j+1,4*t)
    if c>=F(1,14):
        need(c>F(1,14),'strict endpoint firewall')
        safe.append((j,c))
need(len(safe)==335 and safe[0]==(3,F(885,3872)),'literal quarter grid and first strict phase')
need(gcd(4*t,7)==1 and len(safe)>=10,'all counted weak grid phases are strict')
print('PASS actual entry class and full profiles; all separate pair-interval minima vanish.')
print('Physical968: rank11, 231 crossings,4095 profiles, sum',sum(row),'margins',min(fwd),min(rev))
print('Independent individual-danger partition:',len(walls),'walls and',len(midpoints),'cells; active intersection walls',len(active))
print('Exact located overlap minimum15, maximum23; separate minimum sum0; uniform physical safe count>=10.')
print('Literal quarter lifts335; first13/3872 clearance885/3872; first physical denominator7 with full nonzero packet.')
print('SCOPE: located counts at t968 only; bounded class entry and zero separate credits for967..1042; row already safe.')
print('PASS',N,'always-active exact gates; raw LF, no producer imports')
