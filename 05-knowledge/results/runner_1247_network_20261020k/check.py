"""Exact rotation, colouring, voltage, and electrical certificates; no LRC upgrade."""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
CHECKS=0

def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(label)

def integer(n,lo=0):
    if type(n) is not int or n<lo:raise ValueError('literal integer')
    return n

def speeds(ds):
    if type(ds) is not tuple or not ds:raise ValueError('nonempty immutable speeds')
    for d in ds:integer(d,1)
    if tuple(sorted(set(ds)))!=ds:raise ValueError('strict increasing speeds')
    return ds

def rational(q):
    if type(q) not in (int,F):raise ValueError('exact rational')
    return F(q)

def distance(t):
    t=rational(t)%1
    return min(t,1-t)

def clearance(ds,t):
    speeds(ds);rational(t)
    return min(distance(d*t) for d in ds)

def rotation(ds):
    """Complete critical-point enumeration for a finite positive integer source."""
    speeds(ds)
    denominators={2*d for d in ds}
    for a,b in combinations(ds,2):denominators.update((a+b,b-a))
    ts={F(j,q) for q in denominators for j in range(q+1)}
    vals={t:clearance(ds,t) for t in ts}; best=max(vals.values())
    return best,tuple(t for t in sorted(ts) if vals[t]==best)

def safe_intervals(ds,q):
    """Independent closed-interval intersection; preserves singleton witnesses."""
    speeds(ds);q=rational(q)
    if not 0<=q<=F(1,2):raise ValueError('clearance range')
    out=[(F(0),F(1))]
    for d in ds:
        nxt=[]
        for a,b in out:
            for j in range(d):
                lo=max(a,F(j,d)+q/d);hi=min(b,F(j+1,d)-q/d)
                if lo<=hi:nxt.append((lo,hi))
        out=sorted(set(nxt))
    return tuple(out)

def first_denominator(ds,cap):
    speeds(ds);integer(cap,2); threshold=F(1,len(ds)+1)
    for q in range(2,cap+1):
        ts=tuple(F(j,q) for j in range(1,q) if clearance(ds,F(j,q))>=threshold)
        if ts:return q,ts
    return None

def graph_edges(ds,n):
    speeds(ds);integer(n)
    return tuple((i,j) for i,j in combinations(range(n),2) if j-i in ds)

def alpha_brute(ds,n):
    edges=graph_edges(ds,n); masks=tuple((1<<i)|(1<<j) for i,j in edges)
    return max(mask.bit_count() for mask in range(1<<n) if all(mask&e!=e for e in masks))

def proper_period(ds,word):
    speeds(ds)
    if type(word) is not tuple or not word:raise ValueError('literal finite colour word')
    for x in word:integer(x)
    return all(word[r]!=word[(r+d)%len(word)] for r in range(len(word)) for d in ds)

def phase_step(r,z,d):
    integer(r);integer(z);integer(d,1)
    if r>=7 or z>=3:raise ValueError('retained residue/phase range')
    return (r+d)%7,(z+(r+d)//7)%3

def det(a):
    a=[list(row) for row in a];n=len(a);out=F(1)
    for j in range(n):
        p=next((i for i in range(j,n) if a[i][j]),None)
        if p is None:return F(0)
        if p!=j:a[p],a[j]=a[j],a[p];out=-out
        t=a[j][j];out*=t
        for i in range(j+1,n):
            q=a[i][j]/t
            for h in range(j+1,n):a[i][h]-=q*a[j][h]
    return out

def inv(a):
    n=len(a);b=[list(row)+[F(i==j) for j in range(n)] for i,row in enumerate(a)]
    for j in range(n):
        p=next(i for i in range(j,n) if b[i][j]);b[p],b[j]=b[j],b[p]
        t=b[j][j];b[j]=[x/t for x in b[j]]
        for i in range(n):
            if i!=j:
                t=b[i][j];b[i]=[x-t*y for x,y in zip(b[i],b[j])]
    return tuple(tuple(row[n:]) for row in b)

def incidence(i,j,n):return tuple(F(h==i)-F(h==j) for h in range(n-1))
def quad(a,b):return sum(x*a[i][j]*y for i,x in enumerate(b) for j,y in enumerate(b))

@dataclass(frozen=True)
class Network:
    source:tuple
    conductances:tuple
    covariance:tuple

def audit_network(receipt):
    if type(receipt) is not Network:raise ValueError('literal network receipt')
    v=speeds(receipt.source);n=len(v)
    if type(receipt.conductances) is not tuple or type(receipt.covariance) is not tuple:
        raise ValueError('immutable network arrays')
    edges=[]
    for row in receipt.conductances:
        if type(row) is not tuple or len(row)!=3:raise ValueError('edge field shape')
        i,j,w=row;integer(i);integer(j);w=rational(w)
        if not i<j<n or w<=0:raise ValueError('positive ordered edge')
        edges.append((i,j,w))
    if len({(i,j) for i,j,w in edges})!=len(edges):raise ValueError('duplicate edge')
    S=receipt.covariance
    if len(S)!=n-1 or any(type(row) is not tuple or len(row)!=n-1 for row in S):
        raise ValueError('covariance shape')
    for row in S:
        for x in row:rational(x)
    need(all(S[i][j]==S[j][i] for i in range(n-1) for j in range(n-1)),'symmetric covariance')
    L=tuple(tuple(sum(w*incidence(i,j,n)[a]*incidence(i,j,n)[b] for i,j,w in edges)
                  for b in range(n-1)) for a in range(n-1))
    need(all(det(tuple(row[:k] for row in L[:k]))>0 for k in range(1,n)),'positive grounded Laplacian')
    need(inv(L)==S,'exact inverse source network')
    capacities={};active={(i,j) for i,j,w in edges}
    for i,j in combinations(range(n),2):
        R=quad(S,incidence(i,j,n));cap=(F(1,v[i])+F(1,v[j]))**2
        need(R<=cap,'all original pair capacities')
        need((i,j) not in active or R==cap,'every positive edge is a contact')
        capacities[i,j]=(R,cap)
    tree_sum=F(0)
    for es in combinations(edges,n-1):
        B=tuple(tuple(incidence(i,j,n)[a] for i,j,w in es) for a in range(n-1))
        minor=det(B);need(minor in (-1,0,1),'independent incidence minor')
        weight=F(1)
        for i,j,w in es:weight*=w
        tree_sum+=minor*minor*weight
    need(tree_sum==det(L),'weighted spanning-tree determinant')
    return {'laplacian':L,'pair_resistance_capacities':[[i,j,*rc] for (i,j),rc in capacities.items()],
            'det_covariance':det(S),'tree_sum':tree_sum}

def reject(fn):
    try:fn()
    except (ValueError,TypeError,ZeroDivisionError,StopIteration):return True
    return False

def encode(x):
    if type(x) is F:return [x.numerator,x.denominator]
    if type(x) in (list,tuple):return [encode(y) for y in x]
    if type(x) is dict:return {str(k):encode(v) for k,v in x.items()}
    return x

def run():
    d=(1,2,4,7);old=(1,3,4,7);scaled=tuple(3*x for x in d)
    rows=[]
    for ds,expected in ((d,F(1,3)),(old,F(1,5)),(scaled,F(1,3))):
        q,ts=rotation(ds);need(q==expected,'exact source maximum')
        need(safe_intervals(ds,q)==tuple((t,t) for t in ts),'independent closed-witness reader')
        rows.append({'source':ds,'kappa':q,'maximizers':ts,'first_lonely_denominator':first_denominator(ds,20)})
    need(proper_period(d,(0,1,2)),'three-colour periodic upper')
    need(alpha_brute(d,3)==1,'literal triangle lower')
    oldword=tuple(2*(r//4)+(r%2) for r in range(8))
    need(proper_period(old,oldword),'old four-colour upper')
    need(all(abs(a-b) in old for a,b in combinations((0,3,4,7),2)),'old clique lower')
    need(rows[0]['first_lonely_denominator'][0]==3 and rows[2]['first_lonely_denominator'][0]==5,
         'same projected carrier loses absolute denominator')

    # Full source family over {1,...,9} containing 1,2; no empirical universal inference.
    triadic=[]
    for mask in range(1<<7):
        ds=(1,2)+tuple(j for j in range(3,10) if mask>>(j-3)&1)
        k,_=rotation(ds);mult=[x//3 for x in ds if x%3==0]
        need((k==F(1,3))==(not mult),'complete bounded triadic rigidity')
        need(proper_period(ds,(0,1,2))==(not mult),'bounded colour rigidity')
        if mult:need(k<=F(min(mult),3*min(mult)+1),'finite-window fractional barrier')
        triadic.append([ds,k])
    # All proper three-colourings of squared finite paths force period three.
    for n in range(3,9):
        valid=[]
        for word in product(range(3),repeat=n):
            if all(word[i]!=word[i+d] for d in (1,2) for i in range(n-d)):
                need(all(word[i]==word[i+3] for i in range(n-3)),'period three forced locally');valid.append(word)
        need(len(valid)==6,'exact palette gauges only')

    sharp=[]
    for m in range(1,13):
        q=3*m+1;ds=(1,2,3*m);t=F(m,q)
        need(rotation(ds)[0]==t,'sharp triadic-defect rotation')
        colour=tuple(((m*r)%q)//m for r in range(q))
        need(proper_period(ds,colour),'literal four-colour lift')
        I=tuple(3*j for j in range(m));sets=[{(i+j)%q for i in I} for j in range(q)]
        need(all(not any((r+d)%q in s for r in s for d in ds) for s in sets),'fractional independent sets')
        need(all(sum(r in s for s in sets)==m for r in range(q)),'exact fractional multiplicity')
        need(all(min((m*d)%q,q-(m*d)%q)>=m for d in ds),'circular multiplier certificate')
        if m<=5:need(alpha_brute(ds,q)==m,'independent full window bound')
        sharp.append({'m':m,'source':ds,'kappa':t,'fractional_and_circular':F(q,m),
                      'ordinary':4,'independent_set':I,'periodic_colouring':colour})

    # Paley orientation is a quotient with a winding, not the integer source graph.
    residues=(1,2,4)
    need(all(((j-i)%7 in residues)!=((i-j)%7 in residues) for i,j in combinations(range(7),2)),
         'intrinsic seven-vertex tournament')
    r=z=0;walk=[(r,z)]
    for step in (1,2,4):r,z=phase_step(r,z,step);walk.append((r,z))
    need(walk==[(0,0),(1,0),(3,0),(0,1)],'triangle winds once instead of source closure')
    for n in range(147):
        for step in d:
            r,z=phase_step(n%7,(n//7)%3,step)
            need((r,z)==((n+step)%7,((n+step)//7)%3),'exact source lift')
            need((r+z)%3!=(n%7+(n//7)%3)%3,'lifted three-colour predicate')
    for m in range(1,19):
        need((7*m%3==0)==(m%3==0),'minimum fibre-period obstruction')
        if m%3==0:need(proper_period(d,tuple(r%3 for r in range(7*m))),'legal finite phase cover')

    S=((F(64,49),F(0),F(0)),(F(0),F(81,196),F(1,392)),(F(0),F(1,392),F(121,784)))
    net=Network(d,((0,3,F(49,64)),(1,2,F(1,25)),(1,3,F(119,50)),(2,3,F(161,25))),S)
    electrical=audit_network(net)
    need(electrical['det_covariance']==F(200,2401),'exact grounded optimum')
    need(electrical['tree_sum']==F(2401,200),'three spanning trees retain their weights')
    A=F(2*(7+4),7*(4+2));B=F(4*(7+2),7*(4+2));c=F(1,99)
    need(A*A+B*B>1 and A*A+B*B-2*A*B*c==1,'one-unit threshold defect repaired')
    need(49<7*(4+2)+4*2 and 49>=7*(4+1)+4,'sole violating pair')
    need(F(4*(7*(4+2)+8),sum(x*x for x in d))*F(8,7)**2==F(1280,343),'Euclidean determinant conversion')
    L=electrical['laplacian'];b=incidence(1,2,4)
    cut=tuple(tuple(L[i][j]-F(1,25)*b[i]*b[j] for j in range(3)) for i in range(3))
    need(det(cut)>0 and quad(inv(cut),b)>F(9,16),'connected tree deletion loses feasibility')
    scale_net=Network(scaled,tuple((i,j,9*w) for i,j,w in net.conductances),
                      tuple(tuple(x/9 for x in row) for row in S))
    audit_network(scale_net)
    for i in range(4):
        for j in range(4):
            need(F(i==j)-F(d[i]*d[j],sum(x*x for x in d))==
                 F(i==j)-F(scaled[i]*scaled[j],sum(x*x for x in scaled)), 'identical literal orthogonal carrier')
    need(reject(lambda:rotation((1,True,4))),'reject bool speed')
    need(reject(lambda:rotation((1,2,2))),'reject duplicate speed')
    need(reject(lambda:clearance(d,.5)),'reject floating clock')
    need(reject(lambda:audit_network(replace(net,source=old))),'same covariance wrong source rejected')
    need(reject(lambda:audit_network(replace(net,conductances=((False,3,F(49,64)),)+net.conductances[1:]))),
         'reject Boolean edge alias')
    return {'status':'PROVED elementary source-specific structures; FINITE-EXACT controls; general LRC OPEN',
            'rotation_rows':rows,'triadic_source_universe':triadic,'sharp_family':sharp,
            'voltage_triangle':walk,'electrical':electrical,'network_source':net.source,
            'conductances':net.conductances,'covariance':S,'euclidean_determinant':F(1280,343),
            'checks':CHECKS}

if __name__=='__main__':
    data=(json.dumps(encode(run()),sort_keys=True,indent=2)+'\n').encode();p=HERE/'receipt.json'
    if '--write' in sys.argv:p.write_bytes(data)
    else:need(p.read_bytes().replace(b'\r\n',b'\n')==data,'immutable receipt comparison')
    print(json.dumps({'checks':json.loads(data)['checks'],'LF_SHA256':sha256(data).hexdigest()}))
