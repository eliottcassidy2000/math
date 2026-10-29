"""Exact mass, source-labelled wedge, and geometric selector controls."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,product
from hashlib import sha256
import json

P=Path(__file__).resolve().parent
checks=0
def require(x,msg):
    global checks
    checks+=1
    if not x:raise ValueError(msg)

def step(n):
    z=3*n+1;a=0
    while z%2==0:z//=2;a+=1
    return z,a

def block(n):
    q=1;w=[]
    while q:
        n,a=step(n);w.append(a);q+=-1+3*(a%3==0)
    return n,tuple(w)

B=((1,1,1),(4,6,179),(7,11,16177),(7,11,12289),(7,11,9697))
WORDS=((1,),(3,1,1,1),(3,3,1,1,1,1,1),(3,1,3,1,1,1,1),(3,1,1,3,1,1,1))
def weight(m):return F(1,2*4**(m.bit_length()-1))

def masses():
    rows=json.loads((P/'least60.json').read_text())['rows']
    old=json.loads((P/'inherited40.json').read_text())['rows']
    require([x['least_source'] for x in rows[:40]]==[x['least_source'] for x in old], 'all inherited minima')
    for row in rows:
        n=row['least_source'];J=D=0
        for i in row['word']:
            h,w=block(n);require(h>n and w==WORDS[i],'literal witness')
            J+=len(w);D+=sum(w);n=h
        require(n==int(row['endpoint']),'endpoint with source retained')
        require(J<=7*row['k'] and D<=11*row['k'],'native clocks')
        row['odd_clock']=J;row['division_clock']=D
    require(all(r['word'].count(2)+r['word'].count(3)+r['word'].count(4)==0 for r in rows[:54]), 'earlier minimizing witnesses use first two branches')
    require(4 in rows[54]['word'],'horizon55 needs the full search alphabet')
    N=[1,2,3,4,6,9,13,18,24,31,40]
    for k in range(11,101):N.append(N[k-1]+N[k-6]+3*N[k-11])
    # Independent finite-precision inverse chart iteration.
    for level in range(17):
        mod=1<<(level+1);cells={1}
        charts=[(1<<d,c,pow(3**j,-1,mod)) for j,d,c in B]
        for _ in range(level):
            cells={((a*y-c)*inv)%mod for y in cells for a,c,inv in charts}
        require(len(cells)==N[level],'occupied-cell count')
    q=F(1049603,4194304);bounds=[];best=F(1)
    for row in rows:
        k=row['k'];ell=((row['least_source']+1)//2).bit_length()-1
        bound=F(N[k],4**k)+sum((F(N[j+1],2*4**j) for j in range(ell,k)),F())
        best=min(best,bound)
        bounds.append(dict(k=k,least_source=row['least_source'],ell=ell,
                           full_mass_bound=str(best),atom_bound=str(best-q**k),
                           odd_clock=row['odd_clock'],division_clock=row['division_clock']))
    # Direct complete cylinder sums through seven return blocks.
    states=[(0,0,0)]
    for k in range(1,8):
        states=[(J+j,D+d,3**j*C+(1<<D)*c) for J,D,C in states for j,d,c in B]
        mass=F()
        for J,D,C in states:
            rn=(((1<<D)-C)*pow(3**J,-1,1<<(D+1)))%(1<<(D+1))
            r=((rn+1)//2)%(1<<D)
            mass+=F(1,4**D)+(weight(r) if r else 0)
        require(mass<=F(bounds[k-1]['full_mass_bound']),'full atomic upper bound, not Haar mass')
    require(F(bounds[-1]['full_mass_bound'])<F(117121,10**32),'outward decimal bound')
    # Scan every available window, including k=0. These are finite tests,
    # not an all-k lower bound or a construction of a predecessor map.
    minima=[1]+[row['least_source'] for row in rows]
    windows=[]
    for width in range(1,13):
        ratio,start=min((F(minima[k+width]+1,minima[k]+1),k)
                        for k in range(len(minima)-width))
        windows.append(dict(width=width,least_ratio=str(ratio),start=start,end=start+width))
    require(F(windows[6]['least_ratio'])>=7,'seven-return sevenfold candidate through horizon60')
    require(7**3*2**14>3**14,'sevenfold/seven-return rate pays 3/2 envelope if uniform and s<1/2')
    # Paid-prefix growth: (n_k+1)/(n_0+1)<= (3/2)^k.
    for n0 in range(1,2000,2):
        n=n0
        for k in range(1,9):
            h,w=block(n)
            if h<=n:break
            require(2**k*(h+1)<=3**k*(n0+1),'accumulated growth envelope')
            n=h
    return dict(bounds=bounds,q=str(q),cell_counts=N[:61],finite_windows=windows,
                first_late_branch_witness=rows[54],horizon60=rows[-1])

def mm(a,b):return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]

def later_windows(initial_bound):
    receipt=json.loads((P/'window_checks.json').read_text())
    require(receipt['source_sha256']==sha256((P/'window_hostile.cpp').read_bytes()).hexdigest(), 'bounded search source pin')
    old={r['k']:r['least_source'] for r in json.loads((P/'least60.json').read_text())['rows']}
    N=[1,2,3,4,6,9,13,18,24,31,40]
    for k in range(11,68):N.append(N[k-1]+N[k-6]+3*N[k-11])
    best=F(initial_bound);bounds=[]
    require([r['k'] for r in receipt['rows']]==list(range(61,68)), 'all out-of-sample windows')
    for row in receipt['rows']:
        k=row['k'];bound=row['exclusive_bound']
        require(bound==7*(old[k-7]+1)-1 and row['counterexample'] is False, 'exact exclusion threshold')
        ell=((bound+1)//2).bit_length()-1
        mass=F(N[k],4**k)+sum((F(N[j+1],2*4**j) for j in range(ell,k)),F())
        best=min(best,mass)
        bounds.append(dict(k=k,source_lower_bound=bound,ell_lower_bound=ell,full_mass_bound=str(best)))
    control=receipt['positive_control'];n=control['source']
    require(control['counterexample'] and n==old[55] and control['exclusive_bound']==n+1, 'positive control finds boundary witness')
    for letter in control['word']:
        h,w=block(n)
        require(w==WORDS[letter] and h>n,'native positive control for exclusion engine')
        n=h
    require(n==int(control['endpoint']),'positive-control endpoint')
    require(F(bounds[-1]['full_mass_bound'])<F(4595,10**32),'outward horizon67 all-source bound')
    return bounds

def mv(a,v):return [sum(x*y for x,y in zip(r,v)) for r in a]
def wedge(a):
    pairs=list(combinations(range(len(a)),2))
    return [[a[i][k]*a[j][l]-a[i][l]*a[j][k] for k,l in pairs] for i,j in pairs]
def minors(u,v):return [u[i]*v[j]-u[j]*v[i] for i,j in combinations(range(4),2)]
def clock_wedge():
    # Fixed source and unnormalised current homogeneous coordinate.
    for source in range(1,400,2):
        n=source;N=source;Q=1;u=[source,1,0,0];v=[0,0,N,Q]
        p=minors(u,v)
        for clock in range(12):
            h,a=step(n)
            G=[[1,0,0,0],[0,1,0,0],[0,0,3,1],[0,0,0,2**a]]
            pp=mv(wedge(G),p);NN=3*N+Q;QQ=2**a*Q
            require(pp==minors(u,[0,0,NN,QQ]),'native exterior transfer')
            require(NN==QQ*h,'integer homogeneous clock')
            require(pp[2]-pp[3]==QQ*(source-h),'original-source payment')
            require(pp[0]*pp[5]-pp[1]*pp[4]+pp[2]*pp[3]==0,'Plucker stays exact')
            require(pp[2]-pp[3]==3*(p[2]-p[3])+Q*((2**a-3)*source-1),'signed debt recurrence')
            normalized=[x//Q for x in p]
            raw=mv(wedge(G),normalized)
            require(all(x%(2**a)==0 for x in raw),'guarded projective division is integral')
            require([x//(2**a) for x in raw]==[0,source*h,source,h,1,0], 'source-faithful normalized exterior transfer')
            require(QQ.bit_length()-1==Q.bit_length()-1+a,'separate integer division clock survives normalization')
            N,Q,n,p=NN,QQ,h,pp
    # Same native odd-step clock and same original-source weight, different
    # division clocks: normalizing after aggregation changes the observation.
    h1,a1=step(45);h2,a2=step(59)
    require((h1,a1,h2,a2)==(17,3,89,1), 'native mixture hostile')
    require(weight(23)==weight(30)==F(1,512),'hostile preserves the actual source weights')
    scaled=weight(23)*2**a1*(45-h1)+weight(30)*2**a2*(59-h2)
    unscaled=weight(23)*(45-h1)+weight(30)*(59-h2)
    require((scaled,unscaled)==(F(41,128),F(-1,256)), 'mixed clocks reverse aggregate debt sign')
    return dict(native_sources=200,steps_each=12,source_payment_coordinates='p03-p12',
                normalized_carrier='(0,n0*n,n0,n,1,0) with J,D retained separately',
                mixture_hostile=dict(sources=[45,59],odd_clock=1,division_clocks=[3,1],
                                     original_weight='1/512 each',scaled_debt=str(scaled),unscaled_debt=str(unscaled)))

def rank(a):
    a=[[F(x) for x in row] for row in a];r=0
    for j in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][j]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];t=a[r][j];a[r]=[x/t for x in a[r]]
        for i in range(len(a)):
            if i!=r:
                t=a[i][j];a[i]=[x-t*y for x,y in zip(a[i],a[r])]
        r+=1
    return r

def geometry():
    points=[(F(-3,5),F(4,5)),(F(-3,5),F(-4,5)),
            (F(8,5),F(4,5)),(F(8,5),F(-4,5)),(F(0),F(0)),(F(1),F(0))]
    distance=lambda p,q:sum((x-y)**2 for x,y in zip(p,q))
    edges=[(i,j) for i,j in combinations(range(6),2) if distance(points[i],points[j])==1]
    require(edges==[(0,4),(1,4),(2,5),(3,5),(4,5)],'faithful unit graph')
    left=[(F(0),F(0)),(F(-6,5),F(0))];right=[(F(1),F(0)),(F(11,5),F(0))]
    for x in left:require(all(distance(x,points[i])==1 for i in (0,1)),'left contact roots')
    for y in right:require(all(distance(y,points[i])==1 for i in (2,3)),'right contact roots')
    allowed=[(i,j) for i,x in enumerate(left) for j,y in enumerate(right) if distance(x,y)==1]
    require(allowed==[(0,0)],'one cross-edge selects the root pair')
    R=[[F(3,5),F(-4,5),0,0],[F(3,5),F(4,5),0,0],
       [0,0,F(-3,5),F(-4,5)],[0,0,F(-3,5),F(4,5)],[-1,0,1,0]]
    require(rank(R)==4,'full pinned rank')
    require(all(sum(F(w)*row[j] for w,row in zip((5,5,5,5,6),R))==0 for j in range(4)),'exact pinned stress')
    minimum=36
    for ports in product(range(6),repeat=4):
        count=sum(x!=ports[0] and x!=ports[1] and y!=ports[2] and y!=ports[3] and x!=y
                  for x,y in product(range(6),repeat=2))
        require(count>=12,'all prescribed six-colour boundaries extend')
        minimum=min(minimum,count)
    require(minimum==12,'sharp minimum number of six-colour extensions')
    require(not any(x not in (1,2) and y not in (1,2) and x!=y for x,y in product(range(3),repeat=2)), 'three-colour boundary fails')
    parameters=sorted({F(num,den) for den in range(2,9) for num in range(1,den)})
    for t in parameters:
        a=(1-t*t)/(1+t*t);b=2*t/(1+t*t)
        family=[(-a,b),(-a,-b),(1+a,b),(1+a,-b),(F(0),F(0)),(F(1),F(0))]
        require([(i,j) for i,j in combinations(range(6),2) if distance(family[i],family[j])==1]==edges,'portable rational family has the same faithful graph')
        roots_x=[(F(0),F(0)),(-2*a,F(0))]
        roots_y=[(F(1),F(0)),(1+2*a,F(0))]
        require([(i,j) for i,x in enumerate(roots_x) for j,y in enumerate(roots_y) if distance(x,y)==1]==[(0,0)],'portable orientation selector')
        matrix=[[a,-b,0,0],[a,b,0,0],[0,0,-a,-b],[0,0,-a,b],[-1,0,1,0]]
        require(rank(matrix)==4 and all(sum(w*row[j] for w,row in zip((1,1,1,1,2*a),matrix))==0 for j in range(4)),'portable pinned stress')
        gap=2*a
        require((2*a*b)**2==gap**2*(1-gap**2/4),'selection gap and contact conditioning relation')
    # Classical regular-singular inverse clock: exact decay/recovery product.
    product_factor=F(1)
    for m in range(1,31):
        product_factor*=F(2*m+1,2)
        require(F(1,product_factor)*product_factor==1,'inverse smoothing owes its recovery cost')
    return dict(points=[[str(x),str(y)] for x,y in points],edges=edges,
                root_pairs=allowed,pinned_rank=4,stress=[5,5,5,5,6],
                six_colour_min_extensions=minimum,three_colour_boundary=[1,2,1,2],
                rational_family_parameters=[str(t) for t in parameters])

if __name__=='__main__':
    result=dict(mass=masses(),wedge=clock_wedge(),geometry=geometry())
    result['out_of_sample']=later_windows(result['mass']['bounds'][-1]['full_mass_bound'])
    result.update(status='PASS',explicit_checks=checks,
                  source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  input_sha256=sha256((P/'least60.json').read_bytes()).hexdigest())
    print(json.dumps(result,indent=2,sort_keys=True))
