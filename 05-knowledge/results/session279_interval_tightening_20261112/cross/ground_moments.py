"""Ground-state correction and sharp degree-four discrete moment intervals.

The consumer certify(packet) uses no tournament count, partitions, or gcd.
The experimental producer deliberately enumerates odd partitions: this is a
certificate-information experiment, not a faster enumeration claim.
"""
from fractions import Fraction as F
from math import comb, gcd, factorial
from collections import defaultdict, Counter
from functools import lru_cache
import argparse, json, pathlib, sys, time, hashlib

DEGREE=4

def elementary(q, degree=DEGREE):
    es=[F(1)]+[F(0)]*degree
    for x in q:
        for j in range(degree,0,-1): es[j]+=x*es[j-1]
    return es

def parts(n,lo=1):
    if not n:yield ()
    for a in range(lo,n+1,2):
        for p in parts(n-a,a):yield (a,)+p

def sector(p,backbone):
    e=sum((a-1)//2 for a in p); ds=[]
    for i,a in enumerate(p):
        for b in p[i+1:]:
            g=gcd(a,b);base=g if a in backbone or b in backbone else min(a,b)
            e+=base
            if base>g:ds.append(base-g)
    z=1
    for a,c in Counter(p).items():z*=a**c*factorial(c)
    return F(1<<e,z),ds

def independent_exact(p):
    e=sum((a-1)//2 for a in p)+sum(gcd(a,b) for i,a in enumerate(p) for b in p[i+1:])
    z=1
    for a,c in Counter(p).items():z*=a**c*factorial(c)
    return F(1<<e,z)

def produce(n,backbone=(3,),reference=False):
    cells=defaultdict(lambda:[F(0)]*5); oracle=F(0);count=0
    for p in parts(n):
        count+=1;w,ds=sector(p,set(backbone));m=len(ds);delta=min(ds,default=0)
        r=sum(d>delta for d in ds);es=elementary([1-F(1,1<<d) for d in ds])
        for j in range(5):cells[m,delta,r][j]+=w*es[j]
        if reference:oracle+=independent_exact(p)
    packet={'schema':'ground-moments-v1','n':n,'backbone':list(backbone),'degree':4,
            'odd_partitions':count,'cells':[{'M':m,'delta':d,'R':r,'C':[str(x) for x in cs]} for (m,d,r),cs in sorted(cells.items())]}
    if reference and oracle.denominator!=1:raise AssertionError('nonintegral oracle')
    return packet,oracle if reference else None

@lru_cache(None)
def coefficients(roots):
    vals=[]
    for x in range(len(roots)+1):
        v=F(1)
        for r in roots:v*=F(r-x,r)
        vals.append(v)
    cs=[]
    while vals:
        cs.append(vals[0]);vals=[b-a for a,b in zip(vals,vals[1:])]
    return tuple(cs)

@lru_cache(None)
def duals(m):
    if m<=4:return (),()
    # Full-degree vertices suffice, but include lower-degree certificates.
    upper=[(),(m,)]+[(r,r+1) for r in range(1,m)]+[(r,r+1,m) for r in range(1,m-1)]
    upper += [(r,r+1,s,s+1) for r in range(1,m-2) for s in range(r+2,m)]
    lower=[(1,),(1,m)]+[(1,r,r+1) for r in range(2,m)]
    lower += [(1,r,r+1,m) for r in range(2,m-1)]
    return tuple(lower),tuple(upper)

def moment_interval(m,cs):
    if m<=4:
        v=sum(((-1)**j*cs[j] for j in range(m+1)),F(0))
        return v,v,tuple(range(1,m+1)),tuple(range(1,m+1))
    lower,upper=duals(m)
    def val(rs):return sum((a*b for a,b in zip(coefficients(rs),cs)),F(0))
    l,lr=max([(F(0),None)]+[(val(rs),rs) for rs in lower],key=lambda x:x[0])
    u,ur=min((val(rs),rs) for rs in upper)
    return l,u,lr,ur

def integer_interval(lo,hi):
    il=-(-lo.numerator//lo.denominator);ih=hi.numerator//hi.denominator
    if il>ih:raise AssertionError('no integer in certified interval')
    return {'low':str(il),'high':str(ih),'width':str(ih-il),'residue_bits_sufficient':(ih-il).bit_length(), 'integer_exact':il==ih, 'recovered_integer':str(il) if il==ih else None}

def certify(packet):
    """Pure certificate consumer. No exact result or partition oracle is read."""
    if packet['degree']!=4 or packet['schema']!='ground-moments-v1':raise ValueError('schema')
    original=defaultdict(lambda:[F(0)]*5);ground=defaultdict(lambda:[F(0)]*5)
    coarse=defaultdict(lambda:[F(0)]*5)
    for cell in packet['cells']:
        m,d,r=cell['M'],cell['delta'],cell['R'];cs=list(map(F,cell['C']))
        if not (0<=r<=m and d>=0 and len(cs)==5 and all(x>=0 for x in cs)):raise ValueError('invalid cell')
        if m==0 and (d or r):raise ValueError('empty cell')
        if m>0 and d==0:raise ValueError('zero positive defect')
        rho=F(1,1<<d);a=1-1/rho;b=1/rho;scale=rho**m
        for k in range(5):
            original[m][k]+=cs[k]
            z=scale*sum((F(comb(m-j,k-j))*a**(k-j)*b**j*cs[j] for j in range(min(k,m)+1) if k-j<=m-j),F(0)) if k<=m else F(0)
            if z<0 or (k>r and z):raise ValueError('inconsistent transformed moment')
            ground[r][k]+=z;coarse[m][k]+=z
    out={}
    for name,data in [('uncentered_stratified',original),('ground_without_R',coarse),('ground_with_R',ground)]:
        lo=hi=F(0);roots=[]
        for m,cs in sorted(data.items()):
            l,u,lr,ur=moment_interval(m,cs);lo+=l;hi+=u
            roots.append({'support':m,'lower_roots':lr,'upper_roots':ur})
        out[name]={**integer_interval(lo,hi),'rational_low':str(lo),'rational_high':str(hi),'roots':roots}
    # Prior S278 bounds, reconstructed exactly from the same uncentered packet.
    total=[sum((cs[j] for cs in original.values()),F(0)) for j in range(5)]
    lowrema=[sum((F(k+1,m)*cs[k+1] for m,cs in original.items() if m>k),F(0)) for k in range(4)]
    bs=[sum(((-1)**j*total[j] for j in range(k+1)),F(0)) for k in range(4)]
    lo=max(bs[k]+lowrema[k] if k%2 else bs[k]-total[k+1] for k in range(4))
    hi=min(bs[k]+total[k+1] if k%2 else bs[k]-lowrema[k] for k in range(4))
    out['legacy_same_packet']={**integer_interval(lo,hi),'rational_low':str(lo),'rational_high':str(hi)}
    return out

def controls():
    # Every displayed dual is checked independently on its full integer support.
    checks=0
    for m in range(5,41):
        ls,us=duals(m)
        for roots,upper in [(r,False) for r in ls]+[(r,True) for r in us]:
            coeff=coefficients(roots)
            for x in range(m+1):
                v=sum((c*comb(x,j) for j,c in enumerate(coeff) if j<=x),F(0));target=int(x==0)
                if (upper and v<target) or (not upper and v>target):raise AssertionError((m,roots,x))
                checks+=1
    # Affine elementary-moment transform against explicit edge lists.
    for ds in [[],[4],[4,4,6],[4,6,8,4,6],[6,8,12,6,10,8]]:
        m=len(ds);d=min(ds,default=0);rho=F(1,1<<d);q=[1-F(1,1<<x) for x in ds];cs=elementary(q)
        packet={'schema':'ground-moments-v1','degree':4,'cells':[{'M':m,'delta':d,'R':sum(x>d for x in ds),'C':list(map(str,cs))}]}
        exact=F(1)
        for x in ds:exact*=F(1,1<<x)
        # Integer endpoint handling is not appropriate for an arbitrary probability.
        zs=elementary([1-F(1,1<<(x-d)) for x in ds])
        a=1-1/rho;b=1/rho
        for k in range(5):
            trans=sum((comb(m-j,k-j)*a**(k-j)*b**j*cs[j] for j in range(min(k,m)+1)),F(0)) if k<=m else F(0)
            if trans!=zs[k]:raise AssertionError('transform')
    # First exceptional-rank-five cores: exhaustive finite structural controls.
    thresholds=[]
    for backbone,expected in [((3,),39),((3,5),48),((3,5,7),48),((3,5,7,9,11),70)]:
        witness=None
        for n in range(max(backbone)+2,expected+1):
            for p in parts(n,max(backbone)+2):
                _,ds=sector(p,set(backbone));delta=min(ds,default=0);r=sum(d>delta for d in ds)
                if r>=5:witness=(n,p,ds,r);break
            if witness:break
        if witness is None or witness[0]!=expected:raise AssertionError('threshold')
        thresholds.append({'backbone':backbone,'first_R5':witness})
    return {'dual_point_checks':checks,'affine_transform_controls':5,'residual_thresholds':thresholds}

def experiment():
    rows=[]
    for n in [24,32,38,39,40,48,56,64,69,70]:
        for backbone in ([(3,5,7,9,11)] if n>=69 else [(3,),(3,5),(3,5,7),(3,5,7,9,11)] if n>=40 else [(3,)]):
            t=time.perf_counter();packet,oracle=produce(n,backbone,True);result=certify(packet)
            for v in result.values():
                if not F(v['rational_low'])<=oracle<=F(v['rational_high']):raise AssertionError('oracle containment')
            row={'n':n,'backbone':list(backbone),'odd_partitions':packet['odd_partitions'],'moment_cells':len(packet['cells']),
                 'exact_bits':int(oracle).bit_length(),'width_bits':{k:v['residue_bits_sufficient'] for k,v in result.items()},
                 'seconds_produce_and_certify':round(time.perf_counter()-t,6)}
            if n==64 and backbone==(3,):
                row['certificates']=result
            rows.append(row)
    return {'status':'PASS','controls':controls(),'rows':rows}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--consume');ap.add_argument('--write-packet');ap.add_argument('--n',type=int,default=64);ap.add_argument('--backbone',default='3');ap.add_argument('--experiment',action='store_true');args=ap.parse_args()
    if args.consume:result=certify(json.loads(pathlib.Path(args.consume).read_text()))
    elif args.write_packet:
        result,_=produce(args.n,tuple(map(int,args.backbone.split(','))),False)
        pathlib.Path(args.write_packet).write_text(json.dumps(result,indent=2)+'\n');print('wrote moment packet without exact reference');sys.exit()
    else:result=experiment()
    print(json.dumps(result,indent=2))
