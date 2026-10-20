#!/usr/bin/env python3
"""Controls importing no producer: path masks and literal rational orbits."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def orbit(n,stop=None):
    x=n;w=[];xs=[n];A=d=1;B=C=0;first=None
    while True:
        raw=3*x+1;a=(raw&-raw).bit_length()-1;x=raw>>a;w.append(a);xs.append(x)
        B=3*B+d;A*=3;d<<=a;C=3*C+(d if a==1 else 0)
        if first is None and x<n:first=len(w)
        if (stop is None and first is not None) or len(w)==stop:break
        need(len(w)<10000,'orbit cap')
    # Literal compositions and reverse discounted emissions.
    actual=F(n-1);upper=F(n-1)
    for a in w:
        actual=3*actual/2**a+F(4,2**a)-1
        upper=3*upper/2**a+(1 if a==1 else 0)
    suffix=F(1);P=D=F(0)
    for a in reversed(w):
        if a==1:P+=suffix
        if a>=3:D+=(1-F(4,2**a))*suffix
        suffix*=F(3,2**a)
    need(actual==x-1 and upper==F(A*(n-1)+C,d),'literal affine check')
    need(upper-actual==D and upper==suffix*(n-1)+P,'suffix decomposition')
    reciprocal=sum((1-F(4,2**a))*F(x,xs[i+1]) for i,a in enumerate(w) if a>=3)
    if len(w)==first:need(D<=reciprocal,'first-exit reciprocal bound')
    return {'source':n,'endpoint':x,'steps':len(w),'clock':sum(w),'word':w,
            'first_descent_step':first,'paid':upper<n-1,'dividend':str(D),
            'dividend_over_gap':str(D/(n-x)),
            'reciprocal_budget':str(reciprocal),
            'threshold':str(1+F(C,d-A)) if d>A else None}

def graph_audit():
    rows=[]
    for n in range(1,7):
        es=list(combinations(range(n),2));pos={e:i for i,e in enumerate(es)};full=(1<<len(es))-1
        paths=set()
        for p in permutations(range(n)):
            paths.add(sum(1<<pos[tuple(sorted((p[i],p[i+1])))] for i in range(n-1)))
        # Superset closure of path masks, as opposed to subset/endvertex DP.
        traceable=bytearray(full+1)
        for p in paths:traceable[p]=1
        for bit in range(len(es)):
            b=1<<bit
            for m in range(full+1):
                if m&b and traceable[m^b]:traceable[m]=1
        connected=[]
        for mask in range(full+1):
            parent=list(range(n))
            def root(a):
                while a!=parent[a]:a=parent[a]
                return a
            for i,(a,b) in enumerate(es):
                if mask>>i&1:parent[root(a)]=root(b)
            connected.append(len({root(a) for a in range(n)})==1)
        stats=Counter()
        for m in range(full+1):
            h=full^m
            if not traceable[m] and not traceable[h]:
                stats['HH']+=1
                if connected[m] and connected[h]:stats['HH_both_connected']+=1
            if not traceable[m] and not connected[h]:stats['HW']+=1
            if not connected[m] and not traceable[h]:stats['WH']+=1
            if not connected[m] and not connected[h]:stats['WW']+=1
        need(stats['HH_both_connected']==0 and stats['WW']==0,'small-game boundary')
        rows.append({'n':n,'counts':dict(stats),'undirected_path_masks':len(paths)})
    return rows

def sharp_controls():
    checks=0;maxratio=(F(0),None)
    for m in range(1,7):
        for w in product(range(1,6),repeat=m):
            K=sum(w);A=3**m;d=2**K
            if d<=A:continue
            C=0;H=0
            for i,a in enumerate(w):
                H+=a
                if a==1:C+=3**(m-i-1)*2**H
            R=1+F(C,d-A)
            need(R/d<=F(161,256),'all-word sharp bound')
            need(C<=2**(K-m+2)*(3**(m-1)-2**(m-1)),'fixed-clock maximizer')
            if R/d>maxratio[0]:maxratio=R/d,w
            checks+=1
    need(maxratio==(F(161,256),(4,1,1,1,1)),'hostile maximum')
    return {'contracting_words_checked':checks,'universe':'length1..6, every exponent1..5',
            'maximum_ratio':str(maxratio[0]),'maximizer':maxratio[1]}

def hinge_controls():
    # a+b sqrt(-3), listed independently of the four-basis constructor.
    raw=[(-3,-1),(-1,-1),(-2,0),(-4,0),(-5,-1),(-4,-2),(-2,-2),(0,0),(-6,0),(-3,-3)]
    pts=[(F(a,2),F(b,2)) for a,b in raw]
    own=[];cross=[]
    for i,(a,b) in enumerate(pts):
        for j,(c,d) in enumerate(pts):
            if i<j and (a-c)**2+3*(b-d)**2==1:own.append((i,j))
            # Norm difference after rho rotation: rational part below;
            # coefficient of sqrt(33) vanishes exactly when b*c-a*d=0.
            if b*c==a*d and a*a+3*b*b+c*c+3*d*d-F(5,3)*(a*c+3*b*d)==1:
                cross.append((i,j))
    need(len(own)==18 and cross==[(0,0),(1,7),(2,7),(7,1),(7,2)],'separate hinge formula')
    total=((5**6+5)*2//3)*16
    single0=(4**6+4-5*(3**6+3)+10*(2**6+2)-10*2)*16
    singleother=single0*3//5
    blocked=single0**2+5*singleother**2
    need(total==166720 and blocked==836075520,'separate mask census')
    return {'copy_edges':own,'cross_edges':cross,'single_copy_port_count':total,
            'singleton_missing0':single0,'singleton_missing_other':singleother,
            'new_joint_blocks':blocked}

def analytic_constants():
    need(3**5<2**8,'log2(3)<8/5')
    need(sum(F(3,4)**i/__import__('math').factorial(i) for i in range(4))>2,'log2<3/4')
    need(F(144,25)>F(11,2) and F(7,5)**2<2,'radical lower bound')
    need(F(80,27)/F(56,15)==F(50,63),'uniform cone constant')
    need(63<13*65536 and F(25,9)<8,'base and monotonicity boundary')
    return {'source_floor':65536,'depth_bound':'m<=4 log(n)',
            'integer_depth_sufficient':'m<=2 floor(log2(n))',
            'threshold_upper':'R<1+50n/63<n',
            'cited_input':'Ellison Theorem3, x>=28; no numerical approximation used'}

if __name__=='__main__':
    samples=[orbit(n) for n in (3,27,63,4575,123027327)]+[orbit(55,39)]
    need(all(x['paid'] for x in samples[:-1]),'first-exit controls')
    need(not samples[-1]['paid'] and samples[-1]['endpoint']==53,'later-descent hostile')
    need(F(samples[3]['dividend'])>2 and F(samples[4]['dividend'])>4,'constant-dividend hostiles')
    print(json.dumps({'status':'FINITE-EXACT separate controls; same-agent implementation',
        'literal_orbits':samples,'games':graph_audit(),'sharp_bound':sharp_controls(),
        'hinge':hinge_controls(),'analytic_constants':analytic_constants()},indent=2,sort_keys=True))
