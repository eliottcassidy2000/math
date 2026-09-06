#!/usr/bin/env python3
"""Independent physical-predicate and error-polytope audit of all-parity 6/55.

The producer enumerates raw carriers and computes widths by fractional
knapsack. This verifier instead tests every literal sheet predicate between
rational breakpoints and takes widths over explicit error-cube vertices.
No producer imports, floating point, sampling, or optimization-disabled gates.
"""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import gcd, lcm
from pathlib import Path
import csv
import argparse
import sys

sys.stdout.reconfigure(newline="\n")
GATES=0
BASE=Path(__file__).resolve().parent


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def literal_tail(speeds):
    denominator=42*lcm(*speeds)
    end=denominator//3
    points={0,end}
    for w in speeds:
        unit=denominator//(14*w)
        for k in range(w):
            for sign in (-1,1):
                for j in range(3):
                    point=((14*k+sign)*unit-j*end)%denominator
                    if point<=end:
                        points.add(point)
    points=sorted(points)
    twice=2*denominator
    length=0
    for left,right in zip(points,points[1:]):
        x=left+right
        blocked=True
        for j in range(3):
            one=False
            for w in speeds:
                residue=w*(x+2*j*end)%twice
                if 14*min(residue,twice-residue)<twice:
                    one=True
                    break
            if not one:
                blocked=False
                break
        if blocked:
            length+=right-left
    return Q(3*length,denominator)


def head_audit(head_path):
    given=None
    if head_path:
        with Path(head_path).open(newline='',encoding='utf-8') as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        given={}
        for row in rows:
            w=tuple(int(row[x]) for x in ('a','b','c'))
            need(w not in given,('no repeated head row',w))
            given[w]=Q(int(row['mass_numerator']),int(row['denominator']))
    universe={w for w in combinations([x for x in range(1,64) if x%3],3) if gcd(*w)==1}
    need(len(universe)==10074,'complete all-parity primitive height63 universe')
    if given is not None:
        need(set(given)==universe,'raw table has exactly the complete all-parity universe')
    results={}
    for w in sorted(universe):
        observed=literal_tail(w)
        if given is not None:
            need(observed==given[w],('literal physical predicate versus raw carrier sum',w,observed,given[w]))
        need(observed<=Q(6,55),('head ceiling',w))
        results[w]=observed
    equalities=[w for w,x in results.items() if x==Q(6,55)]
    violations=[w for w,x in results.items() if x>Q(6,77)]
    first=min(violations,key=lambda w:(w[2],w[0],w[1]))
    need(equalities==[(1,10,11)],'unique head equality')
    need(len(violations)==151 and first==(2,5,7),'old ceiling failure count and first hostile')
    controls={(1,2,4):Q(0),(2,5,7):Q(22,245),(1,7,8):Q(31,392),
              (1,5,11):Q(6,77),(1,10,11):Q(6,55),(2,11,20):Q(11,140)}
    for w,x in controls.items():
        need(results[w]==x,('positive and hostile control',w))
    second=max(x for x in results.values() if x<Q(6,55))
    need(second==Q(383,3640),'second physical level')
    second_rows=[w for w,x in results.items() if x==second]
    digest=sha256()
    for w,x in results.items():
        digest.update(repr((w,x)).encode())
    need(digest.hexdigest()=='10dc0483d65ccec0c4999932983c387bcbb5c081a4c47ed3d76d625df2a24c0e',
         'frozen independent full literal-head certificate')
    print('LITERAL_HEAD',len(results),'max',Q(6,55),'equality',equalities,
          'old_bound_violations',len(violations),'first',first,'second',second,second_rows)
    print('LITERAL_HEAD_SEMANTIC_SHA256',digest.hexdigest())
    print('HEAD_PARITY_COUNTS',sorted(Counter(sum(x%2==0 for x in w) for w in universe).items()))


def cross(u,v):
    return (u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])


def owner_audit():
    cases=0
    for w in product((1,2),repeat=3):
        for v in product(range(3),repeat=3):
            if not any(v) or sum(x*y for x,y in zip(v,w))%3:
                continue
            need(sum(x==0 for x in v)<=1,'at most one residue-zero relation coefficient')
            images=[set() for _ in range(3)]
            live=[set() for _ in range(3)]
            for n in product(range(3),repeat=3):
                delta=sum(x*y for x,y in zip(v,n))%3
                C=tuple(x%3 for x in cross(w,n))
                images[delta].add(C)
                owners={(-ni*wi)%3 for ni,wi in zip(n,w)}
                if len(owners)==3:
                    live[delta].add(C)
                need((len(owners)==3)==all(C),'literal owners versus carrier units')
            unit=all(v)
            for delta in range(3):
                need(len(images[delta])==3,'full address-image affine fiber')
                expected=(2 if delta==0 else 0) if unit else (0 if delta==0 else 1)
                need(len(live[delta])==expected,'owner-word density and allowed integer defect')
                cases+=1
    need(cases==192,'complete owner-fiber universe')
    print('OWNER_FIBERS',cases,'method=literal nearest-integer owner words')


def plane_box_vertices(v,delta,low,high):
    """A 2D plane section vertex fixes two cube facets, solving the third."""
    vertices=set()
    for free in range(3):
        if not v[free]:
            continue
        other=[j for j in range(3) if j!=free]
        for values in product((low,high),repeat=2):
            z=[Q(0)]*3
            for j,value in zip(other,values):
                z[j]=Q(value)
            z[free]=Q(delta-sum(v[j]*z[j] for j in other),v[free])
            if low<=z[free]<=high:
                vertices.add(tuple(z))
    return vertices


def polygon_pattern(p):
    S=sum(p)
    unit=all(x%3 for x in p)
    radius=(3*S)//14
    ds=tuple(d for d in range(-radius,radius+1)
             if 14*abs(d)<3*S and (d%3==0 if unit else d%3!=0))
    B=Q(4*len(ds),3) if unit else Q(len(ds))
    best=Q(0)
    rho=Q(2 if unit else 1,3)
    for negative in range(3):
        if not p[negative]:
            continue
        v=tuple(-x if j==negative else x for j,x in enumerate(p))
        speed_vertices=plane_box_vertices(v,0,0,1)
        need(bool(speed_vertices),('nonempty closed normalized speed polytope',v))
        error_vertices={d:plane_box_vertices(v,14*d,-3,3) for d in ds}
        for d,vertices in error_vertices.items():
            need(bool(vertices),('nonempty strict error slice',v,d))
        pivot=next(j for j in range(3) if v[j])
        for w in speed_vertices:
            width_sum=Q(0)
            for vertices in error_vertices.values():
                values=[cross(w,z)[pivot]/(14*v[pivot]) for z in vertices]
                width=max(values)-min(values)
                need(width>=0,('error polygon width',v,w))
                width_sum+=width
            best=max(best,rho*width_sum)
    return ds,best,B


def box_audit():
    excluded={(0,1,1),(1,1,1),(1,1,2),(1,2,2)}
    patterns=[p for p in combinations_with_replacement(range(19),3)
              if sum(x!=0 for x in p)>=2 and gcd(*p)==1
              and sum(x%3==0 for x in p)<=1 and p not in excluded]
    need(len(patterns)==747,'parity-free coefficient box cardinality')
    need(Counter(sum(x!=0 for x in p) for p in patterns)=={2:49,3:698},'actual-zero coefficient sidecar')
    digest=sha256()
    maxima=[]
    empty=[]
    for p in patterns:
        ds,slope,B=polygon_pattern(p)
        need(slope<=Q(15,98),('global convex width slope',p,slope))
        need(B<=Q(2*sum(p),7)+Q(4,3),('defect intercept',p,B))
        if slope==Q(15,98):
            maxima.append(p)
        if not ds:
            empty.append(p)
            need(slope==0 and B==0,'empty defect list must be separately discharged')
        digest.update(repr((p,ds,slope,B)).encode())
    need(maxima==[(1,7,8)],'complete slope equality list')
    need(empty==[(0,1,2)],'unique empty-defect parity corner')
    need(digest.hexdigest()=='cf808062354debbefc1d8ead8ad0d10e9da5427cb42b8f083b6af24d0059c87c',
         'independent error polygon certificate matches knapsack table')
    exception_rows=[]
    for p,expected in (((1,1,1),Q(2,7)),((1,1,2),Q(2,7)),((1,2,2),Q(3,14))):
        row=polygon_pattern(p)
        need(row==((0,),expected,Q(4,3)),('low circuit count exception',p,row))
        exception_rows.append((p,row))
    print('ERROR_POLYGON_BOX',len(patterns),'supports',dict(Counter(sum(x!=0 for x in p) for p in patterns)),
          'slope',Q(15,98),'equality',maxima,'empty_defects',empty)
    print('BOX_SEMANTIC_SHA256',digest.hexdigest())
    print('LOW_CIRCUITS',exception_rows)


def analytic_gates():
    gap=Q(14,55)-Q(15,98)
    A=Q(2,7)/gap
    B=Q(4,3)/gap
    need((gap,A,B)==(Q(547,5390),Q(1540,547),Q(21560,1641)),'exact count-gap conversion')
    need(A*18+B==Q(104720,1641)<64,'S<=18 tail')
    h19=Q(3*19**2,16)-A*19-B
    need(h19==Q(27763,26256)>0 and Q(3*19,8)-A>0,'S>=19 monotone quadratic tail')
    need(Q(6,49)+Q(4,7*19)==Q(142,931)<Q(15,98),'large coefficient complementary case')
    need(Q(3,14)+Q(4,3*34)<Q(14,55),'norm5 c34 tail')
    need(Q(9,98)+Q(6,7*50)<Q(6,55),'norm3 c50 tail')
    need(Q(3,49)+Q(6,7*18)<Q(6,55),'norm4 c18 tail')
    print('STRICT_TAILS general=64 norm3=50 norm4=18 norm5=34')
    print('REPAIRED_CORNER (0,1,2): no admissible defect, hence N=physical_mass=0; do not assert 0<0')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--head-tsv',help='Optional independently produced raw table for rowwise comparison')
    args=parser.parse_args()
    print('SCOPE physical Haar failure, primitive distinct positive ternary units, all parities')
    print('PRODUCER_SOURCE_SHA256',sha256((BASE/'overnight4_20260906_lrc_parityfree_probe.py').read_bytes()).hexdigest())
    if args.head_tsv:
        print('RAW_HEAD_SHA256',sha256(Path(args.head_tsv).read_bytes()).hexdigest())
    else:
        print('RAW_HEAD_COMPARISON omitted; entire literal head is still independently verified')
    head_audit(args.head_tsv)
    owner_audit()
    box_audit()
    analytic_gates()
    print('PASS_OPTIMIZATION_LIVE_GATES',GATES)


if __name__=='__main__':
    main()
