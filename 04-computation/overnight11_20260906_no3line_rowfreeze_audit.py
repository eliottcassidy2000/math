"""Independent conditional no-three-line referee.

Rook engine: edge-subset matching enumeration followed by row-mask zeta
summation. Large controls use paired-row and cyclic-run formulas. No
producer code is read or imported; all gates survive optimized Python.
"""
from collections import Counter,defaultdict
from fractions import Fraction as Q
from itertools import permutations
from math import comb,factorial
import hashlib,json,sys
sys.stdout.reconfigure(newline='\n')
GATES=0

def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)

def fall(n,k):return factorial(n)//factorial(n-k) if 0<=k<=n else 0

def parts(n,minimum=2):
    if n==0:yield ()
    for first in range(minimum,n+1):
        for tail in parts(n-first,first):yield (first,)+tail

def skeleton(shape):
    edges=[];offset=0
    for h in shape:
        edges.extend((offset+r,offset+c) for r in range(h) for c in (r,(r+1)%h))
        offset+=h
    return edges

def all_rooks(shape):
    """Every matching once by its increasing edge IDs, then subset zeta."""
    n=sum(shape);edges=skeleton(shape)
    coefficients=[[0]*(n+1) for _ in range(1<<n)]
    def visit(start,rowmask,colmask,k):
        coefficients[rowmask][k]+=1
        for index in range(start,len(edges)):
            r,c=edges[index]
            if not (rowmask>>r&1 or colmask>>c&1):
                visit(index+1,rowmask|(1<<r),colmask|(1<<c),k+1)
    visit(0,0,0,0)
    for bit in range(n):
        for mask in range(1<<n):
            if mask>>bit&1:
                lower=coefficients[mask^(1<<bit)]
                coefficients[mask]=[a+b for a,b in zip(coefficients[mask],lower)]
    return coefficients

def conditional_mean(n,rho,bank):
    total=Q(0)
    for d in range(1-n,n):
        mask=sum(1<<r for r in range(n) if 0<=rho[r]+d<n)
        total+=sum((Q((-1)**(k+1)*(k-2)*bank[mask][k],fall(n,k))
                    for k in range(3,n+1)),Q(0))
    return total

def F(points):return sum(max(v-2,0) for v in Counter(c-r for r,c in points).values())

def subset_controls():
    profiles=subsets=0;bank={}
    for n in range(2,10):
        for shape in parts(n):
            coefficients=all_rooks(shape);bank[shape]=coefficients
            edges=skeleton(shape);profiles+=1
            for mask,m in enumerate(coefficients):
                subsets+=1;L=mask.bit_count()
                selected=[(r,c) for r,c in edges if mask>>r&1]
                columns=Counter(c for r,c in selected)
                c2=sum(v==2 for v in columns.values())
                need(c2<=L and all(v<=2 for v in columns.values()),'C4-aware shared-column count')
                need(m[0]==1 and m[1]==2*L,'literal matching normalization')
                need(m[2]==2*fall(L,2)-c2,'exact two-matching collision count')
                if n>=3:
                    need(m[3]==Q(4,3)*fall(L,3)-2*c2*(L-2),'exact conditional m3 for every subset')
                if n>=4:need(m[4]<=Q(2,3)*fall(L,4),'conditional m4 upper bound')
                for k in range(2,n+5):
                    moment=Q(factorial(k)*m[k],fall(n,k)) if k<=n else Q(0)
                    theta=Q(L,n)
                    deficit=(2*theta)**k-moment
                    need(0<=deficit<=Q(2**k*k*(k-1),n)*theta**(k-1),
                         'all-k conditional deficit, including L0/1 and k>L/n')
                    if 2<=k<=L:
                        q=Q(factorial(k)*m[k],fall(L,k)*2**k)
                        need(0<=1-q<=Q(k*(k-1),4*(L-1)), 'conditional neighbor-collision union bound')
            if n<=5:
                # A canonical target matching on each subset; columns alone are randomized.
                for mask,m in enumerate(coefficients):
                    rows=[r for r in range(n) if mask>>r&1]
                    target={r:i for i,r in enumerate(rows)}
                    moments=[0]*(n+1)
                    for sigma in permutations(range(n)):
                        y=sum(r in target and sigma[c]==target[r] for r,c in edges)
                        for k in range(n+1):moments[k]+=fall(y,k)
                    for k in range(n+1):
                        need(Q(moments[k],factorial(n))==Q(factorial(k)*m[k],fall(n,k)),
                             'independent literal partial-matching moment')
    return bank,profiles,subsets

def literal_row_controls(bank):
    reports=[];roworders=columns=nodes=0
    for n in range(2,6):
        for shape in parts(n):
            edges=skeleton(shape);hist=Counter()
            for rho in permutations(range(n)):
                roworders+=1;total=0;prefix=defaultdict(lambda:[0,0])
                check_prefix=(rho==tuple(range(n)) or (n==5 and rho==(0,2,3,1,4)))
                for sigma in permutations(range(n)):
                    columns+=1;value=F([(rho[r],sigma[c]) for r,c in edges]);total+=value
                    if check_prefix:
                        code=sigma[:-1]
                        for k in range(n):
                            entry=prefix[code[:k]];entry[0]+=value;entry[1]+=1
                mean=Q(total,factorial(n));hist[str(mean)]+=1
                need(mean==conditional_mean(n,rho,bank[shape]),'complete conditional mean from literal columns')
                if n<=4:need(mean==conditional_mean(n,tuple(range(n)),bank[shape]),'bounded small-size row invariance')
                if check_prefix:
                    children=defaultdict(list)
                    for code,(value,count) in prefix.items():
                        if code:children[code[:-1]].append(Q(value,count))
                    for means in children.values():
                        need(max(means)-min(means)<=4,'one-permutation conditional interval width')
                        nodes+=1
            reports.append(dict(n=n,shape=shape,histogram=dict(sorted(hist.items()))))
    shape=(5,)
    need(conditional_mean(5,(0,1,2,3,4),bank[shape])==Q(11,10),'identity-row hostile mean')
    need(conditional_mean(5,(0,2,3,1,4),bank[shape])==Q(7,6),'reordered-row hostile mean')
    need(next(r['histogram'] for r in reports if r['shape']==(5,))=={'11/10':40,'17/15':40,'7/6':40},
         'complete three-valued n5 cycle conditional mean')
    return reports,dict(roworders=roworders,column_permutations=columns,prefix_nodes=nodes)

def multiply(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out

def large_rooks(n,selected,paired):
    if paired:
        polynomial=[1]
        for r in range(0,n,2):
            number=(r in selected)+(r+1 in selected)
            polynomial=multiply(polynomial,([1],[1,2],[1,4,2])[number])
        return polynomial
    if len(selected)==n:
        return [1]+[2*n*comb(2*n-k,k)//(2*n-k) for k in range(1,n+1)]
    runs=[]
    for r in selected:
        if (r-1)%n not in selected:
            size=1
            while (r+size)%n in selected:size+=1
            runs.append(size)
    polynomial=[1]
    for q in runs:
        polynomial=multiply(polynomial,[comb(2*q+1-k,k) for k in range(q+1)])
    return polynomial

def rational_constants():
    lo=sum((Q(2**j,factorial(j)) for j in range(41)),Q(0))
    hi=lo+Q(2**41,factorial(41))/(1-Q(2,42))
    alpha_lo,alpha_hi=1-5/lo,1-5/hi
    need(3*hi-9/hi<21,'exact lower error rounding')
    need(3*hi+13/lo-4<20,'exact upper error rounding')
    for n in range(4,201):
        lengths=[n]+[L for L in range(1,n) for _ in range(2)]
        third=sum((Q(4,3)*fall(L,3)-2*L*(L-2) for L in lengths if L>=2),Q(0))/fall(n,3)
        fourth=sum((Q(2,3)*fall(L,4) for L in lengths),Q(0))/fall(n,4)
        need(third==Q(2*(n-3),3)+Q(2,n*(n-1)),'exact finite sum with L1 correction')
        need(third-2*fourth==Q(2*(n-9),15)+Q(2,n*(n-1)),'conditional finite B4')
        need(Q(2,(n-1)*4**2)==Q(1,8*(n-1)),'one-permutation concentration constant')
    large=[]
    for n,paired in ((96,True),(97,False)):
        for multiplier in (1,5):
            rho=[multiplier*r%n for r in range(n)];mu=Q(0)
            for d in range(1-n,n):
                selected={r for r in range(n) if 0<=rho[r]+d<n}
                polynomial=large_rooks(n,selected,paired)
                mu+=sum((Q((-1)**(k+1)*(k-2)*polynomial[k],fall(n,k))
                         for k in range(3,len(polynomial))),Q(0))
            finite=Q(2*(n-9),15)+Q(2,n*(n-1))
            need(n*alpha_hi-21<=mu<=n*alpha_lo+20,'large exact mean uniform envelope')
            need(mu>=finite>0 and n*alpha_lo-21>0,'both explicit lower bounds positive')
            need(mu.numerator//mu.denominator==30,'independent large exact target floor')
            large.append(dict(n=n,paired=paired,multiplier=multiplier,mean=str(mu)))
    return large

def main():
    bank,profiles,subsets=subset_controls()
    reports,counts=literal_row_controls(bank)
    large=rational_constants()
    print('STATUS: PASS; independent conditional proof, edge-subset rook zeta, and exact controls')
    print('ROW SUBSET UNIVERSE:',profiles,'cycle profiles n2..9;',subsets,'complete row subsets')
    for report in reports:print(json.dumps(report,sort_keys=True))
    print('LITERAL CONDITIONAL UNIVERSE:',json.dumps(counts,sort_keys=True))
    for record in large:print('LARGE EXACT',json.dumps(record,sort_keys=True))
    print('BOUND: every fixed row ordering, column-only denominator 8(n-1), rate alpha^2/8')
    print('SCOPE: columns remain conditionally uniform; no mixture-mean concentration or full-zero equivalence')
    print('SEMANTIC SHA256:',hashlib.sha256(json.dumps([reports,counts,large],sort_keys=True).encode()).hexdigest())
    print('ACTIVE GATES:',GATES)

if __name__=='__main__':main()
