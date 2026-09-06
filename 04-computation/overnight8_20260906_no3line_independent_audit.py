#!/usr/bin/env python3
"""Independent density referee: unique boards and exact polynomial sums.

Standard library only. No producer, shore-label census, symbolic algebra,
or Poisson approximation is used to generate the finite controls.
"""
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb,factorial
from pathlib import Path
import hashlib,json,sys
sys.stdout.reconfigure(newline='\n')
CHECKS=0
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)
def falling(n,k):
    a=1
    for j in range(k):a*=n-j
    return a
def boards(n):
    counts=[0]*n;chosen=[]
    def rec(r):
        if r==n:
            if all(c==2 for c in counts):yield tuple(chosen)
            return
        for pair in combinations([j for j in range(n) if counts[j]<2],2):
            for j in pair:counts[j]+=1
            if all(2-c<=n-r-1 for c in counts):
                chosen.extend((r,j) for j in pair)
                yield from rec(r+1)
                del chosen[-2:]
            for j in pair:counts[j]-=1
    yield from rec(0)
def cycle_type(board,n):
    adj=[[] for _ in range(2*n)]
    for r,c in board:adj[r].append(n+c);adj[n+c].append(r)
    need(all(len(row)==2 for row in adj),'literal saturated simple board')
    seen=set();parts=[]
    for r in range(n):
        if r in seen:continue
        todo=[r];nr=0
        while todo:
            u=todo.pop()
            if u in seen:continue
            seen.add(u);nr+=u<n
            todo.extend(adj[u])
        parts.append(nr)
    return tuple(sorted(parts))
def unique_board_controls():
    expected={(3,):(6,F(1,3),F(2,9),F(2,3)),
              (4,):(72,F(1),F(49,36),F(25,72)),
              (2,2):(18,F(1),F(19,9),F(5,9)),
              (5,):(1440,F(5,3),F(461,120),F(101,360)),
              (2,3):(600,F(5,3),F(1969,450),F(47,200))}
    data={};witness=None
    for n in (3,4,5):
        accum=defaultdict(lambda:[Counter(),Counter(),F(0)])
        total=0
        for board in boards(n):
            total+=1;part=cycle_type(board,n);hist,xhist,joint=accum[part]
            diag=Counter(c-r for r,c in board)
            s=sum(comb(y,3) for y in diag.values() if y>=3)
            x=sum((b[0]-a[0])*(c[1]-a[1])==(b[1]-a[1])*(c[0]-a[0]) for a,b,c in combinations(board,3))
            need(s<=x,'actual integer collinearity target contains diagonal events')
            hist[s]+=1;xhist[x]+=1;accum[part][2]+=diag[0]*diag[1]
            if n==4 and s==0 and x and witness is None:witness=(board,x)
        for part,(hist,xhist,joint) in sorted(accum.items()):
            N=sum(hist.values());mean=F(sum(y*c for y,c in hist.items()),N)
            var=F(sum(y*y*c for y,c in hist.items()),N)-mean*mean
            zero=F(hist[0],N)
            need((N,mean,var,zero)==expected[part],'independent unique-board exact moments')
            multiplicities=Counter(part);aut=1
            for size,multiplicity in multiplicities.items():aut*=(2*size)**multiplicity*factorial(multiplicity)
            need(N*aut==factorial(n)**2,'constant orbit multiplicity matches uniform shore labels')
            L,M,R,C=n,n-1,n-1,n-1
            p2=F(4*n-6,n*(n-1)**2);q2=F(2,n*(n-1))
            need(joint/N==(L*M-R-C)*p2+(R+C)*q2,'exact path versus disjoint-pair normalization')
            if n==4:need(F(xhist[0],N)==(F(1,36) if part==(4,) else F(1,2)),'old full-event hostile retained')
            data[str(part)]={'boards':N,'mean':str(mean),'variance':str(var),'P_S0':str(zero),'P_X0':str(F(xhist[0],N))}
        need(total=={3:6,4:90,5:2040}[n],'complete unique board count')
    need(witness is not None,'S zero does not imply no-three-in-line')
    print('UNIQUE_BOARDS',json.dumps(data,sort_keys=True))
    print('S_ZERO_X_POSITIVE',json.dumps(witness))
    return data
def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c
def power(a,n):
    c=[F(1)]
    for _ in range(n):c=mul(c,a)
    return c
def add(a,b):return [(a[j] if j<len(a) else 0)+(b[j] if j<len(b) else 0) for j in range(max(len(a),len(b)))]
def scale(a,c):return [c*x for x in a]
def integrate01(a):return sum(c/F(j+1) for j,c in enumerate(a))
def interpolation(values,start=1):
    delta=list(map(F,values));out=[F(0)]*len(values);fall=[F(1)]
    for k in range(len(values)):
        for j,c in enumerate(fall):out[j]+=delta[0]*c/factorial(k)
        delta=[delta[j+1]-delta[j] for j in range(len(delta)-1)]
        fall=mul(fall,[-start-k,1])
    return out
def evalpoly(a,n):return sum(c*n**j for j,c in enumerate(a))
def discrete_cov_numerator(n):
    ds=range(1-n,n);L={d:n-abs(d) for d in ds}
    overlap=0
    for d in ds:
        for e in ds:
            if d==e:continue
            om=2*min(L[d],L[e]) if d*e>=0 else 2*max(L[d]+L[e]-n,0)
            overlap+=L[d]**2*L[e]**2*om
    return sum(x**3 for x in L.values())**2-sum(x**6 for x in L.values())-n*overlap
def discrete_local_numerator(n):
    L=[n-abs(d) for d in range(1-n,n)]
    return 24*sum(x**5 for x in L)+24*n*sum(x**4 for x in L)+4*n*n*sum(x**3 for x in L)
def coefficient_controls():
    same=F(2,4*7)
    # Integrate the opposite-sign triangular region using polynomial arithmetic.
    int_y=add(scale(mul([-1,1],add([1],scale(power([1,-1],3),-1))),F(1,3)),
              scale(add([1],scale(power([1,-1],4),-1)),F(1,4)))
    opposite=integrate01(mul([0,0,1],int_y))
    local=2*(F(8,6)+F(8,5)+F(4,12));cross=2-32*(same+opposite)
    need((same,opposite,local,cross)==(F(1,14),F(71,1260),F(98,15),F(-94,45)),'independent elementary integrals')
    need(local+cross==F(40,9) and (local+cross)/F(2,3)**2==10,'variance and Chebyshev leading constants')
    # Each sum is polynomial: split the square into integer triangles and use
    # power sums. Degrees <=8 (covariance numerator), <=6 (local numerator).
    cp=interpolation([discrete_cov_numerator(n) for n in range(1,10)])
    lp=interpolation([discrete_local_numerator(n) for n in range(1,8)])
    need(8*cp[-1]==cross and lp[-1]/3==local,'separate finite-sum leading coefficients')
    for n in range(10,36):
        need(evalpoly(cp,n)==discrete_cov_numerator(n),'fresh exact covariance-kernel sum identity')
        need(evalpoly(lp,n)==discrete_local_numerator(n),'fresh exact individual-kernel sum identity')
    print('COVARIANCE_NUMERATOR_POLYNOMIAL',list(map(str,cp)))
    print('LOCAL_NUMERATOR_POLYNOMIAL',list(map(str,lp)))
    print('COEFFICIENTS single98/15 cross-94/45 total40/9 density10')
def weighted_hostile():
    for n in range(3,10):
        W=[[int(j==i)*2+int(j==(i-1)%n)+int(j==(i+1)%n) for j in range(n)] for i in range(n)]
        total=sum(map(sum,W));rows=sum(sum(r)**2 for r in W)
        cols=sum(sum(W[i][j] for i in range(n))**2 for j in range(n))
        squares=sum(x*x for r in W for x in r)
        second=F(total*total-rows-cols+squares,n*(n-1))+F(squares-total,n)
        need(second==18-F(10,n-1),'weighted permutation factorial moment by entry intersections')
        if n==7:need(second==F(49,3)>16,'binary-union domination hostile')
    print('WEIGHTED_UNION independent matrix identity: E(Y)_2=18-10/(n-1), n7=49/3>16')
def main():
    data=unique_board_controls();coefficient_controls();weighted_hostile()
    primary=Path(__file__).with_name('overnight8_20260906_no3line_diagonal_density.py')
    need(hashlib.sha256(primary.read_bytes()).hexdigest()=='ead04c2f82c7c77a4f94c5b2ca379c7cc56dcdb37e252f4276a3826b2972c200','pinned primary producer')
    print('semantic_sha256',hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest())
    print('PASS',CHECKS,'explicit gates; no primary imports or shore-label enumeration')
if __name__=='__main__':main()
