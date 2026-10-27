"""Independent exact terminal-state audit. No production imports."""
from fractions import Fraction as Q
from math import comb, factorial
import json

def require(ok, label):
    if not ok: raise RuntimeError(label)

def rank(matrix):
    basis={}
    for row in matrix:
        row=list(map(Q,row))
        for p,v in sorted(basis.items()):
            c=row[p]
            row=[x-c*y for x,y in zip(row,v)]
        pivot=next((i for i,x in enumerate(row) if x),None)
        if pivot is not None:
            c=row[pivot]; basis[pivot]=[x/c for x in row]
    return len(basis)

def h(d): return 1 if d==0 else (d//3+d%2)//2

def w(k,t): return Q(2**(k*(k-1)//2+k*t),factorial(k))

def transformed(k,t,d):
    # [x^(d+k)] (x/(1+x))^d F(2^t x/(1+x))
    total=Q(0)
    for j in range(k+1):
        power=d+j; delta=k-j
        coefficient=(int(delta==0) if power==0 else (-1)**delta*comb(power+delta-1,delta))
        total+=coefficient*w(j,t)
    return total

def main():
    maxd=60
    counts=[set() for _ in range(maxd+1)]; counts[0].add(0)
    for part in range(3,maxd+1,2):
        for d in range(part,maxd+1):
            counts[d].update(t+1 for t in counts[d-part])
    rankchecks=0
    for d in range(maxd+1):
        ts=sorted(counts[d])
        expected=[0] if d==0 else [t for t in range(1,d//3+1) if t%2==d%2]
        require(ts==expected and len(ts)==h(d),('support',d))
        for r in range(13):
            for shift in (None,0,d):
                matrix=[[w(k,t) if shift is None else transformed(k,t,shift) for t in ts] for k in range(r+1)]
                require(rank(matrix)==min(h(d),r+1),('rank',d,r,shift))
                rankchecks+=1
    witness=[sum(c*w(k,t) for c,t in zip([4,-5,1],[2,4,6])) for k in range(3)]
    require(witness==[0,0,2880], 'witness')
    for m in range(1,13):
        d=6*m; N=7*m-1
        ts=list(range(2,2*m+1,2))
        require(all(d>=3*t and (d-t)%2==0 for t in ts),'lowerbound attainable')
        require(rank([[w(k,t) for t in ts] for k in range(N-d+1)])==m,'lowerbound rank')
    table=[]
    for N in [0,1,2,19,31,100,200,1000,10000]:
        old=sum(h(d) for d in range(N+1))
        word=sum(min(h(d),1) if d==N else h(d) for d in range(N+1))
        new=sum(min(h(d),N-d+1) for d in range(N+1))
        table.append(dict(N=N,original_degree_cyclecount_states=old,residual_wordwise_states=word,blockwise_linear_dimension=new))
    print(json.dumps(dict(status='PASS',independent_rank_checks=rankchecks,support_degrees=[0,maxd],residual_budgets=[0,12],lower_bound_m=[1,12],witness=list(map(str,witness)),dimensions=table),indent=2))
if __name__=='__main__':main()
