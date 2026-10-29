"""Port-preserving interval transfer for the inherited tribonacci tournament."""
import argparse
from collections import defaultdict
from itertools import combinations, permutations, product
from math import comb, gcd
import json
from random import Random

from atom_precision import endpoints, determinant
from complement_kernel import HERE, need
from three_run import native
from two_run import two_runs
from paley_pencil import diagonal


def tribonacci_rows(n):
    return [(1 << (i+1) if i+1<n else 0) | ((1 << max(0,i-1))-1) for i in range(n)]


def rational_series(numerator, limit):
    values=[]
    for n in range(limit+1):
        values.append((numerator[n] if n<len(numerator) else 0)
                      +sum(values[n-d] for d in (1,2,3) if n>=d))
    return values


def endpoint_statistics(limit):
    h=rational_series([1,0,-1],limit)
    a=rational_series([1,-1],limit)
    b=rational_series([0,0,1],limit)
    ph,pa=[0],[0]
    for x,y in zip(h,a): ph.append(ph[-1]+x);pa.append(pa[-1]+y)
    return h,a,b,ph,pa


def endpoint_entry(n,s,t,stats):
    h,a,b,_,_=stats
    if not(0<=s<n and 0<=t<n): return 0
    if s==0 and t==n-1: return 1
    if s<=t: return 0
    flags=int(t==0)+int(s==n-1)
    return (h if flags==0 else a if flags==1 else b)[s-t-1]


def start_count(n,s,stats):
    _,a,b,ph,pa=stats
    if not 0<=s<n: return 0
    if s==0: return 1
    return b[s-1]+pa[s-1] if s==n-1 else a[s-1]+ph[s-1]


def endpoint_bank(limit):
    stats=endpoint_statistics(limit)
    h=stats[0]
    matrices=[[]]
    for n in range(1,limit+1):
        E=[[endpoint_entry(n,s,t,stats) for t in range(n)] for s in range(n)]
        matrices.append(E)
        need(sum(map(sum,E))==h[n], 'inherited Hamilton mass from exact endpoint formula')
        need(list(map(sum,E))==[start_count(n,s,stats) for s in range(n)], 'constant-time endpoint marginal')
    return h,matrices


def color_counts(j,limit):
    f=[]
    for n in range(limit+1):
        value=1 if n==0 else -1 if n==2 else 0
        value+=(j*f[n-1] if n>=1 else 0)+(f[n-2] if n>=2 else 0)+(j*f[n-3] if n>=3 else 0)
        f.append(value)
    return f


def cover_counts(k,limit):
    values=[0]*(limit+1)
    for j in range(k+1):
        for n,v in enumerate(color_counts(j,limit)):
            values[n]+=(-1)**(k-j)*comb(k,j)*v
    return values


def recurrence_polynomial(k):
    # Low degree first: product (X^3-jX^2-X-j).
    polynomial=[1]
    for j in range(1,k+1):
        factor=[-j,-1,-j,1]
        updated=[0]*(len(polynomial)+3)
        for a,x in enumerate(polynomial):
            for b,y in enumerate(factor): updated[a+b]+=x*y
        polynomial=updated
    return polynomial


def joint_entry(n,starts,ends,stats=None):
    """One actual M_k entry without building its (n)_k-square response matrix."""
    k=len(starts)
    need(len(ends)==k and len(set(starts))==len(set(ends))==k, 'distinct ordered ports')
    need(all(0<=v<n for v in starts+ends), 'ports within source')
    stats=endpoint_statistics(n) if stats is None else stats
    h=stats[0]
    powers=[3**i for i in range(k)]
    # Three local coordinates: empty e, terminal-only G, joint F.
    states=[defaultdict(int) for _ in range(n+1)]
    states[0][(-1,0)]=1
    for left in range(n):
        for right in range(left+1,n+1):
            length=right-left
            for color in range(k):
                s,t=starts[color]-left,ends[color]-left
                start=start_count(length,s,stats)
                end=start_count(length,length-1-t,stats)
                both=endpoint_entry(length,s,t,stats)
                for (last,state),weight in states[left].items():
                    if color==last: continue
                    digit=state//powers[color]%3
                    if digit==0:
                        if end: states[right][(color,state+powers[color])]+=weight*end
                        if both: states[right][(color,state+2*powers[color])]+=weight*both
                    elif digit==1:
                        states[right][(color,state)]+=weight*h[length]
                        if start: states[right][(color,state+powers[color])]+=weight*start
                    # Old F is overwritten, not multiplied by an independent start law.
    return sum(states[n].get((color,3**k-1),0) for color in range(k))


def brute_cover_masses(n):
    rows=tribonacci_rows(n)
    result=[0]*(n+1)
    for P in permutations(range(n)):
        required=sum(not(rows[u]>>v&1) for u,v in zip(P,P[1:]))
        for k in range(required+1,n+1):
            result[k]+=comb(n-1-required,k-1-required)
    return result


def run():
    limit=100
    h,bank=endpoint_bank(limit)
    stats=endpoint_statistics(limit)
    endpoint_checks=[]
    for n in range(1,12):
        actual=endpoints(tribonacci_rows(n))
        need(actual==bank[n], 'independent subset DP endpoint matrix')
        if n>=3: need(abs(determinant(bank[n]))==1, 'all-characteristic one-run recovery control')
        endpoint_checks.append(dict(n=n,mass=h[n]))
    counts={str(k):cover_counts(k,limit) for k in range(1,7)}
    polynomials={str(k):recurrence_polynomial(k) for k in range(1,7)}
    for k in range(1,7):
        values,poly=counts[str(k)],polynomials[str(k)]
        need(all(sum(c*values[n+i] for i,c in enumerate(poly))==0
                 for n in range(1,limit-len(poly)+2)), 'order3k recurrence after the constant polynomial part')
    for n in range(1,9):
        direct=brute_cover_masses(n)
        need(all(direct[k]==cover_counts(k,n)[n] for k in range(1,n+1)), 'all arities by permutation/cut enumeration')
    rng=Random(20261029)
    query_records=[]
    for n in range(3,10):
        _,ports,M=native(tribonacci_rows(n))
        pairs=list(product(range(len(ports)),repeat=2)) if n<=5 else [(rng.randrange(len(ports)),rng.randrange(len(ports))) for _ in range(96)]
        for i,j in pairs:
            value=joint_entry(n,ports[i],ports[j],stats)
            need(value==M[i][j], 'joint port transfer versus both native enumerations')
            if value or len(query_records)<10: query_records.append(dict(n=n,start=ports[i],end=ports[j],value=value))
        need(sum(map(sum,M))==counts['3'][n], 'joint matrix mass equals scalar recurrence')
    for n in range(2,8):
        ports,M,_=two_runs(tribonacci_rows(n))
        for _ in range(32):
            i,j=rng.randrange(len(ports)),rng.randrange(len(ports))
            need(joint_entry(n,ports[i],ports[j],stats)==M[i][j], 'independent two-run port checks')
    # Long queries permit many colour intervals, with all six endpoint names retained.
    long=[]
    for n in (20,40,100):
        starts=tuple(range(n-3,n));ends=(0,1,2)
        value=joint_entry(n,starts,ends,stats)
        need(value>0,'nonempty specified long port consumer')
        long.append(dict(n=n,starts=starts,ends=ends,value=value))
    clock=[]
    clock_images=[]
    from sympy import Matrix
    for j in range(1,13):
        C=Matrix([[j,1,j],[1,0,0],[0,1,0]])
        A=Matrix.eye(3)
        for r in range(13):
            need(diagonal(A.tolist())==[1,1,j**r], 'all-prime clock power Smith formula')
            clock.append(dict(j=j,r=r,smith=[1,1,j**r]))
            if j<=6 and 1<=r<=4:
                rows=A.tolist()
                for modulus in range(2,9):
                    image={tuple(sum(a*x for a,x in zip(row,v))%modulus for row in rows)
                           for v in product(range(modulus),repeat=3)}
                    expected=modulus**3//gcd(modulus,j**r)
                    need(len(image)==expected,'independent full finite clock image')
                    clock_images.append(dict(j=j,r=r,modulus=modulus,size=expected))
            A=A*C
    from complement_kernel import complement
    from complement_verify import strong
    source=tribonacci_rows(6)
    path=(2,3,4,5,0,1)
    target=complement(source,path)
    def profile(rows):
        degrees=[r.bit_count() for r in rows]
        return sorted((degrees[i],sum(degrees[j] for j in range(len(rows)) if rows[i]>>j&1))
                      for i in range(len(rows)))
    need(strong(source) and strong(target) and profile(source)!=profile(target), 'strong actual complement exits interval family')
    need(sum(map(sum,endpoints(target)))==h[6]==17,'equal scalar counts do not retain the target relation')
    return dict(endpoint_checks=endpoint_checks,cover_counts=counts,recurrence_polynomials=polynomials,
                joint_queries=query_records,long_queries=long,clock_checks=clock,clock_images=clock_images,
                nonclosure=dict(source=source,path=path,target=target,profiles=[profile(source),profile(target)]))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=json.loads(json.dumps(run()))
    path=HERE/'run_recurrence.json'
    if args.write: path.write_text(json.dumps(result,separators=(',',':'))+'\n')
    else: need(result==json.loads(path.read_text()),'unchanged receipt')
    print('PASS endpoint sizes1..11, all cover arities through8, joint entries through9 and specified queries20/40/100')
