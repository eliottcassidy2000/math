"""Rational Burnside sums and SCC convolution, independently of graph generation."""
from fractions import Fraction
from collections import Counter
from math import gcd, factorial
from rooted_growth import need,save
import sys

def odd_partitions(n,minimum=1):
    if n==0:yield ()
    for k in range(minimum,n+1,2):
        for tail in odd_partitions(n-k,k):yield (k,)+tail

def census(N):
    a=[1];r=[0];b=[0];sroot=[0];rows=[]
    for n in range(1,N+1):
        aa=Fraction(0);rr=Fraction(0)
        for parts in odd_partitions(n):
            mult=Counter(parts);z=1
            for k,m in mult.items():z*=k**m*factorial(m)
            exponent=sum((k-1)//2 for k in parts)+sum(gcd(p,q) for i,p in enumerate(parts) for q in parts[i+1:])
            term=Fraction(1<<exponent,z);aa+=term;rr+=mult[1]*term
        need(aa.denominator==rr.denominator==1,'Burnside integrality')
        a.append(int(aa));r.append(int(rr));b.append(r[n]-sum(a[i]*b[n-i] for i in range(1,n)))
        sroot.append(b[n]-sum(a[i]*sroot[n-i] for i in range(1,n)))
        need(b[n]>=0 and sroot[n]>=0,'SCC positivity')
        rows.append(dict(n=n,ordinary=a[n],all_rooted=r[n],end_rooted=b[n],strong_rooted=sroot[n]))
    return rows

def main(N=30):
    rows=census(N)
    need([x['end_rooted'] for x in rows[:8]]==[1,1,2,7,33,236,2681,50709],'independent end-state census')
    need([x['ordinary'] for x in rows[:12]]==[1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168],'A000568 published control')
    save('cycle_index',dict(scope=f'Exact odd-partition Burnside and formal SCC division through{N}',rows=rows))
    print(rows[:10])
if __name__=='__main__':main(int(sys.argv[1]) if len(sys.argv)>1 else 30)
