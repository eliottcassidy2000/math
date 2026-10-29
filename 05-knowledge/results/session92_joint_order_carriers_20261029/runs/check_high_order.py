"""Exact parity ensembles with identical proper marginals and different expected H."""
from fractions import Fraction
from itertools import product
from check_seam_carrier import hp,need
import json

def restore(t,mask):
    n=len(t);r=[((1<<n)-1)^(1<<v)^t[v] for v in range(n)]
    for v in range(n-1):
        if mask>>v&1:
            r[v]^=1<<(v+1);r[v+1]^=1<<v
    return tuple(r)
def main():
    report=[]
    for n in range(3,9):
        # Transitive and deterministic alternating-edge orientations.
        for specimen in range(2):
            t=[0]*n
            for i in range(n):
                for j in range(i+1,n):
                    u,v=(i,j) if specimen==0 or (i+j)%3 else (j,i)
                    t[u]|=1<<v
            k=n-1;values=[hp(restore(t,m)) for m in range(1<<k)]
            coefficient=sum((-1)**(k-m.bit_count())*v for m,v in enumerate(values))
            sign=1
            for i in range(k):sign*=1 if t[i]>>(i+1)&1 else -1
            need(coefficient==(1+(-1)**k)*sign,'top path coefficient')
            totals=[sum(v for m,v in enumerate(values) if m.bit_count()%2==parity) for parity in (0,1)]
            difference=Fraction(totals[0]-totals[1],1<<(k-1))
            need(difference==Fraction((-1)**k*coefficient,1<<(k-1)),'parity ensemble gap')
            # Every subset omitting one coordinate has identical projection counts;
            # summing implies all smaller coordinate marginals agree as well.
            for omitted in range(k):
                hist=[{},{}]
                for m in range(1<<k):
                    key=m&~(1<<omitted);d=hist[m.bit_count()%2]
                    d[key]=d.get(key,0)+1
                need(hist[0]==hist[1],'every maximal proper marginal agrees')
            report.append(dict(n=n,specimen=specimen,top_coefficient=coefficient,expectation_gap=str(difference)))
            if n%2==0:
                k=n-2
                vals=[hp(restore(t,m)) for m in range(1<<k)]
                c=sum((-1)**(k-m.bit_count())*v for m,v in enumerate(vals))
                sign=1
                for i in range(k):sign*=1 if t[i]>>(i+1)&1 else -1
                need(c==2*sign,'even-order path-plus-singleton repair')
                report[-1]['even_repaired_coefficient']=c
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
