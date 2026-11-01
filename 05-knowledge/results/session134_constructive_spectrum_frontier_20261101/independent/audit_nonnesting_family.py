#!/usr/bin/env python3
"""Independent literal source and adjacency-DP checks of adjacent attained counts.
Uses only sibling auditor functions, never producer code.
"""
import json
from audit_migration import marked
from audit_migration_compiler import dp

def need(c,msg):
    if not c:raise RuntimeError(msg)
def main():
    A=[20,1,2,7,14];cut=17;n=len(A)
    adjacency=[[bool(A[i]>>j&1) for j in range(n)] for i in range(n)]
    source=marked(adjacency,cut)
    need(source==(9,33,32,2),'literal family source')
    rows=[]
    for k in range(6,13):
        r=k-5;size=n+r+1;x=size-1;out=[0]*size
        for i in range(n):
            out[i]=A[i]
            for j in range(n,n+r):out[i]|=1<<j
            if cut>>i&1:out[x]|=1<<i
            else:out[i]|=1<<x
        for i in range(n,n+r):
            out[i]|=1<<x
            for j in range(i+1,n+r):out[i]|=1<<j
        high=dp(out)[0][-1]
        low_size=k+2;low=[sum(1<<j for j in range(i+1,low_size)) for i in range(low_size)]
        low[0]^=1<<(low_size-1);low[-1]|=1
        lowH=dp(low)[0][-1]
        need(high==2**k+3 and lowH==2**k+1,'adjacent family counts')
        rows.append({'k':k,'lower_H':lowH,'upper_H':high,'lower_order':low_size,'upper_order':size})
    print(json.dumps({'status':'PASS','literal_source_h_F_q_g':source,'actual_family_pairs':rows,'scope':'Adjacency positives are independently verified; all-order construction follows native zero-tail recurrence. Non-nesting of any witnesses uses separately proved all-order +2 induced-gap classification.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
