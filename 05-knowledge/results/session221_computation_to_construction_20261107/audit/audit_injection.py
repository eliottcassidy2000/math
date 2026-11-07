#!/usr/bin/env python3
from itertools import product
from pathlib import Path
import json
from audit_physical import graph
ROOT=Path(__file__).resolve().parent

def paths(out):
    n=len(out);full=(1<<n)-1;result=set()
    def walk(path,mask):
        if mask==full:result.add(path);return
        rem=out[path[-1]]&~mask
        while rem:
            b=rem&-rem;rem-=b;walk(path+(b.bit_length()-1,),mask|b)
    for v in range(n):walk((v,),1<<v)
    return result

def main():
    count=0
    for n in range(2,8):
        for s in product((0,1),repeat=n):
            before=graph(s);old=paths(before)
            for i in range(n-1):
                if s[i:i+2]!=(1,0):continue
                after=before[:];after[i]^=1<<(i+1);after[i+1]^=1<<i;new=paths(after);images=set()
                for p in old:
                    pos=p.index(i)
                    if pos<len(p)-1 and p[pos+1]==i+1:q=p[:pos]+(i+1,i)+p[pos+2:]
                    else:q=p
                    if q not in new:raise ValueError('not a legal image')
                    images.add(q)
                if len(images)!=len(old):raise ValueError('not injective')
                first={p for p in new if any(p[j:j+3]==(n,i+1,i) for j in range(n-1))}
                second={p for p in new if any(p[j:j+3]==(i+1,i,n) for j in range(n-1))}
                if first&second or first|second!=new-images:raise ValueError('unmatched classification')
                a=(1 if i==0 else s[i-1])+sum(s[j]*2**(j-i-2) for j in range(i+2,n))
                b=(1 if i==n-2 else 1-s[i+2])+sum((1-s[j])*2**(i-j-1) for j in range(i))
                if (a,b)!=(len(first),len(second)):raise ValueError('unmatched counts')
                count+=1
    r={'status':'PASS','independent_exact_path_set_injections':count,'maximum_chain_length':7,'checks':['legal injective image','exact two unmatched block languages','separate block cardinality formulas']}
    (ROOT/'injection_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
if __name__=='__main__':main()
