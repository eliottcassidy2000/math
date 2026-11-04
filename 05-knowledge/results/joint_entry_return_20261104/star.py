"""Actual two-complement exposure of the inherited bridge--sink star."""
from pathlib import Path
from itertools import permutations
import json
from repair import native,need,strong
from next_returns import gap
HERE=Path(__file__).resolve().parent

def action(T,coeff):
    out={};n=len(T)
    for p in permutations(range(n)):
        for k in range(1,n):
            c=coeff.get((p[k-1],p[-1]),0)
            if c and all(T[a]>>b&1 for a,b in zip(p[:k],p[1:k])) and all(T[a]>>b&1 for a,b in zip(p[k:],p[k+1:])):
                row=(p[0],p[k]);out[row]=out.get(row,0)+c
    return [[*k,v] for k,v in sorted(out.items()) if v]

def main():
    receipts=[]
    for m in (1,2,4,9,20):
        D=(2,4,9,3)
        for _ in range(m):D=native.extend(D,{0,1,2})
        n=len(D);w,z,a,b,c=3,4,n,n+1,n+2
        X=D
        for parent in (w,a,b):X=native.extend(X,set(range(len(X)))-{parent})
        cycle=(w,)+tuple(range(4,n))+(0,1,2)
        P=(b,c,a,2,w)+cycle[1:-1]
        U,V,W,p,q,case=native.traced_gadget(X,c,a,P)
        Y=native.restrict(W,n+3);need(Y==native.flip(X,c,a) and strong(Y),'Actual paid target')
        coeff={(z,a):1,(a,z):-1,(w,a):-1,(a,w):1};fb=(0,w,a)
        vals=[]
        for T in (X,Y):
            terms=[gap.count(T,fb,2,starts=((0,z),(1,b)),ends=((0,u),(1,v))) for u,v in coeff]
            vals.append(dict(terms=terms,action=sum(v*k for v,k in zip(terms,coeff.values()))))
        need(vals[0]['action']==0 and vals[1]['action']==1,'Uniform star exposed with unit amplitude')
        literal=None
        if m==1:
            before=action(X,coeff);after=action(Y,coeff)
            need(not before and [z,b,1] in after,'Independent full eight-vertex permutation/cut action')
            literal=dict(before=before,after=after)
        receipts.append(dict(m=m,source_rows=list(X),target_rows=list(Y),star_roles=dict(z=z,w=w,a=a),
                             test_row=[z,b],actions=vals,helper_out=U[-1],paths=[p,q],literal=literal))
    (HERE/'star.json').write_text(json.dumps(dict(status='FINITE_EXACT_NATIVE_EXPOSURE',receipts=receipts),indent=2)+'\n')
    print('Native star exposure at module sizes',[r['m'] for r in receipts])

if __name__=='__main__':main()
