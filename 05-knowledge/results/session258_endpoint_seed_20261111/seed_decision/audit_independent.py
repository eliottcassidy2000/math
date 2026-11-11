"""Independent permutation/split states and exhaustive saturated closure."""
from itertools import permutations
from collections import deque
import json
from pathlib import Path
import check as production
HERE=Path(__file__).resolve().parent

def require(x,m):
    if not x:raise RuntimeError(m)

def literal(T,O):
    n=len(T);v=[0]*9
    for p in permutations(range(n)):
        for k in range(n+1):
            if not all(T[p[j]]>>p[j+1]&1 for j in range(n-1) if j!=k-1):continue
            a=0 if k==0 else 1+((O>>p[k-1])&1)
            b=0 if k==n else 1+((O>>p[k])&1)
            v[a*3+b]+=1
    return tuple(v)

def transfer(left,right,cap):
    ans=[0]*9
    for i in range(9):
        le,ls=divmod(i,3)
        for j in range(9):
            re,rs=divmod(j,3)
            end=re if re else le
            start=ls if ls else rs
            ans[3*end+start]+=left[i]*right[j]
    return tuple(min(a,cap) for a in ans)

def main():
    bank,meta=production.make_bank(13)
    for L,b in bank:require(literal(b['rows'],b['selected'])==L,('literal_block',b))
    receipt=json.loads((HERE/'receipt.json').read_text())
    results=[]
    for q in range(1,14):
        alphabet=[L for L,b in bank if b['h']<=q+1]
        start=(1,0,0,0,0,0,0,0,0);seen={start};todo=deque([start]);accepted=[]
        while todo:
            x=todo.popleft()
            if x[1]+x[2]==x[8]+1 and x[2]+x[5]+x[8]==q:accepted.append(x)
            for L in alphabet:
                y=transfer(x,L,q+2)
                if y[1]+y[2]>q+1 or y[2]+y[5]+y[8]>q:continue
                if y not in seen:seen.add(y);todo.append(y)
        wanted=receipt['results'][q-1]
        require(bool(accepted)==(wanted['status']=='SEED'),('outcome',q))
        if accepted:
            direct=literal(tuple(wanted['rows']),wanted['selected_mask'])
            require(direct[1]+direct[2]==direct[8]+1 and direct[2]+direct[5]+direct[8]==q,('winner',q))
        results.append({'q':q,'status':'SEED' if accepted else 'NO_SEED','full_closure_states':len(seen),'accepting_states':len(accepted)})
    out={'status':'PASS','independent_literal_transfer_checks':len(bank),'results':results,'scope':'Fixed U=q in the selected-set transitive-tail constructor; not general ray mechanisms or eventual odd-part seeds.'}
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
