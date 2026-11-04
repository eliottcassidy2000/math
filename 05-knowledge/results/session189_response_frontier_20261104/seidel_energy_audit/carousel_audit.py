#!/usr/bin/env python3
"""Independent integer switching/second-moment audit of K_(0,(n-1)/2)."""
from pathlib import Path
import json,hashlib
import audit as A
P=Path(__file__).resolve().parent

def carousel(n):return tuple(sum(1<<((i+d)%n) for d in range(1,(n+1)//2)) for i in range(n))
def main():
    rows=[];subset_count=0;switch_pairs=0;inverse_entries=0
    for n in range(3,20,2):
        d=(n-1)//2;q=carousel(n);order=[v for i in range(d) for v in [i,i+d+1]]+[d];ell=n-2
        A.require(len(set(order))==n and order[0]==0 and order[-1]==d,'physical switch order')
        for i in range(n):
            for j in range(i+1,n):
                A.require((1 if q[order[i]]>>order[j]&1 else -1)*(-1)**(i+j)==1,'switch transitive matrix');switch_pairs+=1
        # Explicit inverse identity for every even size through n-1.
        for k in range(2,n,2):
            T=[[0 if i==j else 1 if i<j else -1 for j in range(k)] for i in range(k)]
            V=[[0 if i==j else (-1)**(i+j) if i<j else -(-1)**(i+j) for j in range(k)] for i in range(k)]
            for i in range(k):
                for j in range(k):A.require(sum(T[i][h]*V[h][j] for h in range(k))==int(i==j),'transitive inverse');inverse_entries+=1
        moments=[0]*4;reflection=True
        for mask in range(1<<ell):
            if mask.bit_count()%2:continue
            ids=[i+1 for i in range(ell) if mask>>i&1];u=sum((-1)**(rank+i) for rank,i in enumerate(ids,1));rev=sorted(ell+1-i for i in ids);ur=sum((-1)**(rank+i) for rank,i in enumerate(rev,1));A.require(ur==-u,'odd moments cancel by reflection');subset_count+=1
            for p in range(4):moments[p]+=u**p
        A.require(moments==[2**(n-3),0,(n-3)*2**(n-3),0],'all-order moment formula controls')
        numerator=moments[3]-moments[2]-moments[1]+moments[0];coefficient=A.F(64,2**n)*numerator;A.require(coefficient==-8*(n-4),'shared cubic closed form')
        # Independently use exact full induced-subset inverses through nine names.
        inverse_cubic=None
        if n<=9:
            ev,od,_,_,_=A.banks(q,False);ec,oc,_,_,_=A.banks(A.opposite(q),False);inverse_cubic=ev[1|(1<<d),1]+od[1|(1<<d),1]-ec[1|(1<<d),1]-oc[1|(1<<d),1];A.require(inverse_cubic==coefficient,'exact rational named coefficient')
        rows.append({'n':n,'physical_rows':q,'switch_order':order,'moments':moments,'coefficient':int(coefficient),'inverse_coefficient':None if inverse_cubic is None else int(inverse_cubic)})
    # Explicit TT4 normalized profile, including every nonzero coordinate.
    counts=[1,0,0,0,0]
    for _ in range(4):counts=[0]+[k*(counts[k]+counts[k-1]) for k in range(1,5)]
    A.require(counts[1:]==[1,14,36,24],'TT4 profile')
    for k,v in enumerate(counts[1:],1):A.require(v%16==pow(-2,k-1,16),'TT4 geometric profile')
    (P/'carousel_controls.json').write_text(json.dumps(rows,sort_keys=True,indent=2)+'\n');receipt={'status':'PASS','odd_carousels':len(rows),'even_subsets_enumerated':subset_count,'physical_switch_pairs':switch_pairs,'inverse_identity_entries':inverse_entries,'exact_rational_bank_crosschecks':4,'sources':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['carousel_audit.py','audit.py']},'artifacts':{'carousel_controls.json':hashlib.sha256((P/'carousel_controls.json').read_bytes()).hexdigest()}}
    (P/'carousel_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
