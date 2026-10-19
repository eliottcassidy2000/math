#!/usr/bin/env python3
"""Independent exact controls for the appended odd-order recovery section."""
from collections import Counter
from itertools import combinations_with_replacement, product
from pathlib import Path
import hashlib
import json

checks=0

def need(ok,context):
    global checks
    checks+=1
    if not ok:
        raise RuntimeError(context)

rows=[]
for moduli in [(1,),(3,),(5,),(7,),(9,),(15,),(21,),(25,),(3,3),(5,5),(3,3,3)]:
    elements=list(product(*(range(m) for m in moduli)))
    index={x:i for i,x in enumerate(elements)}
    zero=(0,)*len(moduli)
    add=[[index[tuple((a+b)%m for a,b,m in zip(x,y,moduli))] for y in elements] for x in elements]
    sub=[[index[tuple((a-b)%m for a,b,m in zip(x,y,moduli))] for y in elements] for x in elements]
    cases=0
    for k in range(4):
        for ts in combinations_with_replacement(range(len(elements)),k):
            c=[0]*len(elements)
            total=0
            for t in ts:
                total=add[total][t]
            for mask in range(1<<k):
                v=0
                for i,t in enumerate(ts):
                    if mask>>i&1:
                        v=add[v][t]
                c[v]+=1
            periods=[a for a in range(len(elements)) if all(c[x]==c[add[x][a]] for x in range(len(elements)))]
            centers=[a for a in range(len(elements)) if all(c[x]==c[sub[a][x]] for x in range(len(elements)))]
            need(periods==[0],('odd-order period',moduli,ts))
            need(centers==[total],('odd-order unique center',moduli,ts))
            cases+=1
    rows.append({'cyclic_factor_orders':moduli,'group_order':len(elements),'multisets_length_at_most_3':cases})

# Explicit q=7 recovery in the multiplicative signed quotient.
q=7
classes=[1,2,3]
canon=lambda x:min(x%q,q-x%q)
prime_rows=[]
for ps in [[2,3,5,11],[2,3,5,17]]:
    n=1
    ds=[1]
    for p in ps:
        n*=p
        ds=ds+[p*d for d in ds]
    c=Counter(canon(d) for d in ds)
    centers=[a for a in classes if all(c[x]==c[canon(a*pow(x,-1,q))] for x in classes)]
    need(centers==[canon(n)]==[1],('signed product unique center',ps))
    prime_rows.append({'source':n,'kernel':[c[x] for x in classes],
                       'signed_product':canon(n),'unit_product':n%q,
                       'negative_product_legendre':1 if pow(-n%q,3,q)==1 else -1})
need(prime_rows[0]['kernel']==prime_rows[1]['kernel']==[6,5,5],'same nonuniform signed kernel')
need([r['unit_product'] for r in prime_rows]==[1,6],'opposite full unit lifts')
need([r['negative_product_legendre'] for r in prime_rows]==[-1,1],'opposite candidate phases')
for q in [3,7,11,19,23,31,43]:
    need(pow(q-1,(q-1)//2,q)==q-1,('Legendre does not descend',q))

out={'status':'PASS','scope':'Odd-order finite abelian groups and q=7 sign-loss hostile',
     'group_controls':rows,'total_multisets':sum(r['multisets_length_at_most_3'] for r in rows),
     'q7_actual_sources':prime_rows,'always_active_checks':checks}
payload=json.dumps(out,sort_keys=True,indent=2)+'\n'
Path(__file__).with_name('odd_order_receipt.json').write_text(payload,encoding='utf-8')
print('PASS:',checks,'always-active checks;',out['total_multisets'],'odd-order group multisets')
print('receipt sha256:',hashlib.sha256(payload.encode()).hexdigest())
