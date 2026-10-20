#!/usr/bin/env python3
"""Certify the n=11 maximum witness using an explicit word enumeration."""
from itertools import permutations, product
from pathlib import Path
import json
from verify_extrema import adjacency, dp_h, frame_aut

HERE=Path(__file__).resolve().parent
k=5
external=7208
matching=2
a=adjacency(k,external)
n=2*k+1;f=2*k
rows=[sum(x<<j for j,x in enumerate(row)) for row in a]+[0]
for i in range(k):
    x,y=2*i,2*i+1
    rows[x]|=1<<f
    rows[f]|=1<<y
    if matching>>i&1:rows[x]|=1<<y
    else:rows[y]|=1<<x
rho=[i^1 for i in range(2*k)]+[f]
assert all(bool(rows[i]>>j&1)!=bool(rows[j]>>i&1)
           for i in range(n) for j in range(i+1,n))
assert all(bool(rows[i]>>j&1)==bool(rows[rho[j]]>>rho[i]&1)
           for i in range(n) for j in range(n) if i!=j)
words=[]
for p in permutations(range(k)):
    for signs in product(range(2),repeat=k-1):
        half=[2*i+s for i,s in zip(p,signs+(0,))]
        word=half+[f]+[rho[v] for v in reversed(half)]
        if all(rows[u]>>v&1 for u,v in zip(word,word[1:])):words.append(word)
assert len(words)==dp_h(k,external)==195
aut=frame_aut(k,external)
stab=[p for p in aut if all((matching>>i&1)==(matching>>p[i]&1)for i in range(k))]
assert len(aut)==3 and len(stab)==1
pairs=[(i,j) for i in range(k) for j in range(i+1,k)]
m=len(pairs)
images=set()
for p in permutations(range(k)):
    code=sum(a[2*p[i]][2*p[j]]<<z for z,(i,j) in enumerate(pairs))
    code+=sum(a[2*p[i]][2*p[j]+1]<<(m+z) for z,(i,j) in enumerate(pairs))
    images.add(code)
assert images==set(json.loads((HERE/"exact_k5.json").read_text())["max_masks"])
seen=set();completion_classes=[]
for mask in range(1<<k):
    if mask in seen:continue
    orbit={sum((mask>>i&1)<<p[i] for i in range(k)) for p in aut}
    seen.update(orbit)
    completion_classes.append({"matching_representative":mask,"orbit":sorted(orbit),
                               "blue":195*len(orbit)//len(aut),
                               "triangle_parity":(k-mask.bit_count())%2})
assert sum(r["blue"]==195 for r in completion_classes)==8
assert [sum(r["blue"]==195 and r["triangle_parity"]==p for r in completion_classes)
        for p in range(2)]==[4,4]
scores=[r.bit_count() for r in rows]
c3=n*(n-1)*(n-2)//6-sum(d*(d-1)//2 for d in scores)
data={"n":n,"k":k,"external_mask":external,"matching_mask":matching,
      "vertices":"x_i=2i, y_i=2i+1 for 0<=i<5; f=10",
      "adjacency_row_masks":rows,"rho":rho,
      "frame_automorphisms":aut,"matching_stabilizer":stab,
      "symmetric_HP":len(words),"blue":len(words)//len(stab),
      "scores":scores,"c3":c3,"symmetric_HP_words":words,
      "unique_maximizing_hollow_frame":True,"maximizing_frame_labelled_orbit":sorted(images),
      "completion_classes":completion_classes,"maximum_blue_class_multiplicity":8}
(HERE/"n11_witness.json").write_text(json.dumps(data,indent=2)+"\n")
print(json.dumps({key:val for key,val in data.items() if key!="symmetric_HP_words"},indent=2))
