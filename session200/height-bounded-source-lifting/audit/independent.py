import json,itertools
from math import isqrt
T=(1,0,2,1)
def need(q,m):
 if not q:raise RuntimeError(m)
def move(s,a):
 D,E,U,d=s
 return {'I':s,'O':(D,0,2*U,d),'Z':(D-E+U,d,2*U-d,d),'J':(D-E+U,D-E,D+U,D)}[a]
def K(v,w):
 D,E,U,d=v;A,B,V,b=w
 return D*b+d*A+E*B+U*V-E*b-d*B-U*b-d*V+d*b
def bank(B):
 todo=[T];seen={T}
 for s in todo:
  for a in 'OZJ':
   q=move(s,a)
   if q[2]<=B and q not in seen:seen.add(q);todo.append(q)
 return seen
def word(w):
 s=T
 for a in w:s=move(s,a)
 return s
for n in range(3,201):
 B=isqrt(4*(n-1));ell=B.bit_length();R=2**ell-1
 need(R<=2*B-1 and B<=n and 2*R<4*n,('range',n))
 for s in bank(B):
  D,E,U,d=s;G=D+2*U
  for a in 'OZJ':
   A,C,V,b=move(s,a)
   need(2*(A+2*V)>=3*G,('length growth',s,a))
mods=[3,5,7,13];pairs=[('',''),('OZ',''),('OZ','ZJ'),('OZ','OO')]
values=[K(word(a),word(b)) for a,b in pairs]
need(values==[3,11,35,47],'hostile values')
need(all((v-21)%p==0 for v,p in zip(values,mods)),'local solubility')
image=sorted({K(s,t) for s in bank(8) for t in bank(8)})
need(21 not in image and 1365>64,'global absence')
need(910**2%1365==910 and 456**2%1365==456 and 910*456%1365==0 and (910+456)%1365==1,'CRT idempotents')
need(K(word('OZO'),word('OOO'))==201 and all(K(s,t)!=67 for s in bank(16) for t in bank(16)),'terminal-only wrap')
print(json.dumps({'status':'PASS','range_growth_targets_through':200,'crt_pair_values':values,'crt_modulus':1365,'complete_C1_image_cap8':image,'idempotents':[910,456],'terminal_wrap_actual_count':201},sort_keys=True,indent=2))
