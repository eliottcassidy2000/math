from time import monotonic

def autos(p):
 R={x*x%p for x in range(1,p)};full=(1<<p)-1
 rows=[sum(1<<y for y in range(p) if (y-x)%p in R) for x in range(p)]
 inc=[full^r^(1<<x) for x,r in enumerate(rows)]
 count=nodes=0;solutions=[]
 initial=[full]*p;initial[0]=1
 def rec(dom,mapping,used):
  nonlocal count,nodes
  nodes+=1
  left=[x for x in range(p) if mapping[x]<0]
  if not left:
   count+=1;solutions.append(tuple(mapping));return
  x=min(left,key=lambda x:(dom[x]&~used).bit_count());D=dom[x]&~used
  while D:
   bit=D&-D;D-=bit;y=bit.bit_length()-1
   nd=dom.copy();nm=mapping.copy();nm[x]=y;nu=used|bit;ok=True
   for v in left:
    if v==x:continue
    nd[v]&=rows[y] if rows[x]>>v&1 else inc[y]
    nd[v]&=~nu
    if not nd[v]:ok=False;break
   if ok:rec(nd,nm,nu)
 rec(initial,[-1]*p,0)
 # We fixed 0, so affine squares give exactly the expected candidate set.
 expected={tuple(a*x%p for x in range(p)) for a in R}
 if set(solutions)!=expected:raise ValueError('unexpected automorphism')
 return {'p':p,'fix0_automorphisms':count,'full_automorphisms':p*count,'nodes':nodes}
if __name__=='__main__':
 import json
 print(json.dumps([autos(p) for p in (11,31)],indent=2))
