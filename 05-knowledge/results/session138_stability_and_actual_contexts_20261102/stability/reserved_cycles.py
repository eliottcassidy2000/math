from bernstein import *

def triangles(n):
 edges=list(combinations(range(n),2));m=len(edges);idx={p:k for k,p in enumerate(edges)};pw=3**np.arange(m);t=np.zeros(3**m,dtype=np.int64)
 for i,j,k in combinations(range(n),3):
  for arcs in (((i,j),(j,k),(k,i)),((i,k),(k,j),(j,i))):
   terms={0:2**(n-3)}
   for x,y in arcs:
    e=idx[min(x,y),max(x,y)];s=1 if x<y else -1;new={}
    for exp,c in terms.items():
     new[exp]=new.get(exp,0)+c;new[exp+pw[e]]=new.get(exp+pw[e],0)+s*c
    terms=new
   for exp,c in terms.items():t[exp]+=c
 return t
if __name__=='__main__':
 n=int(sys.argv[1]);k=F(sys.argv[2]);h,u,v,V=polys(n);t=triangles(n);m=n*(n-1)//2
 for o in representatives(n):
  for s in (-1,1):
   a=(h-u+s*v-t)*k.denominator-k.numerator*V;b=bernstein(a,m,o)
   print(json.dumps({'n':n,'k':str(k),'orientation':o,'sign':s,'min':int(b.min()),'negatives':int((b<0).sum()),'zero':int((b==0).sum()),'scale':2**(n+m)*k.denominator,'argmin':int(b.argmin())}),flush=True)
