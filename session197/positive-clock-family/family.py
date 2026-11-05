from graphs import birth
Q=[14,32,2,22,103,13,47]
def graph(m,p,q,u,v):
 lengths=[m,p,u,1,q,1,v];starts=[sum(lengths[:i])for i in range(7)];rows=[0]*sum(lengths)
 for i in range(7):
  for a in range(lengths[i]):
   for j in range(7):
    for b in range(lengths[j]):
     if (i==j and a<b)or(i!=j and Q[i]>>j&1):rows[starts[i]+a]|=1<<(starts[j]+b)
 return rows

def coefficients(p,q,u,v):
 x,y,z,w=2**p,2**q,2**u,2**v
 return 2*z*w*(x-1)*(y-1),w*(y-1)*(4-x)+z*(x-1)*(4-y),(x-3)*(y-3)
def value(m,p,q,u,v):
 A,B,C=coefficients(p,q,u,v);return A*3**m+B*2**m+C

def returning(m,v):
 rows=graph(1,2,2,1,v);N=len(rows);ports=[1,2,3,4]
 for _ in range(m-1):rows,ports=birth(rows,ports,'R',15,4)
 mapping=[0]+list(range(N,N+m-1))+list(range(1,N))
 return rows,ports,mapping
