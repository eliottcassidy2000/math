from itertools import combinations,product
from sevenfork_profile import *

T=tagged(CENTER,TAG)
EP=[]
for i,j in combinations(range(8),2):
 r=tagged(toggle(CENTER,(i,j)),TAG) if j<7 else tagged(CENTER,TAG^(1<<i))
 EP.append(((i,j),tuple(a-b for a,b in zip(r,T))))

def delta(e,x,y):
 a,q,c,t,e,r,k=e
 return a*(f(x)+f(y))+q*x*y+c*y+t*x+e*(y&1)+r*(x&1)+k

def verify():
 finite=set();rows=[]
 tags=[d for _,d in EP if d[1]==0]
 non=[(ij,d) for ij,d in EP if d[1]!=0]
 # Each non-tag edge is covered by bounded x or bounded y, including axes.
 for edge,d in non:
  a,q,c,t,e,r,k=d
  M=(2*abs(a)+abs(c)+abs(t)+abs(e)+abs(r)+abs(k)+2)//abs(q)
  for sw in (0,1):
   for x in range(M+1):
    s=q*x+(t if sw else c)
    C=a*f(x)+(c if sw else t)*x+(e if sw else r)*(x&1)+k
    E=r if sw else e
    if s*(s+a)<=0:
     assert not sw
     rows.append((edge,x))
     continue
    limit=(abs(C)+abs(E)+2)//min(abs(s),abs(s+a))
    for y in range(limit+1):
     X,Y=(y,x) if sw else (x,y)
     if abs(delta(d,X,Y))==2:finite.add((X,Y))
 # If only tag arcs are used, the two unit equations have no binary-height term.
 for d,g in combinations(tags,2):
  for px,py,j in product((0,1),(0,1),(-2,2)):
   _,_,c,t,e,r,k=d;_,_,C,T,E,R,K=g
   u=j-e*py-r*px-k;v=-j-E*py-R*px-K
   det=t*C-T*c
   assert det
   X=u*C-v*c;Y=t*v-T*u
   if X%det or Y%det:continue
   x,y=X//det,Y//det
   if min(x,y)>=0 and x%2==px and y%2==py:finite.add((x,y))
 # Pair any remaining unbounded row with every other arc.
 for (edge,x) in rows:
  d=dict(EP)[edge]
  for edge2,g in EP:
   if edge2==edge:continue
   for py,j in product((0,1),(-2,2)):
    def line(p,j):
     a,q,c,t,e,r,k=p
     return (q*x+c,a,j-a*f(x)-t*x-e*py-r*(x&1)-k)
    a,b,u=line(d,j);c,d1,v=line(g,-j)
    det=a*d1-b*c
    if not det:
     # Parallel rows with a common solution can occur only for the designated old pair.
     if a*v==c*u and b*v==d1*u:
      assert x==0 and {edge,edge2}=={(0,2),(2,4)}
     continue
    Y=u*d1-v*b;F=a*v-c*u
    if Y%det or F%det:continue
    y,z=Y//det,F//det
    if y>=0 and y%2==py and z==f(y):finite.add((x,y))
 for x,y in finite:
  ds={delta(d,x,y) for _,d in EP}
  assert not ({-2,2}<=ds) or (x==0 and y>0 and y&(y+1)==0),(x,y)
 print("all-address classification verified;",len(finite),"finite candidate positions, exceptional rows",rows)
if __name__=='__main__':verify()
