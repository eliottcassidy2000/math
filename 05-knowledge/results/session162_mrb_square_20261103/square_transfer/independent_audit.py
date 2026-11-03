from math import isqrt
from collections import Counter
p15=[8,1,15,10,6,3,13,12,4,5,11,14,2,7,9]
A=[3,22,14,2,23,13,12,4,21,15,10,6,19,17,8,1,24,25,11,5,20,16,9,7,18]
B=[18,7,9,16,20,5,11,25,24,1,8,17,19,6,10,15,21,4,12,13,3,22,14,2,23]
def sq(x): return isqrt(x)**2==x
def build(path,n):
 if len(path)!=n or sorted(path)!=list(range(1,n+1)):raise Exception('bad permutation')
 if not all(sq(u+v) for u,v in zip(path,path[1:])):raise Exception('bad path')
 prev=5; edges=[]; ds=[]
 for u,v in zip(path,path[1:]):
  d=abs(u-v);ds.append(d); s=1
  while s*s%2!=d%2 or s*s<2*(prev+1)+d:s+=1
  x=(s*s-d)//2;y=x+d
  if not (prev<x<y and y-x==d and sq(x+y)):raise Exception('bad embedding')
  edges.append((x,y));prev=y
 return ds,edges,prev
def pr(ds, edges, end):
 g=end-5-sum(y-x+1 for x,y in edges)
 if g<0:raise Exception('negative gap')
 P=R=1
 for d in ds:
  P*=1 if d==1 else 1+2**(d-1)
  R*=2**(d+1)+2*(3**(d-1)-1)
 return g,P,R*2**g
D,E,end=build(p15,15)
if end!=243 or E!=[(9,16),(25,39),(58,63),(70,74),(83,86),(93,103),(112,113),(124,132),(144,145),(159,165),(179,182),(194,206),(218,223),(241,243)]:raise Exception('path15 intervals')
g,P,R=pr(D,E,end)
if (g,P,R)!=(143,464858633025922867875,2286430711049262759319343471689019920235263576776074977652877326066121533571465216):raise Exception('path15 source counts')
da,ea,enda=build(A,25);db,eb,endb=build(B,25)
if A==B or A==B[::-1] or Counter(da)!=Counter(db) or set(ea)==set(eb) or (enda,endb)!=(846,623):raise Exception('collision inputs')
pa=pr(da,ea,846);pb=pr(db,eb,846)
if pa!=pb or pa[:2]!=(592,10918691908903704879820836527116367337011525929450009908515625):raise Exception('collision source counts')
print('PASS 15-path exact intervals, all square sums and arithmetic; 25-path equal P,R after padding and unequal physical edge sets')
