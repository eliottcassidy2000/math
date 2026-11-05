from collections import defaultdict
import heapq,json

def worddata(v):
 w=bin(v)[2:];h=1;a=0;p=0;c=0;q=1
 for b in map(int,w):
  if b==0:h+=a+p;c+=q
  a=2*a+b;p=b;q*=2
 return h,a,c,q,1-p

def census(X):
 found={(1,1,1):None};heap=[(1,1,1)]
 while heap:
  h,s,l=st=heapq.heappop(heap)
  children=[]
  if h>1:children.append(((h+l-s,h-s,h+l),0))
  for v in range(1,X//l):
   g,a,c,q,e=worddata(v)
   f=g*h-e*s+(a+1)*l
   if f<=X:children.append(((f,h,(1+c)*h+q*l),v))
  for ns,v in children:
   if ns[0]<=X and ns not in found:
    if ns[0]<=h:raise Exception('nonincreasing')
    found[ns]=(st,v);heapq.heappush(heap,ns)
 return found
if __name__=='__main__':
 import sys
 X=int(sys.argv[1]) if len(sys.argv)>1 else 1000
 F=census(X);counts=sorted(set(s[0] for s in F));missing=sorted(set(range(1,X+1,2))-set(counts))
 print('X',X,'states',len(F),'values',len(counts),'missing',missing[:100]);
 print('states13',[s for s in F if s[0]==13]);print('counts',counts[:100])
 json.dump({'X':X,'states':[{'state':s,'parent':p} for s,p in sorted(F.items())],'values':counts,'missing':missing},open('census.json','w'))
