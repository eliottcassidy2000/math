from collections import deque
import json
S=((0,0),(0,1),(1,0),(1,1),(2,0),(2,1),(2,2))
def monitor(state,pair):
 if pair==(0,0):return state
 if pair==(2,0):return None if state==2 else 1
 if pair==(2,2):return None if state==1 else 2
 return 0
out=[]
for b in(0,1):
 q=deque([(0,0,b)]);seen=set(q);edges=[]
 while q:
  im,om,c=q.popleft()
  for a in S:
   ni=monitor(im,a)
   if ni is None:continue
   u=abs(c-a[0]);v=abs(u-a[1]);pair=(u,v);no=monitor(om,pair)
   if pair not in S or no is None:raise RuntimeError((im,om,c,a,pair))
   dest=(ni,no,v);edges.append(((im,om,c),a,pair,dest))
   if dest not in seen:seen.add(dest);q.append(dest)
 if len(seen)!=7 or len(edges)!=44:raise RuntimeError('unexpected graph')
 out.append(dict(boundary=b,states=sorted(seen),edges=edges))
print(json.dumps(dict(status='PASS',claim='All finite pair words and hence every finite prefix of an infinite word',cases=out),indent=2))
