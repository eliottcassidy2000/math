from collections import deque
import json,hashlib
from verify import hamilton,need

def transitions(st):
 h,x,y,d,m=st
 yield 'O',(h,2*x+d,d,d,m)
 yield 'Z',(h+x+y,2*x,0,d,m)
 if m<2:
  for j in range(4):
   yield str(j),(3*h+(9-3*j)*x+(3-j)*y+((3,5,5,3)[j]-3)*d,12*x+3*j*d,j*d,3*d,m+1)
def build(word):
 rows=[0]
 for c in word:
  old=len(rows)
  if c in 'OZ':
   rows.append(0)
   for i in range(old):
    if c=='Z' and i==0:rows[old]|=1
    else:rows[i]|=1<<old
  else:
   j=int(c);rows.extend([0]*3)
   for i in range(1,old):rows[i]|=7<<old
   for k in range(3):
    rows[old+k]|=1<<(old+(k+1)%3)
    if k<j:rows[0]|=1<<(old+k)
    else:rows[old+k]|=1
 return tuple(rows)
N=67;seed=(1,0,0,1,0);seen={seed:''};queue=deque([seed]);values={1:''}
while queue:
 st=queue.popleft();word=seen[st]
 for letter,z in transitions(st):
  w=word+letter
  if z[0]<=N:values.setdefault(z[0],w)
  if z[0]+z[1]<=N and z not in seen:seen[z]=w;queue.append(z)
gaps=[n for n in range(1,N+1,2) if n%3 and n not in values]
need(gaps==[7,55,67],'complete coprime-three gap list')
controls=[]
for n in (13,41,53,65):
 word=values[n];rows=build(word);need(hamilton(rows)==n,'literal positive witness')
 controls.append({'target':n,'word':word,'rows':rows})
state_data=json.dumps(sorted(seen),separators=(',',':'))
print(json.dumps({'status':'PASS complete target-bounded inverse in specified grammar','target_bound':N,'expanded_states':len(seen),'state_sha256':hashlib.sha256(state_data.encode()).hexdigest(),'coprime_three_odd_gaps':gaps,'positive_literal_controls':controls},indent=2,sort_keys=True))
