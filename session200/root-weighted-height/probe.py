import json
T=(1,0,2,1)
def step(s,c):
 D,E,U,d=s
 if c=='O':return D,0,2*U,d
 if c=='Z':return D-E+U,d,2*U-d,d
 return D-E+U,D-E,D+U,D
X=50001
states={(T,0,0):''};queue=list(states);support={};both=[]
for s,r,k in queue:
 N=s[0]+2*s[2]
 width=N.bit_length()-1
 if not (s[2]-s[3]>=2**k and k<= (N-3).bit_length()-2 and len(states[s,r,k])<= (N-3).bit_length()-2+r):raise RuntimeError('charged height')
 if k in (width-2,width-1):support.setdefault(N,states[s,r,k])
 for c in 'OZJ':
  if c=='J' and r==3:continue
  t=step(s,c);tkey=t,r+(c=='J'),k+(c!='J')
  if t[0]+2*t[2]<=X and tkey not in states:states[tkey]=states[s,r,k]+c;queue.append(tkey)
hostile=[]
for (s,r,k),word in states.items():
 if s[0]+2*s[2]==49937:hostile.append({'word':word,'profile':s,'J_count':r,'retained_count':k,'total_length':len(word)})
print(json.dumps({'hostile49937':hostile,'cap':X,'state_cost_pairs':len(states),'support_size':len(support),'missing':[n for n in range(3,X+1,2) if n not in support]},indent=2))
