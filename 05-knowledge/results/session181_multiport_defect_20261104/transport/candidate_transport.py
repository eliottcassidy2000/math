import json,hashlib,sys
from pathlib import Path
from common import response,strong
p=Path(__file__).resolve().parent.parent/'multiport'/'candidate.json';c=json.loads(p.read_text());names=c['names'];rows=c['rows'];r=names.index(c['root'])
indices=[names.index(v)for v in c['suffix_names']]
u=[sum(1<<j for j,b in enumerate(indices)if rows[a]>>b&1)for a in indices]
w=[c['suffix_names'].index(v)for v in c['in_ports']]
h,fs=response(u);F=sum(fs[i]for i in w)
# Independent endpoint recurrences, retaining the two directed cross-port responses.
n=len(u);N=1<<n;end=[[0]*n for _ in range(N)];start=[[0]*n for _ in range(N)]
for s in range(1,N):
 for v in range(n):
  if s>>v&1:
   t=s^(1<<v)
   end[s][v]=1 if not t else sum(end[t][a]for a in range(n)if u[a]>>v&1)
   start[s][v]=1 if not t else sum(start[t][a]for a in range(n)if u[v]>>a&1)
joints={(a,b):sum(end[s][a]*start[(N-1)^s][b]for s in range(N))for a in w for b in w}
J=sum(joints.values());H,Ef=response(rows)
if (h,F,J,H,Ef[r])!=(18583,124383,3690,512425,516115):raise RuntimeError('candidate mismatch')
if not strong(rows):raise RuntimeError('not strong')
receipt={'status':'PASS','candidate_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'h':h,'F':F,'J':J,'H':H,'E_root2':Ef[r],'joint_suffix':[{'end':c['suffix_names'][a],'start':c['suffix_names'][b],'value':v}for (a,b),v in joints.items()],'clock_alpha':h-J,'clock_F':F,'constant_root_defect':J,'shared_clock_offset_constant':h-J-3*2195,'shared_clock_offset_power2_coefficient':F-3*42153,'plus4_unique_positive_clock':2,'plus2_positive_clock':None}
Path(__file__).with_name('candidate_transport.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
