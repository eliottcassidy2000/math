from pathlib import Path
DATA_DIR = Path(__file__).resolve().parent
import json
from itertools import combinations
for line in open(DATA_DIR / 'target_witnesses.jsonl'):
 r=json.loads(line);n,code,cut=r['seed'];word=r['word'];N=n+len(word)+1;out=[0]*N
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a,b=(i,j)if code>>k&1 else(j,i);out[a]|=1<<b
 for i in range(n+len(word)):
  for j in range(max(n,i+1),n+len(word)):out[i]|=1<<j
 for i,b in enumerate([(cut>>i)&1 for i in range(n)]+word):
  a,c=(N-1,i)if b else(i,N-1);out[a]|=1<<c
 # flattened subset endpoint DP, independent of ear response profile.
 d=[[0]*N for _ in range(1<<N)]
 for i in range(N):d[1<<i][i]=1
 for S in range(1,1<<N):
  for v,c in enumerate(d[S]):
   if c:
    rem=out[v]&~S
    while rem:
     b=rem&-rem;rem-=b;u=b.bit_length()-1;d[S|b][u]+=c
 actual=sum(d[-1])
 if actual!=r['target']:raise RuntimeError((r,actual))
 print('literal_tournament_DP',N,actual,'PASS')
