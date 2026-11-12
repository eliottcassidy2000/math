from pathlib import Path
import json,itertools
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
def counts(rows,O):
 n=len(rows);N=1<<n;full=N-1;dp=[[0]*n for _ in range(N)];sp=[[0]*n for _ in range(N)]
 for v in range(n):dp[1<<v][v]=sp[1<<v][v]=1
 for s in range(1,N):
  for v in range(n):
   if not s>>v&1:continue
   t=s^(1<<v)
   if not t:continue
   dp[s][v]=sum(dp[t][u] for u in range(n) if t>>u&1 and rows[u]>>v&1)
   sp[s][v]=sum(sp[t][u] for u in range(n) if t>>u&1 and rows[v]>>u&1)
 H=sum(dp[full]);tax=sum(sum(dp[s][v] for v in O)*sum(sp[full^s][v] for v in O) for s in range(1,full))
 U=sum(sp[full][v] for v in O)+sum(sum(dp[s])*sum(sp[full^s][v] for v in O) for s in range(1,full))
 return H-tax,U,H

def expand(rows,O,ports,lens):
 labels=[]
 for v in range(len(rows)):labels += [v]*(lens[ports.index(v)] if v in ports else 1)
 child=[sum(1<<j for j,w in enumerate(labels) if i!=j and ((v==w and i<j) or (v!=w and rows[v]>>w&1))) for i,v in enumerate(labels)]
 return child,[i for i,v in enumerate(labels) if v in O]

def ev(rows,O,p,a,b):return counts(*expand(rows,O,p,[a,b]))

def load_bank():
 import subprocess,gzip,hashlib
 meta=json.loads((HERE/'SEARCH_PROVENANCE.json').read_text())
 source=ROOT/meta['source']
 raw=source.read_bytes() if source.exists() else subprocess.check_output(['git','show',meta['source_commit']+':'+meta['source']],cwd=ROOT)
 if hashlib.sha256(raw).hexdigest()!=meta['source_sha256']:raise RuntimeError('inherited census source changed; use pinned commit')
 return [json.loads(l) for l in gzip.decompress(raw).splitlines()]
