import json,sys,hashlib
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'ports'))
from verify import subset,strong,need
P=Path(__file__).resolve().parent

def tt(n):return [(1<<n)-(1<<(i+1)) for i in range(n)]
def join(a,b):return [r|(((1<<len(b))-1)<<len(a)) for r in a]+[r<<len(a) for r in b]
def attach(a,b,W):
 rows=join(a,b);n=len(rows);root=(1<<len(a))-1
 for i in range(len(b)):
  if i in W:rows[len(a)+i]|=1<<n
  else:root|=1<<(len(a)+i)
 return rows+[root]
def born(names,rows,x,incoming):
 n=len(names);new=0
 for i,v in enumerate(names):
  if v in incoming:rows[i]|=1<<n
  else:new|=1<<i
 rows.append(new);names.append(x)
def canonical(names,rows,wanted):
 return [sum(1<<j for j,v in enumerate(wanted) if rows[names.index(u)]>>names.index(v)&1) for u in wanted]
checks=0;short=0
for n in range(4):
 pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
 for code in range(1<<len(pairs)):
  L=[0]*n
  for k,(i,j) in enumerate(pairs):
   a,b=(i,j) if code>>k&1 else (j,i);L[a]|=1<<b
  h=subset(L)[2][-1]
  for mask in range(1<<n):
   U={i for i in range(n) if mask>>i&1}
   for A in [tt(1),tt(2),[2,4,1]]:
    an=[('a',i) for i in range(len(A))];ln=[('l',i) for i in range(n)];bn=[('b',i) for i in range(5)];r=('r',0)
    names=an+ln+bn[:3]+[r];rows=attach(A,join(L,tt(3)),U|{n,n+1,n+2})
    born(names,rows,bn[4],set(an+ln+bn[:2]+[r]));need(strong(rows),'first birth strong')
    born(names,rows,bn[3],set(an+ln+bn[:3]+[r]));need(strong(rows),'second birth strong')
    expected=attach(A,join(L,join(tt(2),[2,4,1])),U|{n,n+1,n+2})
    need(canonical(names,rows,an+ln+bn+[r])==expected,'literal final adjacency');checks+=1
    # Independent endpoint DP tests the unexpanded shorter family too.
    src=attach(A,join(L,tt(2)),U|{n,n+1})
    dst=attach(A,join(L,join(tt(1),[2,4,1])),U|{n,n+1})
    hs=subset(src)[2][-1];hd=subset(dst)[2][-1];ha=subset(A)[2][-1]
    need(hd-3*hs==2*ha*h and strong(src) and strong(dst),'short +2 law');short+=1
out={'status':'PASS','literal_two_birth_schedules':checks,'independent_short_family_full_graph_pairs':short,'source_sha256':{f:hashlib.sha256((P.parent/'mining'/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','FRONTIER.md']}}
(P/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
