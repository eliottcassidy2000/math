from pathlib import Path
from itertools import permutations
import json,hashlib
P=Path(__file__).resolve().parent
r=json.loads((P/'source/retained_sewing_normal.json').read_text())
sourcebank={tuple(x['signature'])for x in r['inventory']}
bank=[]
for line in(P/'independent-bank.csv').read_text().splitlines():
 z=tuple(map(int,line.split(',')));bank.append(z)
if{z[:6]for z in bank}!=sourcebank:raise RuntimeError('full signature set')
statehashes={}
for word in('four','five'):
 states=[list(map(int,line.split(',')))for line in(P/f'states-{word}.csv').read_text().splitlines()]
 h=hashlib.sha256(json.dumps(states,separators=(',',':')).encode()).hexdigest()
 if h!=r[f'order_{word}_census']['state_sha256']:raise RuntimeError(('exact reachable state set',word))
 statehashes[word]=h
# Reconstruct every retained-root word from the literal bank without importing
# either author's helper. The old root remains label0 after every operation.
arcs=0
for w in r['witnesses']:
 A=[0]
 for k in w['word']:
  record=r['inventory'][k];B=record['rows'];b=record['root'];n=len(A);m=len(B)
  C=[0]*(n+m)
  for i in range(n+m):
   for j in range(i+1,n+m):
    if j<n:forward=bool(A[i]>>j&1)
    elif i>=n:forward=bool(B[i-n]>>(j-n)&1)
    else:forward=not(i==0 and j==n+b)
    C[i if forward else j]|=1<<(j if forward else i)
  A=C
 if A!=w['rows']or w['root']!=0:raise RuntimeError('literal retained-root word')
 arcs+=len(A)*(len(A)-1)//2
M=r['continuation_matrix'];det=0
for p in permutations(range(4)):
 inversions=sum(p[i]>p[j]for i in range(4)for j in range(i+1,4));v=(-1)**inversions
 for i in range(4):v*=M[i][p[i]]
 det+=v
if det!=16:raise RuntimeError('determinant')
print(json.dumps(dict(status='PASS',signature_set_size=len(sourcebank),complete_state_hashes=statehashes,word_arcs=arcs,witnesses=len(r['witnesses']),determinant=det),sort_keys=True,indent=2))
