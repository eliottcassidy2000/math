from collections import Counter,deque
import argparse,json,math,sys
from pathlib import Path
from forest_runtime import introduce,forget,literal,BULK

def m(t):return math.isqrt(4*t-1)if t else 0
def width(K):return 1+max(m(t)+m(K-t)for t in range(K+1))
def vector_key(v):return tuple(sorted((p,n)for p,n in v.items()if n))
def saturate(v,c):return Counter({p:min(c,n)for p,n in v.items()if n})
def bounded_read(v,r,c):
 B=s=0
 for parts,n in v.items():
  if len(parts)==1:
   if parts[0][0]==r:s+=n
   else:B+=n
  if len(parts)in(1,2)and any(e==r for a,e in parts):B+=n
 return min(c,B),min(c,s)

def need(ok,message):
 if not ok:raise RuntimeError(message)
def replay(parents,stateid,root):
 word=[]
 while parents[stateid] is not None:
  stateid,op=parents[stateid];word.append(op)
 word.reverse();active={};rows=[]
 for op in word:
  if op[0]=='forget':
   need(op[1] in active,'forget inactive witness port');del active[op[1]]
  else:
   _,u,A=op;need(u not in active and all(a in active for a in A),'invalid witness introduction')
   selected={active[a]for a in A};row=0;n=len(rows)
   for v in range(n):
    if v in selected:row|=1<<v
    else:rows[v]|=1<<n
   rows.append(row);active[u]=n
 need(root in active,'missing witness root')
 return rows,active[root],word

def search(K,X,limit):
 W=width(K);c=max(K,X)+1;v0=Counter({():1});key=(0,vector_key(v0));queue=deque([key]);seen={key:0};parents=[None];states=[key];outputs={};transitions=0
 def finish(status):
  witnesses=[]
  for (B,s),(stateid,r) in sorted(outputs.items()):
   rows,root,word=replay(parents,stateid,r)
   exact=literal(rows,root);need((exact[3],exact[1])==(B,s),'literal witness profile mismatch')
   witnesses.append({'B':B,'s':s,'root':root,'rows':rows,'word':word})
  return {'status':status,'K':K,'X':X,'W':W,'cutoff':c,'states':len(seen),'transitions':transitions,'pairs':sorted(outputs),'literal_witnesses_checked':len(witnesses),'witnesses':witnesses}
 while queue:
  active,keyv=queue.popleft();v=Counter(dict(keyv));stateid=seen[(active,keyv)]
  for r in range(W):
   if active>>r&1:
    B,s=bounded_read(v,r,c)
    if 1<=B<=X and 1<=s<=K:outputs.setdefault((B,s),(stateid,r))
  successors=[]
  A=[u for u in range(W)if active>>u&1]
  for u in range(W):
   if active>>u&1:successors.append((active^(1<<u),forget(v,u,2),('forget',u)))
   else:
    for mask in range(1<<len(A)):
     chosen={A[j]for j in range(len(A))if mask>>j&1}
     successors.append((active|(1<<u),introduce(v,u,chosen),('introduce',u,tuple(sorted(chosen)))))
  for aa,vv,op in successors:
   transitions+=1;kk=(aa,vector_key(saturate(vv,c)))
   if kk in seen:continue
   if len(seen)>=limit:return finish('INCOMPLETE')
   seen[kk]=len(states);states.append(kk);parents.append((stateid,op));queue.append(kk)
 return finish('COMPLETE')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--K',type=int,default=1);p.add_argument('--X',type=int,default=5);p.add_argument('--max-states',type=int,default=100000);a=p.parse_args()
 if a.K<1 or a.X<1 or a.max_states<1:raise SystemExit('positive integer arguments required')
 print(json.dumps(search(a.K,a.X,a.max_states),indent=2))
