"""Exact almost-transitive physical endpoint compiler; no external dependencies."""
from pathlib import Path
import json, itertools, hashlib, sys
P=Path(__file__).resolve().parent
from codec import chart

def need(ok,msg):
 if not ok: raise RuntimeError(msg)

def covers(old,beat,slots=1,ends=None,starts=None):
 # old vertices in transitive order; beat[v] means v -> root, root is -1.
 ends=ends or {}; starts=starts or {}; total=0
 for root_slot in range(slots):
  pre,post=0,1
  bins={s:2+s-(s>root_slot) for s in range(slots) if s!=root_slot}
  allowed={v:set(range(slots+1)) for v in old}; empty=set(); impossible=False
  for is_end,restrictions in [(True,ends),(False,starts)]:
   for slot,w in restrictions.items():
    if slot==root_slot:
     b=post if is_end else pre
     if w==-1: empty.add(b);continue
    else:
     if w==-1:impossible=True;break
     b=bins[slot]
    if w not in allowed:impossible=True;break
    allowed[w]&={b}
    for v in old:
     if (v>w if is_end else v<w):allowed[v].discard(b)
  if impossible:continue
  # a=state of latest pre vertex (0 absent,1 beats root,2 loses);
  # b=whether post is nonempty; its first vertex must lose to root.
  dp={(0,0):1}
  for v in old:
   nd={}
   for (a,b),count in dp.items():
    for dest in allowed[v]-empty:
     if dest==pre:key=(1 if beat[v] else 2,b)
     elif dest==post:
      if not b and beat[v]:continue
      key=(a,1)
     else:key=(a,b)
     nd[key]=nd.get(key,0)+count
   dp=nd
  total+=sum(c for (a,b),c in dp.items() if a!=2)
 return total

def rows_for(old,beat):
 n=len(old)+1; rows=[0]*n
 for i in range(len(old)):
  for j in range(i+1,len(old)):rows[i]|=1<<j
  if beat[old[i]]:rows[i]|=1<<(n-1)
  else:rows[-1]|=1<<i
 return rows

def subset(rows):
 n=len(rows);N=1<<n;E=[[0]*n for _ in range(N)];S=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for m in range(1,N):
  for v in range(n):
   if not m>>v&1:continue
   q=m^(1<<v)
   E[m][v]=sum(E[q][u] for u in range(n) if q>>u&1 and rows[u]>>v&1) if q else 1
   S[m][v]=sum(S[q][u] for u in range(n) if q>>u&1 and rows[v]>>u&1) if q else 1
  H[m]=sum(E[m])
 return E,S,H

def strong(rows):
 n=len(rows)
 for v in range(n):
  seen=1<<v
  while True:
   new=seen
   for u in range(n):
    if seen>>u&1:new|=rows[u]
   if new==seen:break
   seen=new
  if seen!=(1<<n)-1:return False
 return True

def joint(f):
 g=j=e=0
 while f:
  b=f&1;f>>=1
  if b:j+=g+e
  g=2*g+b;e=b
 return j

def test():
 cases=fields=joints=0
 for m in range(0,7):
  old=list(range(m));names=old+[-1]
  for cut in range(1<<m):
   beat={v:bool(cut>>v&1) for v in old};rows=rows_for(old,beat);E,S,H=subset(rows);allm=(1<<len(rows))-1
   need(covers(old,beat)==H[allm],'H');cases+=1
   for i,w in enumerate(names):
    for starts in [False,True]:
     table=S if starts else E;kw={'starts' if starts else 'ends':{0:w}}
     need(covers(old,beat,1,**kw)==table[allm][i],'single endpoint')
     expected=sum(table[z][i]*H[allm^z] for z in range(1,allm+1) if z>>i&1)
     need(covers(old,beat,2,**kw)==expected,'double endpoint');fields+=2
    for j,v in enumerate(names):
     if i==j:continue
     expected=sum(E[z][i]*S[allm^z][j] for z in range(1,allm) if z>>i&1 and not z>>j&1)
     need(covers(old,beat,2,ends={0:w},starts={1:v})==expected,'joint');joints+=1
 return dict(tournaments=cases,endpoint_comparisons=fields,joint_comparisons=joints)

def main():
 tests=test();states=[]
 for t,f in [(1,2710066),(1,3035358),(1,3186353),(1,3568813),(2,1355033),(2,1517679)]:
  m=f.bit_length();old=list(range(t+m));beat={v:bool(v>=t and f>>(v-t)&1) for v in old};names=old+[-1];rows=rows_for(old,beat)
  h=covers(old,beat);need(h==4611841 and strong(rows),'target strong')
  need(h==1+(1<<t)*f-joint(f),'formula');records=[];joints=[]
  for w in names:
   e1=covers(old,beat,ends={0:w});s1=covers(old,beat,starts={0:w})
   e2=covers(old,beat,2,ends={0:w});s2=covers(old,beat,2,starts={0:w})
   deletion=1 if w==-1 else covers([v for v in old if v!=w],beat)
   need((s2-e1)%2==0,'golden integrality');g=[h,(s2-e1)//2,(s2+e1)//2]
   records.append(dict(name=w,E1=e1,S1=s1,deletion_H=deletion,E2=e2,S2=s2,defect=e2-h,golden=g,charts=chart(g)))
  for w in names:
   js=[]
   for v in names:
    if w==v:continue
    j=covers(old,beat,2,ends={0:w},starts={1:v});js.append((v,j));joints.append([w,v,j])
   i=names.index(w);need(sum(j for v,j in js if rows[names.index(v)]>>i&1)==records[i]['defect'],'physical defect')
  need(records[-1]['defect']==joint(f),'root defect')
  cut=f<<t; reversed_complement=int(format(((1<<len(old))-1)^cut,'0'+str(len(old))+'b')[::-1],2)
  need(records[-1]['E2']==1+cut and records[-1]['S2']==1+reversed_complement,'root binary endpoints')
  states.append(dict(t=t,f=f,j=joint(f),H=h,names=names,rows=rows,in_ports=[v for v in old if beat[v]],ports=records,joints=joints))
 out=dict(status='PASS',tests=tests,states=states)
 (P/'receipt.json').write_text(json.dumps(out,indent=2)+'\n')
 (P/'candidate.json').write_text(json.dumps(states[0],indent=2)+'\n')
 print(json.dumps({'status':'PASS','tests':tests,'representatives':len(states),'orders':[len(s['names']) for s in states],'first_root':states[0]['ports'][-1]},sort_keys=True))
if __name__=='__main__':main()
