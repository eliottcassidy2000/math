import json
from pathlib import Path
from functools import cache
G=0

def need(c,label):
 global G
 G+=1
 if not c: raise RuntimeError(label)

def strong(rows):
 for i in range(len(rows)):
  reached={i}; stack=[i]
  while stack:
   v=stack.pop()
   for w in range(len(rows)):
    if rows[v]>>w&1 and w not in reached: reached.add(w);stack.append(w)
  if len(reached)!=len(rows):return False
 return True

@cache
def observe(rows,r):
 n=len(rows);size=1<<n;full=size-1
 dp=[[0]*n for _ in range(size)];P=[0]*size;P[0]=1
 for i in range(n):dp[1<<i][i]=1
 for m in range(1,size):
  for v in range(n):
   value=dp[m][v]
   if not value:continue
   P[m]+=value
   todo=rows[v]&(full^m)
   while todo:
    b=todo&-todo;todo-=b;dp[m|b][b.bit_length()-1]+=value
 return P[full],sum(dp[m][r]*P[full^m] for m in range(size))

def add(rows,cut):
 n=len(rows);return tuple(row | (((cut>>i)&1)<<n) for i,row in enumerate(rows))+(((1<<n)-1)^cut,)

def flip(rows,a,b):
 out=list(rows);out[a]^=1<<b;out[b]^=1<<a;return tuple(out)

def native(rows,word):
 n=len(rows);full=(1<<n)-1
 need(sorted(word)==list(range(n)),'exact guard names')
 need(all(rows[a]>>b&1 for a,b in zip(word,word[1:])),'guard directed')
 out=tuple(full^(1<<i)^row for i,row in enumerate(rows))
 for a,b in zip(word,word[1:]):out=flip(out,a,b)
 return out

receipt=json.loads(Path('05-knowledge/results/nonzero_preparation_20261103/surgery/native.json').read_text())
num_native=0
for case in receipt['cases']:
 names=case['names'];r=names.index('r');need(r==8,'physical root')
 current=tuple(case['source_rows'])
 for stage,record in enumerate(case['program']):
  local=names+[record['helper']];cut=sum(1<<i for i,name in enumerate(names) if name not in record['helper_out_names'])
  states=[add(current,cut)]
  need(list(states[0])==record['birth_rows'],'actual helper birth')
  for i,guard in enumerate(record['guards']):
   states.append(native(states[-1],[local.index(v) for v in guard]));num_native+=1
  need(list(states[1])==record['after_first_rows'] and list(states[2])==record['after_second_rows'],'literal native matrices')
  for state,obs in zip(states,record['native_root_observers']):
   need(strong(state),'independent strong native state')
   need(observe(state,r)==(obs['H'],obs['F_root']) and obs['defect']==obs['F_root']-obs['H'],'independent native H F defect')
  mask=(1<<len(names))-1;after=tuple(v&mask for v in states[-1][:-1])
  a,b=[names.index(v) for v in record['old_pair']]
  need(after==flip(current,a,b),'only traced pair changed')
  need(list(after)==record['after_helper_deletion_rows'],'only helper deleted')
  need(strong(after),'strong completed stage')
  obs=case['completed_adjacency_root_observers'][stage+1]
  need(observe(after,r)==(obs['H'],obs['F_root']),'completed-stage observers')
  current=after
 need(list(current)==case['target_rows'],'entire native target')
 need(observe(tuple(case['original_rows']),r)==(1429,1563),'old source observers')
 need(observe(current,r)==(case['joint']['H'],case['joint']['F']),'ready target observers')
 need(observe(tuple(case['child_rows']),r)==(case['child_H'],case['child_F']),'literal child observers')
 # Independent canonical funnel construction and explicit physical mapping.
 if case['gap']==2:core=(2,4,9,3);ports={1,2};word=(0,1,1,1,1);mapping=[2,1,9,10,5,4,6,3,0,8,7]
 else:core=(0,);ports={0};word=(0,0,1,1,1,1,1,0,1);mapping=[10,7,9,2,1,4,5,6,3,0,8,11]
 n=len(core);t=len(word);rr=n+t
 funnel=[row|(((1<<t)-1)<<n)|(0 if i in ports else 1<<rr) for i,row in enumerate(core)]
 funnel += [sum(1<<(n+j) for j in range(i+1,t))|((1<<rr) if bit else 0) for i,bit in enumerate(word)]
 funnel += [sum(1<<v for v in ports)+sum(1<<(n+i) for i,bit in enumerate(word) if not bit)]
 canonical=add(tuple(funnel),1<<rr);full=(1<<len(canonical))-1
 canonical=tuple(full^(1<<i)^row for i,row in enumerate(canonical))
 mapped=[0]*len(canonical)
 for i,row in enumerate(canonical):mapped[mapping[i]]=sum(1<<mapping[j] for j in range(len(canonical)) if row>>j&1)
 need(tuple(mapped)==current,'actual canonical-to-physical funnel isomorphism')
 need(mapping[rr]==r and mapping[-1]==case['aux_index'],'root and auxiliary mapping')
need(num_native==24,'complete guard universe')
print(json.dumps({'gates':G,'native_complements':num_native,'native_states':36,'physical_root':8,'canonical_maps':2},sort_keys=True))
