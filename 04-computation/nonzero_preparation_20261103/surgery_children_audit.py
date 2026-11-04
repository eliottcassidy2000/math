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

receipt=json.loads(Path('05-knowledge/results/nonzero_preparation_20261103/surgery/children.json').read_text())
num_native=0; weak_completed=[]
for ci,case in enumerate(receipt['cases']):
 names=case['names'];r=names.index('r');need(r==8,'physical root')
 current=tuple(case['source_rows'])
 for stage,record in enumerate(case['program']):
  local=names+[record['helper']];cut=sum(1<<i for i,name in enumerate(names) if name not in record['helper_out_names'])
  states=[add(current,cut)]
  need(list(states[0])==record['birth_rows'],'actual helper birth')
  for guard in record['guards']:
   states.append(native(states[-1],[local.index(v) for v in guard]));num_native+=1
  need(list(states[1])==record['after_first_rows'] and list(states[2])==record['after_second_rows'],'literal native matrices')
  for state,obs in zip(states,record['native_root_observers']):
   need(strong(state),'independent strong native state')
   need(observe(state,r)==(obs['H'],obs['F_root']) and obs['defect']==obs['F_root']-obs['H'],'independent native observers')
  mask=(1<<len(names))-1;after=tuple(v&mask for v in states[-1][:-1])
  a,b=[names.index(v) for v in record['old_pair']]
  need(after==flip(current,a,b),'only traced pair changed')
  need(list(after)==record['after_helper_deletion_rows'],'only helper deleted')
  if not strong(after):weak_completed.append((ci,stage))
  obs=case['completed_adjacency_root_observers'][stage+1]
  need(observe(after,r)==(obs['H'],obs['F_root']),'completed observers')
  current=after
 need(list(current)==case['target_rows'],'complete target')
 need(observe(current,r)==(case['joint']['H'],case['joint']['F']),'target observers')
 cfg=case['construction'];core=tuple(cfg['core']);ports=set(cfg['ports']);word=cfg['word'];mapping=cfg['mapping']
 n=len(core);t=len(word);rr=n+t
 funnel=[row|(((1<<t)-1)<<n)|(0 if i in ports else 1<<rr) for i,row in enumerate(core)]
 funnel += [sum(1<<(n+j) for j in range(i+1,t))|((1<<rr) if bit else 0) for i,bit in enumerate(word)]
 funnel += [sum(1<<v for v in ports)+sum(1<<(n+i) for i,bit in enumerate(word) if not bit)]
 canonical=add(tuple(funnel),1<<rr);full=(1<<len(canonical))-1
 canonical=tuple(full^(1<<i)^row for i,row in enumerate(canonical))
 if cfg['padding']:
  canonical=((1<<(len(canonical)+1))-2,)+tuple(row<<1 for row in canonical);rr+=1
 mapped=[0]*len(canonical)
 for i,row in enumerate(canonical):mapped[mapping[i]]=sum(1<<mapping[j] for j in range(len(canonical)) if row>>j&1)
 need(tuple(mapped)==current,'canonical-to-physical reconstruction')
 need(mapping[rr]==r and mapping[-1]==case['aux_index'],'retained root and auxiliary')
 need(current[case['aux_index']]==1<<r,'actual sole exit')
 first=add(current,((1<<len(current))-1)^(1<<case['aux_index']))
 child=add(first,(1<<len(current))-1)
 need(list(first)==case['C2_first_rows'] and list(child)==case['child_rows'],'actual permanent C2 births')
 need(observe(child,r)==(case['child_H'],case['child_F']),'child observers')
 need(case['child_H']==3*case['source_H']+2*case['gap'] and case['child_F']==3*case['joint']['F'],'joint child response')
 for state in (current,first,child):
  need(strong(state)==(not cfg['padding']),'exact strongness scope')
  if cfg['padding']:
   z=mapping[0];full=(1<<len(state))-1
   need(state[z]==full^(1<<z),'retained source prefix')
   ids=[j for j in range(len(state)) if j!=z]
   reduced=tuple(sum(1<<i for i,j in enumerate(ids) if state[v]>>j&1) for v in ids)
   need(strong(reduced),'only source singleton outside strong core')
need(num_native==46,'complete 46 complements')
need(weak_completed==[(3,5)],'unique weak completed stage')
print(json.dumps({'gates':G,'native_complements':num_native,'native_states':69,'physical_root':8,'canonical_maps':4,'weak_completed':weak_completed},sort_keys=True))
