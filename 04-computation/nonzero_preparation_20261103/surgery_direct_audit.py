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

def reverse(T):
 full=(1<<len(T))-1
 return tuple(full^(1<<i)^row for i,row in enumerate(T))
def induced(T,ids):return tuple(sum(1<<j for j,w in enumerate(ids) if T[v]>>w&1) for v in ids)
def tables(T):
 n=len(T);full=(1<<n)-1;E=[[0]*n for _ in range(full+1)];H=[0]*(full+1);H[0]=1
 for i in range(n):E[1<<i][i]=1
 for m in range(1,full+1):
  for v in range(n):
   value=E[m][v]
   if not value:continue
   H[m]+=value;q=T[v]&(full^m)
   while q:
    b=q&-q;q-=b;E[m|b][b.bit_length()-1]+=value
 return H,E

def seven(T,ports):
 H,E=tables(T);_,A=tables(reverse(T));full=(1<<len(T))-1
 e=[sum(row[v] for v in ports) for row in E];a=[sum(row[v] for v in ports) for row in A]
 return (H[full],sum(H[m]*H[full^m] for m in range(full+1)),sum(e[m]*H[full^m] for m in range(full+1)),sum(a[m]*H[full^m] for m in range(full+1)),sum(e[m]*a[full^m] for m in range(full+1)),e[full],a[full])

def joint(T,r,a):
 h,f=observe(T,r);s=observe(reverse(T),a)[1]
 ids=[v for v in range(len(T)) if v!=a];L=observe(induced(T,ids),ids.index(r))[1]
 ids=[v for v in ids if v!=r];H,_=tables(induced(T,ids));full=len(H)-1;K=sum(H[m]*H[full^m] for m in range(full+1))
 return dict(H=h,F=f,S=s,L=L,K=K,G=L+K,root_defect=f-h,sigma=s-h,tau=L+K-f)

receipt=json.loads(Path('05-knowledge/results/nonzero_preparation_20261103/surgery/direct.json').read_text())
old=(1790,1668,1672,1778,1670,1942,1718,1280,1631,128,512)
first=add(old,((1<<11)-1)^(1<<9));source=add(first,(1<<11)-1)
need(tuple(receipt['common_prefix'][0]['rows'])==old,'inherited literal1429 source')
for state,record in zip((old,first,source),receipt['common_prefix']):
 need(tuple(record['rows'])==state,'literal permanent zero-C2 births')
 need(strong(state),'strong zero prefix')
 need(observe(state,8)==(record['H'],record['F']) and record['defect']==record['F']-record['H'],'independent zero-prefix observers')
need(observe(source,8)==(4287,4689),'zero-C2 full response')
count=0
for ci,case in enumerate(receipt['cases']):
 names=case['names'];r=names.index('r');need(r==8,'physical root')
 current=source;need(tuple(case['source_rows'])==source,'same actual common source')
 for stage,record in enumerate(case['program']):
  local=names+[record['helper']];cut=sum(1<<i for i,name in enumerate(names) if name not in record['helper_out_names'])
  states=[add(current,cut)];need(tuple(record['birth_rows'])==states[0],'actual helper birth')
  for guard in record['guards']:
   states.append(native(states[-1],[local.index(v) for v in guard]));count+=1
  need(tuple(record['after_first_rows'])==states[1] and tuple(record['after_second_rows'])==states[2],'literal native states')
  for state,obs in zip(states,record['native_root_observers']):
   need(strong(state),'independent native strongness')
   need(observe(state,r)==(obs['H'],obs['F_root']) and obs['defect']==obs['F_root']-obs['H'],'independent native observers')
  mask=(1<<len(names))-1;after=tuple(row&mask for row in states[-1][:-1])
  a,b=[names.index(v) for v in record['old_pair']]
  need(after==flip(current,a,b),'exact old-pair change and helper-only deletion')
  need(tuple(record['after_helper_deletion_rows'])==after and strong(after),'strong completed state')
  obs=case['completed_adjacency_root_observers'][stage+1]
  need(observe(after,r)==(obs['H'],obs['F_root']),'completed observers')
  current=after
 need(tuple(case['target_rows'])==current,'literal ready target')
 cfg=case['construction'];core=tuple(cfg['core']);ports=set(cfg['ports']);word=cfg['word'];mapping=cfg['mapping'];n=len(core);t=len(word);rr=n+t
 funnel=[row|(((1<<t)-1)<<n)|(0 if i in ports else 1<<rr) for i,row in enumerate(core)]
 funnel += [sum(1<<(n+j) for j in range(i+1,t))|((1<<rr) if bit else 0) for i,bit in enumerate(word)]
 funnel += [sum(1<<v for v in ports)+sum(1<<(n+i) for i,bit in enumerate(word) if not bit)]
 canonical=reverse(add(tuple(funnel),1<<rr));mapped=[0]*len(canonical)
 for i,row in enumerate(canonical):mapped[mapping[i]]=sum(1<<mapping[j] for j in range(len(canonical)) if row>>j&1)
 need(tuple(mapped)==current and mapping[rr]==8 and mapping[-1]==case['aux_index'],'complete marked canonical map')
 a=case['aux_index'];need(current[a]==1<<8,'literal sole exit to retained root')
 need(joint(current,8,a)==case['joint'],'independent all joint observers')
 first=add(current,((1<<13)-1)^(1<<a));child=add(first,(1<<13)-1)
 need(tuple(case['next_first_rows'])==first and tuple(case['next_child_rows'])==child,'two subsequent permanent births')
 need(strong(first) and strong(child),'strong immediate next C2 steps')
 need(observe(child,8)==(case['next_H'],case['next_F']) and case['next_defect']==case['next_F']-case['next_H'],'independent next child observers')
 unif=receipt['uniformity'][ci];ids=[v for v in range(1,13) if v!=8]
 suffix=induced(current,ids);ports={i for i,v in enumerate(ids) if current[v]>>8&1}
 sig=seven(suffix,ports)
 need(tuple(unif['signature'])==sig and tuple(unif['suffix_rows'])==suffix and unif['suffix_ports']==sum(1<<v for v in ports),'independent full suffix signature')
 need(sig[2]-2268==unif['error_coefficient'] and 2*sig[2]+sig[0]-sig[4]==cfg['H'],'exact dyadic affine error')
 k=2;m=len(suffix);root=k+m
 inflated=[sum(1<<j for j in range(i+1,k))|(((1<<m)-1)<<k) for i in range(k)]
 inflated += [(row<<k)|((1<<root) if i in ports else 0) for i,row in enumerate(suffix)]
 inflated += [((1<<k)-1)+sum(1<<(k+i) for i in range(m) if i not in ports)]
 inflated=tuple(inflated);ia=k+ids.index(a)
 need(strong(inflated) and joint(inflated,root,ia)==unif['clock2_observers'],'literal inflated-clock hostile')
need(count==44,'complete44 guard universe')
print(json.dumps({'gates':G,'native_complements':count,'native_states':66,'canonical_maps':3,'suffix_signatures':3,'physical_root':8},sort_keys=True))
