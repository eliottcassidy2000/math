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

receipt=json.loads(Path('05-knowledge/results/nonzero_preparation_20261103/surgery/direct_children.json').read_text())
prior=json.loads(Path('05-knowledge/results/nonzero_preparation_20261103/surgery/direct.json').read_text())
old=tuple(prior['cases'][0]['target_rows']);a=prior['cases'][0]['aux_index']
first=add(old,((1<<len(old))-1)^(1<<a));source=add(first,(1<<len(old))-1)
for T,rec in zip((old,first,source),receipt['common_prefix']):
 need(tuple(rec['rows'])==T and strong(T),'actual strong shared C2 prefix')
 need(observe(T,8)==(rec['H'],rec['F']) and rec['defect']==rec['F']-rec['H'],'independent shared prefix observers')
need(observe(old,8)==(4289,7653) and observe(source,8)==(12871,22959),'actual common marked source')
count=0;states_count=0

def verify_program(rec,original):
 global count,states_count
 names=rec['names'];r=names.index('r');need(r==8,'fixed physical root')
 current=original;need(tuple(rec['source_rows'])==original,'same actual source graph')
 for stage,record in enumerate(rec['program']):
  local=names+[record['helper']];cut=sum(1<<i for i,name in enumerate(names) if name not in record['helper_out_names'])
  states=[add(current,cut)];need(tuple(record['birth_rows'])==states[0],'literal helper birth')
  for guard in record['guards']:
   states.append(native(states[-1],[local.index(v) for v in guard]));count+=1
  need(tuple(record['after_first_rows'])==states[1] and tuple(record['after_second_rows'])==states[2],'literal native states')
  for T,obs in zip(states,record['native_root_observers']):
   states_count+=1;need(strong(T),'every native state strong')
   need(observe(T,r)==(obs['H'],obs['F_root']) and obs['defect']==obs['F_root']-obs['H'],'every independent native root observer')
  mask=(1<<len(names))-1;after=tuple(row&mask for row in states[-1][:-1])
  a,b=[names.index(v) for v in record['old_pair']]
  need(after==flip(current,a,b),'only actual old pair changes')
  need(tuple(record['after_helper_deletion_rows'])==after and strong(after),'only helper deleted, completed state strong')
  obs=rec['completed_adjacency_root_observers'][stage+1]
  need(observe(after,r)==(obs['H'],obs['F_root']),'every completed observer')
  current=after
 need(tuple(rec['target_rows'])==current,'actual full program endpoint')
 return current

for rec in receipt['cases']:
 T=verify_program(rec,source);cfg=rec['construction']
 core=tuple(cfg['core']);ports=set(cfg['ports']);word=cfg['word'];mapping=cfg['mapping'];n=len(core);t=len(word);r=n+t
 rows=[row|(((1<<t)-1)<<n)|(0 if i in ports else 1<<r) for i,row in enumerate(core)]
 rows += [sum(1<<(n+j) for j in range(i+1,t))|((1<<r) if bit else 0) for i,bit in enumerate(word)]
 rows += [sum(1<<v for v in ports)+sum(1<<(n+i) for i,bit in enumerate(word) if not bit)]
 canonical=reverse(add(tuple(rows),1<<r));mapped=[0]*len(canonical)
 for i,row in enumerate(canonical):mapped[mapping[i]]=sum(1<<mapping[j] for j in range(len(canonical)) if row>>j&1)
 need(tuple(mapped)==T and mapping[r]==8 and mapping[-1]==rec['aux_index'],'actual canonical-to-physical port transport')
 need(T[rec['aux_index']]==1<<8 and joint(T,8,rec['aux_index'])==rec['joint'],'independent full positive joint guard')
 need(len(T)==15 and rec['joint']['H']==3*4289+rec['digit_from4289'],'retained order and both actual count branches')
B=tuple(prior['cases'][1]['target_rows']);C=tuple(prior['cases'][2]['target_rows']);switch=receipt['same_H_rearm']
need(verify_program(switch,B)==C,'switch uses the single actual B child and reaches actual C')
need(joint(B,8,switch['source_aux'])==switch['source_joint'] and joint(C,8,switch['target_aux'])==switch['target_joint'],'independent changed joint response through same-H switch')
need(switch['source_joint']['H']==switch['target_joint']['H']==4291,'endpoint-H-neutrality only')
expected=[tuple(rec['target_rows']) for rec in receipt['cases']]
for T,a in ((B,switch['source_aux']),(C,switch['target_aux'])):
 mid=add(T,((1<<len(T))-1)^(1<<a));leaf=add(mid,(1<<len(T))-1)
 need(strong(mid) and strong(leaf),'both final C2 stages strong')
 expected.append(leaf)
for T,rec in zip(expected,receipt['four_leaves']):
 need(T==tuple(rec['rows']) and len(T)==15 and strong(T),'exact common-state physical leaf')
 need(observe(T,8)==(rec['H'],rec['F']) and rec['defect']==rec['F']-rec['H'],'all actual final leaf observers')
need({x['H'] for x in receipt['four_leaves']}=={12869,12871,12875,12877},'all four second-generation scalar values')
need([x['positive_guard'] for x in receipt['four_leaves']]==[1,2,None,None],'prepared-vs-C2-leaf scope retained')
need(count==52 and states_count==78,'complete final program universe')
print(json.dumps({'gates':G,'programs':3,'native_complements':count,'native_states':states_count,'canonical_maps':2,'shared_state_switch':True,'physical_root':8,'strong_leaves':4},sort_keys=True))
