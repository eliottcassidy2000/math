import json,pathlib
D=pathlib.Path(__file__).resolve().parent
c=json.load(open(D.parent/'missing44/candidate44.json')); src=json.load(open(D.parent/'missing44/source.json'))
names=c['names']; sn=c['suffix_names']; n=len(sn); ALL=(1<<n)-1
T=[sum(1<<j for j,b in enumerate(sn) if c['rows'][names.index(a)]>>names.index(b)&1)for a in sn]
E=[{}for _ in range(1<<n)]; B=[{}for _ in range(1<<n)]; H=[0]*(1<<n);H[0]=1
for m in range(1,1<<n):
 for v in range(n):
  if m>>v&1:
   sub=m^(1<<v)
   E[m][v]=sum(x for u,x in E[sub].items() if T[u]>>v&1) if sub else 1
   B[m][v]=sum(x for u,x in B[sub].items() if T[v]>>u&1) if sub else 1
 H[m]=sum(E[m].values())
f=[sum(E[m].get(v,0)*H[ALL^m]for m in range(1,ALL+1))for v in range(n)]
j=[[sum(E[m].get(u,0)*B[ALL^m].get(v,0)for m in range(1,ALL+1))for v in range(n)]for u in range(n)]
P=[sn.index(v)for v in c['in_ports']];F=sum(f[v]for v in P);J=sum(j[u][v]for u in P for v in P)
if (H[ALL],F,J,H[ALL]+8*F-J)!=(40085,633897,499420,4611841):raise RuntimeError('Counts mismatch')
if names!=src['names']:raise RuntimeError('Names changed')
changes=[]
for i,a in enumerate(names):
 for k in range(i):
  b=names[k]
  if ((src['rows'][i]^c['rows'][i])>>k)&1:changes.append(sorted([a,b]))
if sorted(changes)!=sorted(sorted(x)for x in c['all_changed_pairs']):raise RuntimeError('Unexpected changed pairs')
prefix=c['prefix'];suffix=c['suffix_names'];root=c['root'];arc=lambda a,b:c['rows'][names.index(a)]>>names.index(b)&1
if not all(arc(a,b)for i,a in enumerate(prefix)for b in prefix[i+1:]):raise RuntimeError('Prefix not transitive')
if not all(arc(root,a)for a in prefix):raise RuntimeError('Root/prefix')
if not all(arc(a,b)for a in prefix for b in suffix):raise RuntimeError('Prefix/suffix')
if sorted(a for a in suffix if arc(a,root))!=sorted(c['in_ports']):raise RuntimeError('Port set mismatch')
result=dict(status='PASS_INDEPENDENT_PORT_DP',h=H[ALL],F=F,J=J,H=H[ALL]+8*F-J,root_E2=H[ALL]+8*F,defect=J,changes=changes,all_old_names_preserved=True,prefix_unchanged=True,responses=dict(zip(sn,f)),joint_names=sn,joint=j)
(D/'ports_audit.json').write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps({k:v for k,v in result.items()if k not in ['joint']}))
