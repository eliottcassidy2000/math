from pathlib import Path
import json, hashlib, itertools, re
root=Path('/workspace/scratch/07d5532f2fce/session26-root')
shared=Path('/workspace/shared')
a=json.loads((root/'ROOT_SAMPLE.json').read_text())
b=json.loads((shared/'session26-tournament/random-manifest.json').read_text())
c=json.loads((shared/'session26-runners/sample-manifest.json').read_text())
d=json.loads((shared/'session26-wildcard/manifest.json').read_text())
paths=a['files']+b+c['paths']+[x['path'] for x in d['samples']]
if len(paths)!=103 or len(set(paths))!=103: raise ValueError('sample counts')
for i,entry in enumerate(a['snapshot_verification']):
 if entry['path']!=a['files'][i]: raise ValueError('path binding')
 if hashlib.sha256((root/f'sample_{i:02}.md').read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('hash binding')
links=[]
for file in ['SYNTHESIS.md','EXPLORATION.md']:
 for target in re.findall(r'\]\(([^)]+)\)',(root/file).read_text()):
  choices=[root/target, shared/target,root.parent/target]
  if not any(p.exists() for p in choices):raise ValueError(target)
  links.append(target)
# Exact scalar identity, averaging, symmetries, and direct Hamiltonian DP.
checks=0
for m in range(1,10):
 weights=lambda d: 2 if d==1 else 2**(d-1)
 total_weight=sum(weights(i-j) for j in range(m) for i in range(j+1,m))
 if total_weight!=2**m-2:raise ValueError('weight sum')
 total=0
 for b in itertools.product((0,1),repeat=m):
  d=[sum(weights(i-j) for i in range(j+1,m))-sum(weights(j-i) for i in range(j)) for j in range(m)]
  energy=sum(weights(i-j)*(b[i]-b[j])**2 for j in range(m) for i in range(j+1,m))
  field=sum(x*y for x,y in zip(d,b))
  h=1+sum(b[j]*(1-b[i])*weights(i-j) for j in range(m) for i in range(j+1,m))
  if 2*h!=2+energy+field: raise ValueError('identity')
  comp=tuple(1-x for x in b)
  hc=1+sum(comp[j]*(1-comp[i])*weights(i-j) for j in range(m) for i in range(j+1,m))
  if h+hc!=2+energy or h-hc!=field:raise ValueError('complement')
  # Direct subset DP, independent of cut recurrence.
  n=m+1; dp={ (1<<v,v):1 for v in range(n)}
  for mask in range(1,1<<n):
   for v in range(n):
    z=dp.get((mask,v),0)
    if not z:continue
    for u in range(n):
     if mask>>u&1:continue
     edge=(v<u if v<m and u<m else (bool(b[u]) if v==m else not b[v]))
     if edge:dp[mask|1<<u,u]=dp.get((mask|1<<u,u),0)+z
  hp=sum(dp.get(((1<<n)-1,v),0) for v in range(n))
  if hp!=h:raise ValueError('DP')
  total+=h;checks+=1
 if 4*total!=(2**m)*(2**m+2):raise ValueError('mean')
rows=[]
for n in (8,9,10):
 A=2**(n-1)-(n-1)*(n-2)//2-1
 V=4**(n-2)-3*2**(n-2)-(n*n-7*n+4)//2
 r=n//2; C=r if n%2==0 else 2**r+r-1
 B=2**r+r-2 if n%2==0 else 2**(2*r-1)+r-1
 rows.append([n,(A+C)//2,V,B,(V-B)//2-(12 if n==8 else 0)])
print(json.dumps({'status':'PASS','source_sha256':{f:hashlib.sha256((root/f).read_bytes()).hexdigest() for f in ['SYNTHESIS.md','EXPLORATION.md']},'sample_paths':len(paths),'unique_paths':len(set(paths)),'root_copy_hashes':len(a['snapshot_verification']),'resolved_links':len(links),'electrical_words_direct_DP':checks,'boundary_rows':rows},indent=2))
