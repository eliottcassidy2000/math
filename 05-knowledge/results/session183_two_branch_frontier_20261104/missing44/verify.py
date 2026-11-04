import json,itertools,pathlib,hashlib,collections
D=pathlib.Path(__file__).parent
checks=[]
for digit in [44]:
 c=json.loads((D/f'candidate{digit}.json').read_text());ns=c['names'];sn=c['suffix_names'];r=[sum(1<<j for j,b in enumerate(sn)if c['rows'][ns.index(a)]>>ns.index(b)&1)for a in sn];n=len(r);N=1<<n
 end=[[0]*n for _ in range(N)];start=[[0]*n for _ in range(N)];h=[1]+[0]*(N-1)
 for m in range(1,N):
  for v in range(n):
   if m>>v&1:
    sub=m^(1<<v)
    end[m][v]=sum(end[sub][u]for u in range(n)if r[u]>>v&1)if sub else 1
    start[m][v]=sum(start[sub][u]for u in range(n)if r[v]>>u&1)if sub else 1
  h[m]=sum(end[m])
 fs=[sum(end[m][v]*h[(N-1)^m]for m in range(1,N))for v in range(n)];ports=[sn.index(v)for v in c['in_ports']]
 alljoints=[dict(end=sn[u],start=sn[v],value=sum(end[m][u]*start[(N-1)^m][v]for m in range(1,N)))for u in range(n) for v in range(n)]
 joints=[z for z in alljoints if z['end'] in c['in_ports'] and z['start'] in c['in_ports']]
 for u in range(n):
  if fs[u]-h[-1] != sum(z['value'] for z in alljoints if z['end']==sn[u] and r[sn.index(z['start'])]>>u&1):raise RuntimeError('Inner defect transport identity')
 rootidx=ns.index(c['root'])
 if set(c['in_ports'])!={v for v in sn if c['rows'][ns.index(v)]>>rootidx&1}:raise RuntimeError('Physical root ports mismatch')
 for i,a in enumerate(c['prefix']):
  if not (c['rows'][rootidx]>>ns.index(a)&1):raise RuntimeError('Root-prefix edge')
  for b in c['prefix'][i+1:]+sn:
   if not (c['rows'][ns.index(a)]>>ns.index(b)&1):raise RuntimeError('Prefix or cross edge')
 F=sum(fs[u]for u in ports);J=sum(z['value']for z in joints);H=h[-1]-J+(1<<len(c['prefix']))*F
 if (h[-1],F,J,H)!=(c['h'],c['F'],c['J'],c['target_H']):raise RuntimeError('Candidate evaluator mismatch')
 # Full tournament orientation, source preservation except recorded edge flips.
 src=json.load(open(D/'source.json'));srcnames=src['names'];expected={frozenset(p)for p in c['all_changed_pairs']};actual=set()
 for i,a in enumerate(srcnames):
  for j,b in enumerate(srcnames):
   if i<j and ((src['rows'][i]>>j)&1)!=((c['rows'][ns.index(a)]>>ns.index(b))&1):actual.add(frozenset((a,b)))
 if actual!=expected:raise RuntimeError('Edge changes mismatch')
 for i in range(len(ns)):
  if c['rows'][i]>>i&1:raise RuntimeError('Loop')
  for j in range(i):
   if ((c['rows'][i]>>j)&1)+((c['rows'][j]>>i)&1)!=1:raise RuntimeError('Tournament')
  reach=1<<i
  while True:
   nxt=reach
   for j in range(len(ns)):
    if reach>>j&1:nxt|=c['rows'][j]
   if nxt==reach:break
   reach=nxt
  if reach!=(1<<len(ns))-1:raise RuntimeError('Not strong')
 counts=collections.Counter(z['s']for z in map(json.loads,(D/'hits.jsonl').read_text().splitlines())if 's' in z)
 checks.append(dict(digit=digit,H=H,h=h[-1],F=F,J=J,c=h[-1]-J,root_E2=h[-1]+(1<<len(c['prefix']))*F,port_responses={str(sn[u]):fs[u]for u in ports},joint=joints,all_joint=alljoints,all_port_responses={str(sn[u]):fs[u]for u in range(n)},inner_defects={str(sn[u]):fs[u]-h[-1]for u in range(n)},hits_by_clock=dict(counts),all_old_names_preserved=True,source_H=src['target_H'],candidate_sha256=hashlib.sha256((D/f'candidate{digit}.json').read_bytes()).hexdigest()))
(D/'verification.json').write_text(json.dumps(dict(status='PASS',candidates=checks),indent=2)+'\n');print(json.dumps(checks,indent=2))
