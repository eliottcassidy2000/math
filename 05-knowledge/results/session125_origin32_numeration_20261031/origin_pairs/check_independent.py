#!/usr/bin/env python3
"""Independent degree-subset solver for original nice-pair onset; no external data."""
import argparse,base64,gzip,hashlib,itertools,json,math,pathlib

DECODED_CERTIFICATE_SHA256 = 'bef547efacc1ae2a4365e44a22200d5028db327262d518c74395cf0deb6fe720'

def read_reference():
 """Read our compressed certificate and authenticate its exact decoded bytes."""
 here=pathlib.Path(__file__).resolve().parent
 encoded=(here/'independent_certificate.json.gz.b64').read_bytes()
 raw=gzip.decompress(base64.b64decode(b''.join(encoded.split()),validate=True))
 if hashlib.sha256(raw).hexdigest()!=DECODED_CERTIFICATE_SHA256:
  raise RuntimeError('decoded certificate SHA256 mismatch')
 return raw

def graph(n):
 return [(a,b) for a in range(1,n+1) for b in range(a+1,n+1) if math.isqrt(a+b)**2==a+b]

def enumerate_by_degrees(n,cycle=False):
 edges=graph(n); incident=[0]*(n+1)
 for i,(a,b) in enumerate(edges):incident[a]|=1<<i;incident[b]|=1<<i
 goal=3 if n%2 else 8
 target=[0]+[2]*n
 if not cycle:target[1]=target[goal]=1
 maskpaths={}; solutions=0; components={}; nodes=0; allcycles=[]; raw_count=0; port_counts={}
 def search(selected,allowed):
  nonlocal solutions,nodes,raw_count
  nodes+=1
  while True:
   old=(selected,allowed)
   for v in range(1,n+1):
    have=(incident[v]&selected).bit_count();avail=(incident[v]&allowed).bit_count()
    if have>target[v] or avail<target[v]:return
    if have==target[v]:allowed&=~(incident[v]&~selected)
    if avail==target[v]:selected|=incident[v]&allowed
   if (selected,allowed)==old:break
  optional=allowed&~selected
  if optional:
   choices=[]
   for v in range(1,n+1):
    opts=incident[v]&optional
    if opts:choices.append((opts.bit_count(),v,opts))
   _,v,opts=min(choices)
   need=target[v]-(incident[v]&selected).bit_count()
   ids=[i for i in range(len(edges)) if opts>>i&1]
   for comb in itertools.combinations(ids,need):
    choose=sum(1<<i for i in comb)
    search(selected|choose,allowed&~(opts&~choose))
   return
  solutions+=1
  adj={v:[] for v in range(1,n+1)}
  for i,(a,b) in enumerate(edges):
   if selected>>i&1:adj[a].append(b);adj[b].append(a)
  sizes=[];seen=set()
  for root in range(1,n+1):
   if root in seen:continue
   stack=[root];seen.add(root);sz=0
   while stack:
    x=stack.pop();sz+=1
    for y in adj[x]:
     if y not in seen:seen.add(y);stack.append(y)
   sizes.append(sz)
  key='+'.join(map(str,sorted(sizes)));components[key]=components.get(key,0)+1
  if len(sizes)>1:return
  path=[1];previous=0;current=1
  while len(path)<n:
   nxt=next(x for x in adj[current] if x!=previous)
   path.append(nxt);previous,current=current,nxt
  if cycle:allcycles.append(path);return
  if path[-1]!=goal:raise RuntimeError('wrong terminal')
  mask=sum(1<<(v-1) for v in path[::2])
  raw_count+=1
  pk=str(int(bool(mask&4)))+str(int(bool(mask&128)));port_counts[pk]=port_counts.get(pk,0)+1
  if not(mask&4) or mask&128:return
  entry=maskpaths.setdefault(mask,{'count':0,'path':path});entry['count']+=1
 search(0,(1<<len(edges))-1)
 return dict(n=n,edges=len(edges),degree_subgraphs=solutions,component_counts=components,nodes=nodes,unfiltered_path_count=raw_count,port_color_counts=port_counts,path_count=sum(e['count'] for e in maskpaths.values()),mask_count=len(maskpaths),masks=sorted(maskpaths),witnesses=maskpaths,cycles=allcycles)

def verify_path(path,n,odd_mask=None):
 if sorted(path)!=list(range(1,n+1)):raise ValueError('not permutation')
 if path[0]!=1 or path[-1]!=(3 if n%2 else 8):raise ValueError('endpoints')
 if not all(math.isqrt(a+b)**2==a+b for a,b in zip(path,path[1:])):raise ValueError('non-square')
 actual=sum(1<<(v-1) for v in path[::2])
 if odd_mask is not None and actual!=odd_mask:raise ValueError('mask')
 if not(actual&4) or actual&128:raise ValueError('ports')
 return actual

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output',default='independent_certificate.json');args=parser.parse_args()
 here=pathlib.Path(__file__).resolve().parent
 producer=[json.loads(l) for l in (here/'path_census.jsonl').read_text().splitlines()]
 results=[]
 for src in producer:
  n=src['n'];r=enumerate_by_degrees(n);print(n,r['path_count'],r['mask_count'],r['degree_subgraphs'],flush=True)
  for p in src['witnesses']:
   verify_path(p['path'],n,p['odd_mask'])
   if r['witnesses'].get(p['odd_mask'],{}).get('count')!=p['count']:raise RuntimeError('mask multiplicity mismatch')
  if r['path_count']!=src['path_count'] or r['masks']!=sorted(p['odd_mask'] for p in src['witnesses']):raise RuntimeError('independent mismatch at '+str(n))
  results.append(r)
 compatibility=[]
 for a,b in zip(results,results[1:]):
  n=a['n'];mb={x&((1<<n)-1):x for x in b['masks']};common=sorted(set(a['masks'])&set(mb));pairs=[]
  for m in common:
   pa=a['witnesses'][m];pb=b['witnesses'][mb[m]]
   pairs.append(dict(shorter=pa['path'],longer=pb['path'],odd_mask_common=m,shorter_count=pa['count'],longer_count=pb['count']))
  compatibility.append(dict(n=n,common_mask_count=len(common),path_pair_count=sum(p['shorter_count']*p['longer_count'] for p in pairs),witnesses=pairs))
 cycles=[enumerate_by_degrees(n,cycle=True) for n in (31,32,33)]
 for r in cycles:print('CYCLES',r['n'],r['component_counts'],flush=True)
 out=dict(status='FINITE-EXACT',method='independent exhaustive degree subsets vs endpoint DFS',path_orders=results,nice_pairs=compatibility,cycles=cycles)
 raw=(json.dumps(out,indent=2)+'\n').encode('utf-8')
 if raw!=read_reference():raise RuntimeError('recomputed certificate differs from pinned reference')
 pathlib.Path(args.output).write_bytes(raw)
 print('PASS',args.output)
if __name__=='__main__':main()
