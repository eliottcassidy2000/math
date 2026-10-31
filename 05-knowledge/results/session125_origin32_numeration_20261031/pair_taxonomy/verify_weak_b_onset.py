#!/usr/bin/env python3
"""Independent degree-subset solver for weak-B pair onset; no external data."""
import argparse,itertools,json,math,pathlib

def graph(n):
 return [(a,b) for a in range(1,n+1) for b in range(a+1,n+1) if math.isqrt(a+b)**2==a+b]

def enumerate_by_degrees(n,cycle=False):
 edges=graph(n); incident=[0]*(n+1)
 for i,(a,b) in enumerate(edges):incident[a]|=1<<i;incident[b]|=1<<i
 goal=10 if n%2 else 8
 target=[0]+[2]*n
 if not cycle:target[1]=target[goal]=1
 maskpaths={}; solutions=0; components={}; nodes=0; allcycles=[]; raw_count=0; port_counts={}; raw_port_witnesses={}
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
  pk=str(int(bool(mask&4)))+str(int(bool(mask&128)))+str(int(bool(mask&512)));port_counts[pk]=port_counts.get(pk,0)+1
  raw_port_witnesses.setdefault(pk,path[:])
  if mask&4 or mask&128 or not(mask&512):return
  entry=maskpaths.setdefault(mask,{'count':0,'path':path});entry['count']+=1
 search(0 if cycle else 1<<edges.index((1,3)),(1<<len(edges))-1)
 return dict(n=n,edges=len(edges),degree_subgraphs=solutions,component_counts=components,nodes=nodes,unfiltered_path_count=raw_count,port_color_counts=port_counts,raw_port_witnesses=raw_port_witnesses,path_count=sum(e['count'] for e in maskpaths.values()),mask_count=len(maskpaths),masks=sorted(maskpaths),witnesses=maskpaths,cycles=allcycles)

def verify_path(path,n,odd_mask=None):
 if sorted(path)!=list(range(1,n+1)):raise ValueError('not permutation')
 if path[:2]!=[1,3] or path[-1]!=(10 if n%2 else 8):raise ValueError('endpoints/prefix')
 if not all(math.isqrt(a+b)**2==a+b for a,b in zip(path,path[1:])):raise ValueError('non-square')
 actual=sum(1<<(v-1) for v in path[::2])
 if odd_mask is not None and actual!=odd_mask:raise ValueError('mask')
 if actual&4 or actual&128 or not(actual&512):raise ValueError('ports')
 return actual


def enumerate_by_paths(n,ports=True):
 # Independent path-growth search; graph constructed by enumerating squares.
 adj=[0]*(n+1)
 for k in range(2,math.isqrt(2*n)+1):
  for a in range(1,n+1):
   b=k*k-a
   if 1<=b<=n and a!=b:adj[a]|=1<<(b-1)
 goal=10 if n%2 else 8
 path=[1,3];records={};calls=0
 def dfs(x,unused,mask):
  nonlocal calls
  calls+=1
  if not unused:
   if x==goal:
    entry=records.setdefault(mask,{'count':0,'path':path[:]});entry['count']+=1
   return
  if x==goal:return
  available=unused|(1<<(x-1));scan=unused
  while scan:
   bit=scan&-scan;scan-=bit;v=bit.bit_length()
   if (adj[v]&available).bit_count()<(1 if v==goal else 2):return
  seen=1<<(x-1);front=seen
  while front:
   bit=front&-front;front-=bit;v=bit.bit_length()
   fresh=adj[v]&available&~seen;seen|=fresh;front|=fresh
  if seen!=available:return
  opts=adj[x]&unused;pos=len(path)+1
  while opts:
   bit=opts&-opts;opts-=bit;v=bit.bit_length()
   if ports and ((v==8 and pos%2) or (v==10 and not pos%2)):continue
   path.append(v);dfs(v,unused^bit,mask|(bit if pos%2 else 0));path.pop()
 dfs(3,((1<<n)-1)^1^4,1)
 return records,calls

def require(claim,message):
 if not claim:raise RuntimeError(message)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output',default='weak_b_certificate.json');args=parser.parse_args()
 reports=[];allrecords={};last_boundary=None
 for n in range(32,44):
  degree=enumerate_by_degrees(n);paths,calls=enumerate_by_paths(n)
  require(set(paths)==set(degree['witnesses']),'phase masks disagree at '+str(n))
  for m,p in paths.items():
   require(p['count']==degree['witnesses'][m]['count'],'multiplicity differs')
   require(verify_path(p['path'],n)==m,'path verifier failed')
   require(verify_path(degree['witnesses'][m]['path'],n)==m,'degree witness verifier failed')
  allrecords[n]=degree['witnesses']
  if n==41:last_boundary=degree['raw_port_witnesses']
  reports.append({'n':n,'unfiltered_prefix_endpoint_paths':degree['unfiltered_path_count'],'eligible_paths':degree['path_count'],'phase_masks':degree['mask_count'],'degree_subgraphs':degree['degree_subgraphs'],'degree_nodes':degree['nodes'],'path_nodes':calls,'port_color_counts':degree['port_color_counts'],'phase_counts':[[m,paths[m]['count']] for m in sorted(paths)]})
  print('AGREE',n,degree['path_count'],degree['mask_count'],flush=True)
 pairs=[]
 for n in range(32,43):
  a=allrecords[n];b=allrecords[n+1];mask=(1<<n)-1;common=[]
  for m,q in sorted(b.items()):
   if m&mask in a:
    p=a[m&mask];common.append({'common_odd_mask':m&mask,'shorter_count':p['count'],'longer_count':q['count'],'shorter':p['path'],'longer':q['path']})
  pairs.append({'n':n,'phase_masks':len(common),'path_pairs':sum(p['shorter_count']*p['longer_count'] for p in common),'witnesses':common})
 require(all(p['phase_masks']==0 for p in pairs[:-1]),'earlier pair')
 require((pairs[-1]['n'],pairs[-1]['phase_masks'],pairs[-1]['path_pairs'])==(42,6,20),'onset42 count')
 require(not allrecords[32],'B32 should be empty')
 require(set(last_boundary)=={'011'},'B41 last boundary port colors')
 last41=last_boundary['011']
 raw41,_=enumerate_by_paths(41,ports=False)
 require(sum(x['count'] for x in raw41.values())==1,'independent raw41 uniqueness')
 require(next(iter(raw41.values()))['path']==last41,'independent raw41 witness')
 require(last41[:2]==[1,3] and last41[-1]==10,'raw B41 normalization')
 require(sorted(last41)==list(range(1,42)),'raw B41 permutation')
 require(all(math.isqrt(a+b)**2==a+b for a,b in zip(last41,last41[1:])),'raw B41 squares')
 require((last41.index(8)+1)%2==1,'raw B41 should put8 in odd slot')
 # Closing a B32 path by1--8 would be a Hamiltonian cycle containing1--3.
 c32=enumerate_by_degrees(32,cycle=True)
 require(len(c32['cycles'])==1,'unique C32 mismatch')
 cyc=c32['cycles'][0]
 require(all({a,b}!={1,3} for a,b in zip(cyc,cyc[1:]+cyc[:1])),'C32 contains1--3')
 p=pairs[-1]['witnesses'][0]['shorter'];q=pairs[-1]['witnesses'][0]['longer']
 require(verify_path(p,42)==verify_path(q,43)&((1<<42)-1),'pair parity')
 require((q.index(43)+1)%2==1,'new43 not odd position')
 # Hostile tests: numerical parity and original-nice port coloring are not substitutes.
 require((verify_path(p,42)&4)==0,'weak B label3 is not even-positioned')
 require((verify_path(p,42)&512)!=0,'weak B label10 is not odd-positioned')
 hostile={}
 for label,mutated in [('duplicate',p[:1]+[1]+p[2:]),('wrong_prefix',p[::-1]),('nonsquare',p[:2]+[p[3],p[2]]+p[4:])]:
  try:verify_path(mutated,42)
  except ValueError:hostile[label]=True
  else:hostile[label]=False
 require(all(hostile.values()),'hostile controls')
 out={'status':'FINITE-EXACT','definition':'weak B: prefix1,3; last8 even/10 odd; common positional parity','first_pair_n':42,'first_pair_count':20,'first_pair_phase_masks':6,'orders':reports,'pairs':pairs,'C32':cyc,'last_boundary_41':{'unique_prefix1_3_endpoint10_path':last41,'position_of8':last41.index(8)+1,'obstruction':'8 is odd-positioned but any even-member endpoint8 is even-positioned'},'first_new43_position':q.index(43)+1,'hostile_controls':hostile}
 pathlib.Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
 print('PASS first weak-B pair42/43,20 pairs,6 phases')
if __name__=='__main__':main()
