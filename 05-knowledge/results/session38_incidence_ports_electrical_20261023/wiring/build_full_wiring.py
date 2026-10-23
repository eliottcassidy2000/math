#!/usr/bin/env python3
"""Certify and enumerate the entire order-eight complement-line metagraph.
No isomorphism oracle, numerical solver, or removable assert is used.
"""
import argparse,array,csv,gzip,hashlib,io,itertools,json,math,struct,time
from collections import Counter,defaultdict,deque
from pathlib import Path
HERE=Path(__file__).resolve().parent
N=8; PAIRS=list(itertools.combinations(range(N),2))
TILES=[(b,a) for b in range(N-2) for a in range(b+2,N)]
FULL=(1<<len(TILES))-1; CLOSE=TILES.index((0,N-1))
REF=[TILES.index((N-1-a,N-1-b)) for b,a in TILES]
STARTBITS=sum(1<<k for k,(b,a) in enumerate(TILES) if a==N-1)
ENDBITS=sum(1<<k for k,(b,a) in enumerate(TILES) if b==0)
checks=Counter(); WRITE=False
def need(ok,why):
 checks[why]+=1
 if not ok:raise ValueError(why)
def decode(code):
 rows=[0]*N
 for k,(i,j) in enumerate(PAIRS):
  u,v=(i,j) if code>>(27-k)&1 else(j,i);rows[u]|=1<<v
 return rows
def canonical(rows):
 scores=[r.bit_count() for r in rows]
 blocks=[[i for i in range(N) if scores[i]==s] for s in sorted(set(scores))]
 best=None;aut=0;count=0
 for chunks in itertools.product(*(itertools.permutations(b) for b in blocks)):
  p=sum(chunks,());value=0
  for i,j in PAIRS:value=2*value+((rows[p[i]]>>p[j])&1)
  count+=1
  if best is None or value<best:best=value;aut=1
  elif value==best:aut+=1
 return best,aut,count

def write_gzip(path,data):
 out=io.BytesIO()
 with gzip.GzipFile(filename='',fileobj=out,mode='wb',mtime=0,compresslevel=9) as z:z.write(data)
 compressed=out.getvalue()
 if WRITE:path.write_bytes(compressed)
 need(path.read_bytes()==compressed,'frozen artifact bytes match')
def tsv_gzip(name,header,rows):
 raw=('\t'.join(header)+'\n'+''.join('\t'.join(map(str,row))+'\n' for row in rows)).encode()
 write_gzip(HERE/name,raw)
 return hashlib.sha256(raw).hexdigest()
def components(count,es):
 adj=[[] for _ in range(count)]
 for (u,v),w in es.items():
  if u!=v and w:adj[u].append(v);adj[v].append(u)
 seen=set();out=[]
 for v in range(count):
  if v in seen:continue
  seen.add(v);q=[v];comp=[]
  while q:
   u=q.pop();comp.append(u)
   for w in adj[u]:
    if w not in seen:seen.add(w);q.append(w)
  out.append(sorted(comp))
 return out,adj

def graph_inventory(count,es):
 cs,adj=components(count,es)
 # Iterative Tarjan, avoiding recursion depth problems.
 disc=[-1]*count;low=[0]*count;par=[-1]*count;timer=0;bridges=[]
 for root in range(count):
  if disc[root]>=0:continue
  disc[root]=low[root]=timer;timer+=1;stack=[(root,0)]
  while stack:
   u,k=stack[-1]
   if k<len(adj[u]):
    v=adj[u][k];stack[-1]=(u,k+1)
    if disc[v]<0:par[v]=u;disc[v]=low[v]=timer;timer+=1;stack.append((v,0))
    elif v!=par[u]:low[u]=min(low[u],disc[v])
   else:
    stack.pop()
    if par[u]>=0:
     p=par[u];low[p]=min(low[p],low[u])
     if low[u]>disc[p]:bridges.append((min(p,u),max(p,u)))
 support=sum(u!=v and w>0 for (u,v),w in es.items())
 loops=sum(w for(u,v),w in es.items() if u==v)
 weighted=[0]*count
 for(u,v),w in es.items():weighted[u]+=w;weighted[v]+=w
 return dict(vertices=count,support_leaves=[i for i,a in enumerate(adj) if len(a)==1],multigraph_leaves=[i for i,d in enumerate(weighted) if d==1],components=len(cs),component_sizes=sorted(map(len,cs),reverse=True),isolated_vertices=[c[0] for c in cs if len(c)==1],nonloop_support_edges=support,nonloop_lines=sum(w for(u,v),w in es.items() if u!=v),loop_support_nodes=sum(u==v and w>0 for(u,v),w in es.items()),loop_lines=loops,support_cycle_rank=support-count+len(cs),multigraph_cycle_rank=sum(es.values())-count+len(cs),support_bridges=[dict(source=u,target=v,multiplicity=es[u,v]) for u,v in sorted(bridges)],multigraph_bridges=[list(e) for e in sorted(bridges) if es[e]==1]),cs

def main():
 global WRITE
 parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args();WRITE=args.write
 start=time.time();inputs=list(csv.DictReader((HERE/'classes_input_n8.tsv').open(),delimiter='\t'))
 inputs.sort(key=lambda r:int(r['canon']));canons=[int(r['canon']) for r in inputs];lookup={c:i for i,c in enumerate(canons)}
 need(len(inputs)==6880,'ordinary class count');need(len(lookup)==6880,'unique canonical codes')
 mats=[decode(c) for c in canons];conv=[];auts=[];perms=0;orbitmass=0
 for ci,(r,A) in enumerate(zip(inputs,mats)):
  code,aut,np=canonical(A);perms+=np;need(code==int(r['canon']),'score canonical minimum');need(aut==int(r['Aut']),'automorphism agreement');auts.append(aut);orbitmass+=math.factorial(N)//aut
  ca,_,_=canonical([((1<<N)-1)^A[i]^(1<<i) for i in range(N)])
  need(ca in lookup,'converse representative exists');conv.append(lookup[ca])
  need((ca==code)==bool(int(r['SC'])),'self-converse agreement')
 need(orbitmass==1<<28,'full labelled orbit mass')
 need(all(conv[conv[i]]==i for i in range(len(conv))),'converse involution')
 print('Canonical certificates:',len(inputs),'classes;',perms,'score-block permutations;',orbitmass,'labelled tournaments',flush=True)
 cid=array.array('H',[65535])*(FULL+1);occ=bytearray(FULL+1);Hs=[];HC=[];fibs=[]
 for ci,(r,A) in enumerate(zip(inputs,mats)):
  path=[];h=0;closed=0;fiber=0
  def visit(v,unused):
   nonlocal h,closed,fiber
   path.append(v)
   if not unused:
    h+=1;closed+=(A[v]>>path[0])&1;p=path[::-1];mask=0
    for k,(b,a) in enumerate(TILES):mask|=((A[p[b]]>>p[a])&1)<<k
    if cid[mask]==65535:cid[mask]=ci;fiber+=1
    else:need(cid[mask]==ci,'Hamilton path class uniqueness')
    occ[mask]+=1
   else:
    nxt=A[v]&unused
    while nxt:
     bit=nxt&-nxt;nxt^=bit;visit(bit.bit_length()-1,unused^bit)
   path.pop()
  for v in range(N):visit(v,((1<<N)-1)^(1<<v))
  need(h==int(r['H']),'Hamilton path count agreement');need(closed%N==0,'cyclic path divisibility');need(h==fiber*auts[ci],'ordinary class Hamilton fiber')
  need(fiber%2==1,'odd ordinary Hamilton fiber');Hs.append(h);HC.append(closed//N);fibs.append(fiber)
 need(65535 not in cid,'all fixed-path masks covered');need(sum(fibs)==FULL+1,'complete fixed-path fiber mass')
 need(all(occ[m]==auts[cid[m]] for m in range(FULL+1)),'every mask has exact automorphism multiplicity')
 need(all(fibs[i]==fibs[conv[i]] for i in range(len(inputs))),'converse fiber agreement')
 print('Path enumeration:',sum(Hs),'Hamilton paths;',sum(fibs),'distinct masks;',round(time.time()-start,2),'seconds',flush=True)
 # Mirror lookup by independent chunks. Every free arc index is explicitly reflected.
 mirrorlo=[sum(((x>>k)&1)<<REF[k] for k in range(11)) for x in range(1<<11)]
 mirrorhi=[sum(((x>>k)&1)<<REF[k+11] for k in range(10)) for x in range(1<<10)]
 nodekeys=sorted(set(min(i,conv[i]) for i in range(len(inputs))));nodeix={key:i for i,key in enumerate(nodekeys)};nodeof=[nodeix[min(i,conv[i])] for i in range(len(inputs))];nv=len(nodekeys)
 bluefiber=Counter();closedfiber=Counter();endports=Counter();openports=Counter();closedports=Counter();edges=defaultdict(lambda:[0,0]);ports=Counter();unoriented={'blue':Counter(),'black':Counter(),'paired_black':Counter(),'all':Counter(),'fundamental':Counter()};current=Counter();c3=[math.comb(N,3)-sum(math.comb(row.bit_count(),2) for row in A) for A in mats]
 for mask in range(FULL+1):
  ci=cid[mask];cj=cid[mask^FULL];mirror=mirrorlo[mask&2047]|mirrorhi[mask>>11];blue=mask==mirror
  need(cid[mirror]==conv[ci],'all-mask reflection converse transport')
  a=N-1-(mask&STARTBITS).bit_count();b=(mask&ENDBITS).bit_count();isclosed=(mask>>CLOSE)&1
  bluefiber[ci]+=blue;endports[ci,a,b]+=1
  need(c3[cj]-c3[ci]==a-b-1,'all-mask triangle current law')
  need((N-1-((mask^FULL)&STARTBITS).bit_count(),((mask^FULL)&ENDBITS).bit_count())==(N-a,N-2-b),'all-mask endpoint transport')
  current[nodeof[ci]]+=a-b-1
  if not isclosed:openports[ci,a,b,blue]+=1;continue
  closedfiber[ci]+=1;closedports[ci,a,b,blue]+=1;u,v=nodeof[ci],nodeof[cj];edges[u,v][0 if blue else 1]+=1;ports[ci,cj,a,b,int(blue)]+=1
 for ci,r in enumerate(inputs):
  need(bluefiber[ci]==int(r['bm']),'blue class fiber agreement');need(closedfiber[ci]*auts[ci]==N*HC[ci],'closed Hamilton-cycle fiber')
  need(not bluefiber[ci] or conv[ci]==ci,'blue only self-converse');need(bluefiber[ci]%2==int(conv[ci]==ci),'ordinary blue parity')
 for (u,v),(blue,black) in edges.items():
  need(black%2==0,'black directed twins including loops');e=min(u,v),max(u,v)
  for name,w in [('blue',blue),('black',black),('paired_black',black//2),('all',blue+black),('fundamental',blue+black//2)]:
   if w:unoriented[name][e]+=w
 need(sum(bluefiber.values())==4096,'blue tiling count');need(sum(v[0] for v in edges.values())==2048,'blue line count');need(sum(v[1] for v in edges.values())==1046528,'black line count');need(sum(unoriented['paired_black'].values())==523264,'paired black line count')
 need(sum(ports.values())==1<<20,'closed port mass');need(sum(openports.values())==1<<20,'open port mass')
 # All degree identities count a loop twice. Never substitute simple-support degree.
 degrees={}
 for name,es in unoriented.items():
  d=[0]*nv
  for(u,v),w in es.items():d[u]+=w;d[v]+=w
  degrees[name]=d
 for u,ci in enumerate(nodekeys):
  t=fibs[ci]*(1 if conv[ci]==ci else 2);s=bluefiber[ci]
  need(degrees['all'][u]==t,'full degree equals merged tiling fiber')
  need(degrees['blue'][u]==s,'blue degree equals fixed tiling fiber')
  need(degrees['paired_black'][u]==(t-s)//2,'paired-black degree equals free reflection orbits')
  need(degrees['fundamental'][u]==(t+s)//2,'fundamental degree identity')
 # Literal actual-port conservation and reflection at ordinary class level.
 targetm=Counter();sourcem=Counter()
 for (ci,cj,a,b,blue),w in ports.items():
  targetm[cj,N-a,N-2-b,bool(blue)]+=w;sourcem[ci,a,b,bool(blue)]+=w
  need(ports[conv[ci],conv[cj],N-1-b,N-1-a,blue]==w,'ordinary port converse-reflection')
 need(targetm==openports,'all target port capacities');need(sourcem==closedports,'all source port capacities')
 nodec3=[c3[i] for i in nodekeys];lapcurrent=Counter()
 for(u,v),w in unoriented['all'].items():
  lapcurrent[u]+=w*(nodec3[v]-nodec3[u]);lapcurrent[v]+=w*(nodec3[u]-nodec3[v])
 need(all(lapcurrent[u]==current[u] for u in range(nv)),'integer merged Laplacian current');need(sum(current.values())==0,'total current cancellation')
 inventories={};compmaps={}
 for name,es in unoriented.items():inventories[name],compmaps[name]=graph_inventory(nv,es)
 # Every support bridge is a true multigraph bridge exactly at multiplicity1.
 scnodes=[nodeof[i] for i in range(len(inputs)) if conv[i]==i]
 blueactive=[c for c in compmaps['blue'] if len(c)>1 or unoriented['blue'].get((c[0],c[0]),0)]
 print('Full graph:',nv,'merged nodes;',len(edges),'directed endpoint pairs;',len(ports),'ordinary port entries;',len(unoriented['all']),'undirected support pairs incl loops;',inventories['all']['components'],'components',flush=True)
 hashes={}
 raw=struct.pack('<'+str(len(cid))+'H',*cid);write_gzip(HERE/'class_of_mask.u16le.gz',raw);hashes['class_of_mask.u16le']=hashlib.sha256(raw).hexdigest()
 hashes['classes.tsv']=tsv_gzip('classes.tsv.gz',['class_id','canon','input_class_id','converse_class','node_id','H','Aut','fiber','blue_fiber','closed_fiber','HC','c3','scores'],((ci,canons[ci],inputs[ci]['class'],conv[ci],nodeof[ci],Hs[ci],auts[ci],fibs[ci],bluefiber[ci],closedfiber[ci],HC[ci],c3[ci],','.join(map(str,sorted(r.bit_count() for r in mats[ci])))) for ci in range(len(inputs))))
 hashes['nodes.tsv']=tsv_gzip('nodes.tsv.gz',['node_id','canon','class_id','converse_class','SC','H','Aut','tilings','blue_tilings','closed_tilings','c3','current_c3'],((u,canons[ci],ci,conv[ci],int(conv[ci]==ci),Hs[ci],auts[ci],fibs[ci]*(1 if conv[ci]==ci else 2),bluefiber[ci],closedfiber[ci]*(1 if conv[ci]==ci else 2),c3[ci],current[u]) for u,ci in enumerate(nodekeys)))
 hashes['edges.tsv']=tsv_gzip('edges.tsv.gz',['source','target','blue','black'],((u,v,*edges[u,v]) for u,v in sorted(edges)))
 hashes['ports.tsv']=tsv_gzip('ports.tsv.gz',['source_class','target_class','start_score','end_score','blue','multiplicity'],((*key,w) for key,w in sorted(ports.items())))
 hashes['endpoint_capacities.tsv']=tsv_gzip('endpoint_capacities.tsv.gz',['class_id','start_score','end_score','blue','closed','open'],((*key,closedports[key],openports[key]) for key in sorted(set(closedports)|set(openports))))
 hashes['components.tsv']=tsv_gzip('components.tsv.gz',['graph','component_id','node_id'],((name,k,u) for name,cs in compmaps.items() for k,c in enumerate(cs) for u in c))
 receipt=dict(scope='FINITE-EXACT: every order-eight tournament class and all 2^21 fixed-path presentations; converse-merged complement lines, including loops and black twins.',n=N,ordinary_classes=len(inputs),self_converse_classes=len(scnodes),merged_nodes=nv,score_block_permutations=perms,labelled_orbit_mass=orbitmass,hamilton_paths_enumerated=sum(Hs),fixed_path_presentations=FULL+1,blue_presentations=sum(bluefiber.values()),closed_to_open_lines=1<<20,blue_lines=2048,black_lines=1046528,paired_black_lines=523264,directed_merged_entries=len(edges),ordinary_port_entries=len(ports),endpoint_capacity_entries=len(set(closedports)|set(openports)),blue_active_components=[len(c) for c in blueactive],graph_inventories=inventories,checks=dict(sorted(checks.items())),uncompressed_sha256=hashes,input_sha256=hashlib.sha256((HERE/'classes_input_n8.tsv').read_bytes()).hexdigest(),file_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob('*.gz'))})
 if args.write:(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
 else:need(receipt==json.loads((HERE/'receipt.json').read_text()),'frozen receipt reproduction')
 print('DONE',round(time.time()-start,2),'seconds; receipt',('written' if args.write else 'matched'),flush=True)
if __name__=='__main__':main()
