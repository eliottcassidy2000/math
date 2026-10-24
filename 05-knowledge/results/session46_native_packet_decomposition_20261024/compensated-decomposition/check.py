#!/usr/bin/env python3
"""Independent exact native, resource, row-current, image, and bottleneck audit.
No assert statements; works unchanged under python -O. Search output is not trusted.
"""
from pathlib import Path
from itertools import combinations,permutations
from collections import Counter,defaultdict,deque
from array import array
import argparse,csv,gzip,hashlib,json,sys,heapq
H=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--out',type=Path);ap.add_argument('--include-eighteen',action='store_true');ap.add_argument('--include-eighteen37',action='store_true');ap.add_argument('--include-variants',action='store_true');args=ap.parse_args()
if args.include_variants:args.include_eighteen=True;args.include_eighteen37=True
R=args.repo
checks=0
def need(v,m):
 global checks
 checks+=1
 if not v:raise RuntimeError(m)
def sha(b):return hashlib.sha256(b).hexdigest()
pins=json.loads((H/'SOURCE_PINS.json').read_text());blobs={}
for p,h in pins['originals'].items():
 b=(R/p).read_bytes();need(sha(b)==h,'original data hash '+p);blobs[Path(p).name]=b
nodes=list(csv.DictReader(blobs['merged_nodes.tsv'].decode().splitlines(),delimiter='\t'));N=len(nodes);need(N==3528,'original node count')
owner=array('H');owner.frombytes(gzip.decompress(blobs['owner_map.u16.gz']))
if sys.byteorder!='little':owner.byteswap()
merged={int(r[k]):int(r['merged_id']) for r in nodes for k in ['class_a','class_b']}
edges=[tuple(int(r[k]) for k in ['source','target','blue','black']) for r in csv.DictReader(gzip.decompress(blobs['merged_edges.tsv.gz']).decode().splitlines(),delimiter='\t')];colored={(i,j):(b,k) for i,j,b,k in edges}
D=json.loads(blobs['receipt.json'])['order8'];need(D['complete_two_primary'] and D['higher_layers'][-1]['level']==8 and D['higher_layers'][-1]['next_dimension']==0,'inherited certified dyadic exponent saturation at256');X=D['higher_layers'][-1]['lifts'][0];Y=D['radical_lifts_mod4'];chars=[[64*y%256 for y in row] for row in Y[:4]]+[X]
def lap(v):
 out=[0]*N
 for i,j,b,k in edges:d=(b+k)*(v[i]-v[j]);out[i]+=d;out[j]-=d
 return out
for chi in chars:need(all(z%256==0 for z in lap(chi)),'authentic higher-dyadic character')
U,V=1965,2234
need(X[U]==X[V]==146,'literal chosen target phase')
need(all(2*(chi[U]-chi[V])%256==0 for chi in chars),'paired target charge equality')
P=json.loads((H/'packet.json').read_text());need(P['targets']==[U,V],'target identifiers')
if args.include_eighteen:P['routes'].extend(json.loads((H/'packet_eighteen.json').read_text())['routes'])
if args.include_eighteen37:P['routes'].extend(json.loads((H/'packet_eighteen37.json').read_text())['routes'])
if args.include_variants:P['routes'].extend(json.loads((H/'packet_variants.json').read_text())['routes'])
pairs=list(combinations(range(8),2));chords=[(i,j) for i,j in pairs if j>=i+2];indices={e:i for i,e in enumerate(chords)};full=(1<<21)-1

def matrix(code):
 A=[[0]*8 for _ in range(8)]
 for bit,(i,j) in enumerate(pairs):A[i][j]=code>>(27-bit)&1;A[j][i]=1-A[i][j]
 return A

def word(A,p):return sum(A[p[i]][p[j]]<<b for b,(i,j) in enumerate(chords))
def reflect(w):return sum(((w>>b)&1)<<indices[7-j,7-i] for b,(i,j) in enumerate(chords))
def orbit(w):return {w,w^full,reflect(w),reflect(w)^full}
def transpose(A):return list(map(list,zip(*A)))
def partner(A,p):
 E={frozenset(e) for e in zip(p,p[1:])}
 return [[0 if i==j else A[i][j] if frozenset((i,j)) in E else A[j][i] for j in range(8)] for i in range(8)]
def autos(A):
 deg=list(map(sum,A));groups=[[j for j,d in enumerate(deg) if d==deg[i]] for i in range(8)];assignment={};out=[]
 def rec():
  if len(assignment)==8:out.append(tuple(assignment[i] for i in range(8)));return
  i=next(i for i in range(8) if i not in assignment)
  for j in groups[i]:
   if j in assignment.values() or any(A[i][u]!=A[j][v] for u,v in assignment.items()):continue
   assignment[i]=j;rec();del assignment[i]
 rec();return out
opts={r['source']:r for r in P['source_options']};need(len(opts)==28,'exact28 distinct owners');labels={};allresources={};domain_reports=[]
for s,r in opts.items():
 A=matrix(int(nodes[s]['canon_a']));need(A==r['canonical_matrix'],'actual source matrix');need(autos(A)==[tuple(range(8))],'actual trivial source kernel / automorphism group');need(nodes[s]['class_a']!=nodes[s]['class_b'],'two distinct ordinary converse source classes')
 native_paths=[]
 def enumerate_paths(pp,left):
  if not left:native_paths.append(pp);return
  for v in sorted(left):
   if not pp or A[pp[-1]][v]:enumerate_paths(pp+[v],left-{v})
 enumerate_paths([],set(range(8)));need(2*len(native_paths)==int(nodes[s]['fibre']),'actual full paired native path inventory')
 target_capacities=Counter(merged[owner[word(A,pp)^full]] for pp in native_paths)
 f=r['local_to_physical'];need(sorted(f)==list(range(8)),'literal owner bijection');C=[[0]*8 for _ in range(8)]
 for i in range(8):
  for j in range(8):C[f[i]][f[j]]=A[i][j]
 need(C==r['physical_matrix'],'actual physical source frame');need(r['physical_owner_symbols']==[{0:'shared:first',7:'shared:last',3:'shared:arc0',4:'shared:arc1'}.get(i,f'private:{s}:{i}') for i in range(8)],'only declared shared owners identified; all remaining owners source-private')
 before,after=r['options'];p=before['path'];q=after['path'];fixed=[[j,u] for j,u in enumerate(p) if q[j]==u];free=[j for j in range(8) if p[j]!=q[j]]
 need(r['fixed_local_clock_owners']==fixed and len(free)==r['cost'],'all common pins, exact owner-clock cost');need(r['fixed_physical_clock_owners']==[[j,f[u]] for j,u in fixed],'physical fixed clocks')
 first,last,k,a,b=r['shared_pins'];need(len({first,last,a,b})==4 and 1<=k<=5,'four actual shared pin roles');need(f[first]==0 and f[last]==7 and f[a]==3 and f[b]==4,'shared owner alignment')
 domains=[]
 for values in permutations([p[j] for j in free]):
  pp=p[:]
  for j,v in zip(free,values):pp[j]=v
  if all(A[u][v] for u,v in zip(pp,pp[1:])):
   w=word(A,pp);t=merged[owner[w^full]];need(owner[w]==int(nodes[s]['class_a']),'pinned path ordinary source owner');rr=orbit(w);need(len(rr)==4,'complete pinned resource pure black');res=min(rr);need(res not in allresources,'global no reuse, including extra pinned domain');allresources[res]=(s,pp)
   domains.append(dict(path=pp,word=w,target=t,resource=res))
 need(sorted(x['path'] for x in domains if x['target'] in [U,V])==sorted([p,q]),'exact pair-target pinned domain is binary')
 for rec in r['options']:
  pp=rec['path'];physical=[f[u] for u in pp];need(physical==rec['physical_path'],'exact physical path');need(sorted(pp)==list(range(8)) and all(A[u][v] for u,v in zip(pp,pp[1:])),'literal Hamiltonian source path')
  need(all(physical[j]==v for j,v in r['shared_physical_clock_owners']),'all shared indexed clocks');need(all(physical[j]==v for j,v in r['fixed_physical_clock_owners']),'all private fixed clocks')
  need(C[3][4]==1,'shared indexed arc');w=word(C,physical);need(w==rec['mask'] and merged[owner[w^full]]==rec['target'],'actual forward resource and target');need(colored[tuple(sorted((s,rec['target'])))]==(0,2*target_capacities[rec['target']]),'actual pure-black paired multiplicity equals distinct native inventory');need(owner[reflect(w)]==int(nodes[s]['class_b']) and merged[owner[reflect(w)^full]]==rec['target'],'literal converse source/target owner and actual paired factor2')
  B=partner(C,physical);rev=physical[::-1];closures=[('I',C,physical),('F',B,physical),('R',transpose(C),rev),('FR',transpose(B),rev)]
  need(partner(B,physical)==C and partner(transpose(C),rev)==transpose(B),'literal complement involution and converse square')
  observed=[]
  for stored,(action,T,ppp) in zip(rec['literal_closure'],closures):
   need(stored==dict(action=action,matrix=T,path=ppp),'retained exact matrix/path kernel lift');need(all(T[u][v] for u,v in zip(ppp,ppp[1:])),'closure Hamiltonian arcs');observed.append(word(T,ppp))
  need(observed==[w,w^full,reflect(w),reflect(w)^full],'literal F/R closure word action')
  need((X[s]-X[rec['target']])%8==(X[s]-X[U])%8,'chosen ordinary per-slot X mod8 phase')
 labels[s]=(sum((Y[j][s]%2)<<j for j in range(4)),X[s]%128);need(r['source_label']==[*(Y[j][s]%2 for j in range(4)),X[s]%128],'native charge label')
 domain_reports.append(dict(source=s,cost=r['cost'],all_pinned_paths=domains,fixed_local_clock_owners=fixed,automorphisms=1))
need(len(set(labels.values()))==28,'all source labels distinct');need(len(allresources)==57,'full common-pin domain resources, including1559 extra');need(Counter(x['target'] for r in domain_reports for x in r['all_pinned_paths'] if x['target'] not in [U,V])=={1559:1},'only extra target has zero initial margin and cannot be selected')

W=json.loads((H/'image_witnesses.json').read_text());need(W['modulus']==256 and W['input_sha256']==sha(blobs['merged_edges.tsv.gz']),'original complete Laplacian image witness pin');charge_list=json.loads((H/'charges.json').read_text());need(len(W['solutions'])==len(charge_list),'exact certificate count')
if args.include_eighteen:
 extra=json.loads((H/'image_witnesses_eighteen.json').read_text());need(extra['modulus']==256 and extra['input_sha256']==W['input_sha256'],'extension image antecedent');W['solutions'].extend(extra['solutions']);charge_list.extend(json.loads((H/'charges_eighteen.json').read_text()));need(len(W['solutions'])==len(charge_list),'extended certificate count')
if args.include_eighteen37:
 extra=json.loads((H/'image_witnesses_eighteen37.json').read_text());need(extra['modulus']==256 and extra['input_sha256']==W['input_sha256'],'cost37 extension image antecedent');W['solutions'].extend(extra['solutions']);charge_list.extend(json.loads((H/'charges_eighteen37.json').read_text()));need(len(W['solutions'])==len(charge_list),'cost37 extended certificate count')
if args.include_variants:
 extra=json.loads((H/'image_witnesses_variants.json').read_text());need(extra['modulus']==256 and extra['input_sha256']==W['input_sha256'],'variant image antecedent');W['solutions'].extend(extra['solutions']);charge_list.extend(json.loads((H/'charges_variants.json').read_text()));need(len(W['solutions'])==len(charge_list),'variant certificate count')
images={};image_reports=[]
for rec,wanted in zip(W['solutions'],charge_list):
 q={int(k):v for k,v in wanted['charge'].items()};need(dict(rec['charge'])==q,'literal requested charge');z=rec['potential'];need(len(z)==N and z[-1]==0 and all(type(x)is int and 0<=x<256 for x in z),'complete grounded integral vector');a=lap(z)
 correction=[[i,(v-q.get(i,0))//256] for i,v in enumerate(a) if v!=q.get(i,0)];need(all((v-q.get(i,0))%256==0 for i,v in enumerate(a)),'full original Lz=charge modulo256');need(correction==rec['integer_correction_nonzero'],'full exact integer correction ledger')
 images[wanted['name']]=q;image_reports.append(dict(name=wanted['name'],charge=rec['charge'],potential_sha256=sha(json.dumps(z,separators=(',',':')).encode())))
need(images['common_target']=={U:2,V:-2},'common target charge image')

def signature(mask,ss):
 n=0;p=0;x=0
 for i,s in enumerate(ss):
  if mask>>i&1:n+=1;p^=labels[s][0];x+=labels[s][1]
 return n,p,x%128

def accepted(ss,start):
 target=signature(start,ss);out=[]
 # Exact finite subset DP, not a search certificate copied from generation.
 states=[(0,0,0)]
 for s in ss:states += [(n+1,p^labels[s][0],(x+labels[s][1])%128) for n,p,x in states]
 for m,t in enumerate(states):
  if t==target:out.append(m)
 return out

def route(states,start,end,budget,weights=None,slot_cap=None):
 prev={start:None};q=deque([start])
 while q:
  a=q.popleft()
  if a==end:
   path=[]
   while a is not None:path.append(a);a=prev[a]
   return path[::-1],sorted(prev)
  for b in states:
   if b not in prev and (slot_cap is None or (a^b).bit_count()<=slot_cap) and ((a^b).bit_count() if weights is None else sum(w for i,w in enumerate(weights) if (a^b)>>i&1))<=budget:prev[b]=a;q.append(b)
 return None,sorted(prev)

def minimum_work(states,start,end,weights,slot_cap,clock_cap=None):
 dist={start:(0,0)};prev={start:None};queue=[(0,0,start)]
 while queue:
  cost,events,a=heapq.heappop(queue)
  if (cost,events)!=dist[a]:continue
  for b in states:
   if a==b or (a^b).bit_count()>slot_cap:continue
   work=sum(w for i,w in enumerate(weights) if (a^b)>>i&1)
   if clock_cap is not None and work>clock_cap:continue
   candidate=(cost+work,events+1)
   if candidate<dist.get(b,(10**9,10**9)):
    dist[b]=candidate;prev[b]=a;heapq.heappush(queue,(*candidate,b))
 need(end in dist,'weighted endpoint reachable');path=[];a=end
 while a is not None:path.append(a);a=prev[a]
 return dict(total_owner_clock_work=dist[end][0],events=dist[end][1],path=path[::-1],distance_labels={str(k):list(v) for k,v in sorted(dist.items())})

reports=[]
for r in P['routes']:
 ss=r['sources'];path=r['path'];n=len(ss);start=path[0];end=path[-1];base=r.get('base_source_order',P['base_source_order']);m=len(base);k=m//2;initial_base=((1<<k)-1)<<k;final_base=(1<<k)-1;need(ss==base+r['buffers'],'exact route source order');need(start==initial_base+(r['initial']<<m) and end==final_base+(r['initial']<<m),'original endpoint and every buffer restored');need(accepted(base,initial_base)==[final_base,initial_base],'underlying original circuit has exactly two states');legal=accepted(ss,start);need(all(m in legal for m in path),'every visited state satisfies complete Q test');need(len(legal)==r['states'],'complete exact guarded fiber size')
 minimum=None;separation={}
 for budget in range(0,m+1,2):
  found,component=route(legal,start,end,budget)
  if found:minimum=budget;break
  separation[str(budget)]=component
 need(minimum==r['budget'],'exact bottleneck budget')
 shortest_path,_=route(legal,start,end,minimum);weights=[opts[s]['cost'] for s in ss];work_optimum=minimum_work(legal,start,end,weights,minimum)
 for clock_cap in range(sum(weights)+1):
  clock_path,_=route(legal,start,end,clock_cap,weights,minimum)
  if clock_path:break
 clock_work_optimum=minimum_work(legal,start,end,weights,minimum,clock_cap)
 accounting=dict(slot_budget=minimum,minimum_events=len(shortest_path)-1,event_minimum_witness=shortest_path,minimum_work=work_optimum,minimum_clock_bottleneck=clock_cap,minimum_events_at_minimum_clock_bottleneck=len(clock_path)-1,minimum_clock_bottleneck_event_witness=clock_path,minimum_work_at_minimum_clock_bottleneck=clock_work_optimum,minimum_work_also_minimizes_events=work_optimum['events']==len(shortest_path)-1)
 step_reports=[];resources_used=set()
 for state in path:
  resources=[min(orbit(opts[s]['options'][state>>i&1]['mask'])) for i,s in enumerate(ss)];need(len(set(resources))==n,'no resource reused within state');resources_used.update(resources)
  need(Counter(opts[s]['options'][state>>i&1]['target'] for i,s in enumerate(ss))==Counter(opts[s]['options'][start>>i&1]['target'] for i,s in enumerate(ss)),'every physical target margin, including outside targets')
 for j,(a,b) in enumerate(zip(path,path[1:])):
  delta={s:(b>>i&1)-(a>>i&1) for i,s in enumerate(ss) if ((a^b)>>i)&1};need(sum(delta.values())==0,'zero source balance');charge={s:2*v for s,v in delta.items()};need(images[r.get('step_image_names',[f"{r['name']}_step_{t+1}" for t in range(len(path)-1)])[j]]==charge,'step has explicit full image lift')
  signed_edges=defaultdict(int)
  for s,v in delta.items():signed_edges[tuple(sorted((s,U)))]+=-2*v;signed_edges[tuple(sorted((s,V)))]+=2*v
  matrix_delta=defaultdict(int)
  for (s,t),w in signed_edges.items():matrix_delta[s,s]+=w;matrix_delta[t,t]+=w;matrix_delta[s,t]-=w;matrix_delta[t,s]-=w
  support=sorted(set(delta)|{U,V});v={U:1,V:-1}
  for i in support:
   need(sum(matrix_delta[i,k] for k in support)==0,'physical row sum zero')
   for k in support:need(matrix_delta[i,k]==2*(delta.get(i,0)*v.get(k,0)+v.get(i,0)*delta.get(k,0)),'complete physical-row signed Laplacian identity')
  currents=[]
  for chi in chars:
   y={i:sum(matrix_delta[i,k]*chi[k] for k in support) for i in support};need(all(t%256==0 for t in y.values()),'every physical row current preserved for authentic generator');currents.append([[i,t] for i,t in y.items() if t])
  clocks=sum(opts[s]['cost'] for s in delta);step_reports.append(dict(step=j+1,old_mask=a,new_mask=b,plus=sorted(s for s,d in delta.items() if d==1),minus=sorted(s for s,d in delta.items() if d==-1),changed_slots=len(delta),changed_owner_clocks=clocks,signed_edge_change=[[s,t,w] for (s,t),w in sorted(signed_edges.items())],integer_generator_currents=currents))
 need(max(s['changed_slots'] for s in step_reports)==minimum,'saved route achieves optimum');need(len(resources_used)==2*n,'every selected resource alternative available and disjoint over union')
 reentries=[]
 for i,s in enumerate(ss):
  seen={start>>i&1};previous=start>>i&1
  for step,state in enumerate(path[1:],1):
   current=state>>i&1
   if current!=previous:
    if current in seen:reentries.append(dict(step=step,source=s,selected_bit=current,resource=min(orbit(opts[s]['options'][current]['mask']))))
    seen.add(current)
   previous=current
 need(len(reentries)==sum(x['changed_slots'] for x in step_reports)-n,'exact released-resource reentry ledger');need(len(reentries)>0,'positive compensated route is NOT strict temporal nonreturn')
 reports.append(dict(name=r['name'],sources=ss,buffers=r['buffers'],initial=r['initial'],path=path,optimal_accounting=accounting,released_resource_reentries=reentries,satisfies_strict_temporal_nonreturn=False,legal_states=legal,minimum_maximum_changed_slots=minimum,below_threshold_components=separation,steps=step_reports,normalized_resources_per_state=n,normalized_resources_union=len(resources_used),closed_ports_per_state=4*n,closed_ports_union=4*len(resources_used),max_owner_clock_changes=max(s['changed_owner_clocks'] for s in step_reports),total_owner_clock_changes=sum(s['changed_owner_clocks'] for s in step_reports)))

src=sorted(opts);pair_buckets=defaultdict(list)
for a,b in combinations(src,2):pair_buckets[(labels[a][0]^labels[b][0],(labels[a][1]+labels[b][1])%128)].append((a,b))
four=[]
for pairs in pair_buckets.values():
 for a,b in combinations(pairs,2):
  need(set(a).isdisjoint(b),'distinct labels force disjoint pair collision');need(sum(s in [1638,1974] for s in a)==sum(s in [1638,1974] for s in b),'integer invariant b1638+b1974 for every4move');four.append([list(a),list(b)])
need(len(four)==27,'complete28 four-move list');need(all(s in P['base_source_order'][:8] for s in [1638,1974]),'endpoint changes invariant fromzero totwo')
# Independent exhaustive one-cost2-buffer bottlenecks: two are necessary for budget6 in20.
base=P['base_source_order'];outside_cost2=sorted(s for s in opts if s not in base and opts[s]['cost']==2);need(outside_cost2==[1394,1921,2905,3369],'complete cost2 outside reservoir')
cost2_reports=[]
for s in outside_cost2:
 for initial in [0,1]:
  ss=base+[s];start=65280+(initial<<16);end=255+(initial<<16);legal=accepted(ss,start);none,component=route(legal,start,end,6);found,_=route(legal,start,end,8);need(none is None and found is not None,'every onecost2 buffer has exact minimum8');cost2_reports.append(dict(source=s,initial=initial,legal_states=len(legal),six_move_start_component=component,eight_move_witness=found))
base_legal=accepted(base,65280);need(base_legal==[255,65280],'zero buffers retain primitive16slot circuit')
clock_reports=[]
for s in sorted(set(opts)-set(base)):
 for initial in [0,1]:
  ss=base+[s];start=65280+(initial<<16);end=255+(initial<<16);legal=accepted(ss,start);weights=[opts[u]['cost'] for u in ss]
  for budget in range(33):
   found,component=route(legal,start,end,budget,weights)
   if found:break
  need(found is not None,'weighted route exists for every single buffer')
  clock_reports.append(dict(source=s,initial=initial,minimum_maximum_owner_clock_changes=budget,witness=found))
need(min(r['minimum_maximum_owner_clock_changes'] for r in clock_reports)==13,'one-buffer sharp global clock budget13')
out=dict(status='PASS',checks=checks,native_domain_reports=domain_reports,normalized_resources_full28_common_pin_union=len(allresources),literal_closed_ports_full28_common_pin_union=4*len(allresources),image_reports=image_reports,routes=reports,all_four_slot_trades=four,integer_obstruction=dict(weights={1638:1,1974:1},endpoint_values=[0,2],applies_to='All28 supplied slots, every buffer initialization, any path length; fixed actual pins and all target margins.'),all_one_cost2_buffer_tests=cost2_reports,all_one_buffer_clock_tests=clock_reports,scope='Native partial observer packets with available distinct resources and compatible unchanged remainder. Every physical row and genuine dyadic bank is preserved, plus selected ordinary X mod8. No complete-atlas rewire, common-host amalgam, universal ordinary-port phase preservation, or all-size supply theorem.')
text=json.dumps(out,indent=2,sort_keys=True)+'\n'
if args.out:args.out.write_text(text)
else:print(text,end='')
