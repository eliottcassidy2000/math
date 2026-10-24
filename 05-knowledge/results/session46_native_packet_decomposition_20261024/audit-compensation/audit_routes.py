#!/usr/bin/env python3
"""Independent route audit using native records and literal nodewise currents."""
import argparse,csv,gzip,hashlib,itertools,json,struct
from collections import Counter,deque
from pathlib import Path
P=argparse.ArgumentParser();P.add_argument('--repo',type=Path,required=True);P.add_argument('--packet',type=Path,required=True);P.add_argument('--foundation',type=Path,required=True);P.add_argument('--out',type=Path,required=True);P.add_argument('--images',type=Path);P.add_argument('--extra-packet',type=Path);P.add_argument('--extra-images',type=Path);P.add_argument('--extra37-packet',type=Path);P.add_argument('--extra37-images',type=Path);a=P.parse_args()
checks=0
def require(ok,why):
 global checks
 checks+=1
 if not ok:raise RuntimeError(why)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
R=a.repo/'05-knowledge/results';G=R/'metagraph_order8_20261023'
paths={'nodes':G/'merged_nodes.tsv','owners':G/'owner_map.u16.gz','edges':G/'merged_edges.tsv.gz','electrical':R/'metagraph_order8_electrical_20261023/receipt.json','packet':a.packet,'foundation':a.foundation}
if a.images:paths['images']=a.images
if a.extra_packet:paths['extra_packet']=a.extra_packet
if a.extra_images:paths['extra_images']=a.extra_images
if a.extra37_packet:paths['extra37_packet']=a.extra37_packet
if a.extra37_images:paths['extra37_images']=a.extra37_images
hashes={k:digest(v)for k,v in paths.items()}
F=json.loads(a.foundation.read_text());packet=json.loads(a.packet.read_text())
require(F['status']=='PASS'and all(F['input_sha256'][k]==hashes[k]for k in('nodes','owners','edges','electrical')),'independent foundation tied to same literal originals')
nodes=[{k:int(v)for k,v in r.items()}for r in csv.DictReader(paths['nodes'].open(),delimiter='\t')];N=len(nodes)
owners=struct.unpack('<'+'H'*(1<<21),gzip.decompress(paths['owners'].read_bytes()));merge={r[k]:r['merged_id']for r in nodes for k in('class_a','class_b')}
edges=[tuple(int(r[k])for k in('source','target','blue','black'))for r in csv.DictReader(gzip.open(paths['edges'],'rt'),delimiter='\t')]
E=json.loads(paths['electrical'].read_text())['order8'];Y=E['radical_lifts_mod4'];X=E['higher_layers'][-1]['lifts'][0];chars=[[64*y for y in yy]for yy in Y[:4]]+[X]
U,V=1965,2234;full=(1<<21)-1;pairs=list(itertools.combinations(range(8),2));chords=[(i,j)for i,j in pairs if j>i+1];idx={e:k for k,e in enumerate(chords)}
def word(B,p):return sum(B[p[i]][p[j]]<<k for k,(i,j)in enumerate(chords))
def reflect(w):return sum(((w>>idx[7-j,7-i])&1)<<k for k,(i,j)in enumerate(chords))
def resource(w):return {w,w^full,reflect(w),reflect(w)^full}
def transpose(B):return [list(r)for r in zip(*B)]
def partner(B,p):
 keep={tuple(sorted(z))for z in zip(p,p[1:])};C=[[0]*8 for _ in range(8)]
 for i,j in pairs:C[i][j]=B[i][j]if(i,j)in keep else B[j][i];C[j][i]=1-C[i][j]
 return C
def lap(x):
 out=[0]*N
 for i,j,b,k in edges:z=(b+k)*(x[i]-x[j]);out[i]+=z;out[j]-=z
 return out
for chi in chars:require(all(z%256==0 for z in lap(chi)),'each separator harmonic on full original Laplacian')
foundation={r['source']:r for r in F['sources']};slots={r['source']:r for r in packet['source_options']};require(set(slots)==set(foundation),'packet retains exact entire authenticated supply')
resources={};slot_report=[]
for s,row in slots.items():
 fr=foundation[s];B=fr['matrix'];f=row['local_to_physical'];M=row['physical_matrix']
 require(row['canonical_matrix']==B,'actual canonical matrix')
 expected_symbols=[{0:'shared:first',7:'shared:last',3:'shared:arc0',4:'shared:arc1'}.get(i,f'private:{s}:{i}')for i in range(8)]
 require(row['physical_owner_symbols']==expected_symbols and len(set(expected_symbols))==8,'only four declared owner roles shared, all other physical symbols source-private')
 require(sorted(f)==list(range(8))and all(M[f[i]][f[j]]==B[i][j]for i in range(8)for j in range(8)),'physical matrix actual relabeling, not a substituted host')
 first,last,k,c,d=row['shared_pins'];require([f[z]for z in(first,last,c,d)]==[0,7,3,4],'literal common physical endpoint and arc owners')
 require(row['fixed_local_clock_owners']==list(map(list,zip(fr['fixed_clocks'],fr['fixed_owner_values']))),'all fixed local clocks exactly recovered')
 require(row['fixed_physical_clock_owners']==[[j,f[u]]for j,u in row['fixed_local_clock_owners']],'all fixed physical clocks exactly recovered')
 require(row['shared_physical_clock_owners']==[[0,0],[7,7],[k,3],[k+1,4]],'complete declared common pins')
 require(len(row['options'])==2,'two explicit in-margin alternatives')
 for bit,opt in enumerate(row['options']):
  t=(U,V)[bit];require(opt['target']==t,'bit convention independently checked')
  native=next(z for z in fr['complete_pinned_domain']if z['target']==t)
  require(opt['path']==native['path']and opt['mask']==native['word'],'actual option recovered by independent all-permutation census')
  pp=[f[v]for v in opt['path']];require(opt['physical_path']==pp,'exact retained owner lift')
  require(all(pp[j]==u for j,u in row['fixed_physical_clock_owners']),'all physical pins preserved by alternative')
  K=partner(M,pp);rev=pp[::-1];expected=[('I',M,pp),('F',K,pp),('R',transpose(M),rev),('FR',transpose(K),rev)]
  require(len(opt['literal_closure'])==4,'full literal closure schema')
  actualwords=[]
  for declared,(action,C,path)in zip(opt['literal_closure'],expected):
   require(declared=={'action':action,'matrix':C,'path':path},'literal closure matrices and paths verified, not just word labels')
   require(sorted(path)==list(range(8))and all(C[i][j]for i,j in zip(path,path[1:])),'actual closure Hamiltonian path')
   actualwords.append(word(C,path))
  w=opt['mask'];require(actualwords==[w,w^full,reflect(w),reflect(w)^full],'literal four-resource words')
  require(merge[owners[w]]==s and merge[owners[w^full]]==t,'actual original owner and partner target')
  require(len(set(actualwords))==4 and set(actualwords)==set(native['closed_ports']),'actual black normalized resource')
  resources[s,bit]=set(actualwords)
 require(sum(x!=y for x,y in zip(*[z['physical_path']for z in row['options']]))==row['cost'],'actual ordinary owner-clock change cost')
 slot_report.append({'source':s,'resource_ids':[min(resources[s,i])for i in range(2)],'cost':row['cost']})
require(len(set.union(*resources.values()))==4*len(resources)==224,'entire selected two-target supply is globally injective')

base=packet['base_source_order'];require(len(base)==len(set(base))==16,'original sixteen distinct sources')
no_buffer={'name':'no_buffer','sources':base,'buffers':[],'initial':0,'path':[65280,255],'states':2,'budget':16,'step_sizes':[16]}
routes=[no_buffer]+packet['routes']
for extra_path in(a.extra_packet,a.extra37_packet):
 if not extra_path:continue
 extra=json.loads(extra_path.read_text())
 for r in extra['routes']:
  b=r['base_source_order'];k=len(b)//2;last=(1<<k)-1;first=last<<k
  routes.append({'name':'no_buffer_'+str(len(b))+'_cost'+str(sum(slots[s]['cost']for s in b)),'base_source_order':b,'sources':b,'buffers':[],'initial':0,'path':[first,last],'states':2,'budget':len(b),'step_sizes':[len(b)]})
  routes.append(r)
route_reports=[];required_charges={tuple([(U,2),(V,-2)])}

def state_current(sources,mask,chi):
 # Every literal normalized black tile contributes two actual undirected black
 # edges. This assembles ALL N physical metagraph rows, including zero rows.
 out=[0]*N
 for j,s in enumerate(sources):
  t=slots[s]['options'][(mask>>j)&1]['target'];v=2*(chi[s]-chi[t]);out[s]+=v;out[t]-=v
 return out

def bfs(states,start,end,budget,costs=None):
 seen={start:None};q=deque([start])
 while q:
  x=q.popleft()
  for y in states:
   weight=(x^y).bit_count()if costs is None else sum(c for j,c in enumerate(costs)if(x^y)>>j&1)
   if y not in seen and weight<=budget:seen[y]=x;q.append(y)
 if end not in seen:return len(seen),None
 path=[end]
 while path[-1]!=start:path.append(seen[path[-1]])
 return len(seen),path[::-1]

for route in routes:
 sources=route['sources'];n=len(sources);start,end=route['path'][0],route['path'][-1]
 routebase=route.get('base_source_order',base);b=len(routebase);k=b//2;last=(1<<k)-1;first=last<<k
 require(b%2==0 and sources==routebase+route['buffers']and len(sources)==len(set(sources)),'base packet plus genuinely distinct compensator sources')
 require(start==first+(route['initial']<<b)and end==last+(route['initial']<<b),'exact original endpoints and literal restored buffer bits')
 require(not set(sources).intersection((U,V)),'source and target roles disjoint')
 critical=[sum((mask>>j)&1 for j,s in enumerate(sources)if s in(1638,1974))for mask in(start,end)]
 require(critical[0]!=critical[1],'full28 four-move invariant separates these endpoints, even with restored extra buffers')
 baseline=[state_current(sources,start,chi)for chi in chars]
 # Enumerate the complete binary domain by actual nodewise current contributions,
 # not by the producer's source-code subset-sum helper. Target margins fix popcount.
 source_unchanged=all(2*(chi[U]-chi[V])%256==0 for chi in chars)
 require(source_unchanged,'every possible source row has identical residues for both alternatives')
 contributions=[tuple((-2*(chi[s]-chi[V]))%256 for chi in chars)for s in sources]
 target=tuple(row[V]%256 for row in baseline);sums=[(0,)*5]*(1<<n);accepted=[];balanced=0;separators=Counter()
 for mask in range(1<<n):
  if mask:
   low=mask&-mask;prev=mask^low;j=low.bit_length()-1;sums[mask]=tuple((x+y)%256 for x,y in zip(sums[prev],contributions[j]))
  if mask.bit_count()!=start.bit_count():continue
  balanced+=1
  # Source-row invariance plus zero column sum guarantees the U row once V agrees.
  if sums[mask]==target:accepted.append(mask)
  else:
   k=next(k for k in range(5)if sums[mask][k]!=target[k]);residue=(sums[mask][k]-target[k])%256
   require(residue!=0,'every rejected margin-correct state has an actual harmonic separator at physical row V')
   separators[str(k)]+=1
 require(len(accepted)==route['states'],'exact complete accepted-fibre size')
 state_reports=[]
 for mask in accepted:
  choices=[slots[s]['options'][mask>>j&1]for j,s in enumerate(sources)]
  require(Counter(z['target']for z in choices)==Counter(slots[s]['options'][start>>j&1]['target']for j,s in enumerate(sources)),'all accepted literal target margins')
  T=[resources[s,mask>>j&1]for j,s in enumerate(sources)]
  require(len(set.union(*T))==4*n and len({min(t)for t in T})==n,'every accepted state normalized-resource and closed-port injection')
  currents=[state_current(sources,mask,chi)for chi in chars]
  delta=[[value-old for value,old in zip(row,baseline[k])]for k,row in enumerate(currents)]
  require(all(all(v%256==0 for v in row)for row in delta),'direct full physical row current predicate on every accepted state')
  require(all(all(z['physical_path'][clock]==owner for clock,owner in slots[s]['fixed_physical_clock_owners'])for s,z in zip(sources,choices)),'every accepted state all physical-clock pins')
  state_reports.append({'mask':mask,'target_counts':dict(Counter(z['target']for z in choices)),'normalized_resources':[min(t)for t in T],'integer_current_delta_nonzero':[[[i,v]for i,v in enumerate(row)if v]for row in delta]})
 require(all(mask in accepted for mask in route['path']),'every advertised intermediate literal state accepted')
 lower=[];budget=None;optimal=None
 for b in range(2,n+1,2):
  count,path=bfs(accepted,start,end,b)
  if path is not None:budget=b;optimal=path;break
  lower.append({'budget':b,'reachable_from_start':count})
 require(budget==route['budget'],'exact selected complete-pinned-menu minimax slot budget')
 clock_minimum=None;clock_path=None
 for clock_budget in range(1,1+sum(slots[s]['cost']for s in sources)):
  _,p=bfs(accepted,start,end,clock_budget,[slots[s]['cost']for s in sources])
  if p is not None:clock_minimum=clock_budget;clock_path=p;break
 step_reports=[];ever_occupied={min(resources[s,start>>j&1])for j,s in enumerate(sources)};reentry_events=[]
 for old,new in zip(route['path'],route['path'][1:]):
  old_set={min(resources[s,old>>j&1])for j,s in enumerate(sources)};new_set={min(resources[s,new>>j&1])for j,s in enumerate(sources)}
  reentries=sorted((new_set-old_set)&ever_occupied);reentry_events.append({'from':old,'to':new,'reentered_normalized_resources':reentries});ever_occupied.update(new_set)
  changed=[j for j in range(n)if(old^new)>>j&1];sign={sources[j]:((new>>j)&1)-((old>>j)&1)for j in changed}
  require(sum(sign.values())==0,'each route step balanced')
  q=tuple(sorted((s,2*v)for s,v in sign.items()));required_charges.add(q)
  delta_edges=Counter()
  for j in changed:
   s=sources[j];old_t=slots[s]['options'][old>>j&1]['target'];new_t=slots[s]['options'][new>>j&1]['target']
   delta_edges[tuple(sorted((s,old_t)))]-=2;delta_edges[tuple(sorted((s,new_t)))]+=2
  for k,chi in enumerate(chars):
   literal=[0]*N
   for(i,j),w in delta_edges.items():v=w*(chi[i]-chi[j]);literal[i]+=v;literal[j]-=v
   difference=[x-y for x,y in zip(state_current(sources,new,chi),state_current(sources,old,chi))]
   require(literal==difference and all(v%256==0 for v in literal),'each step independently assembled full physical current difference')
  owner_cost=sum(sum(x!=y for x,y in zip(slots[sources[j]]['options'][old>>j&1]['physical_path'],slots[sources[j]]['options'][new>>j&1]['physical_path']))for j in changed)
  require(owner_cost==sum(slots[sources[j]]['cost']for j in changed),'literal owner-clock support agrees with actual changed slots')
  step_reports.append({'from':old,'to':new,'changed_sources':[sources[j]for j in changed],'slot_support':len(changed),'ordinary_owner_clock_changes':owner_cost,'integer_X_target_debt_at_U':2*sum(sign[s]*X[s]for s in sign),'source_charge':[list(z)for z in q]})
 require([r['slot_support']for r in step_reports]==route['step_sizes']and max(r['slot_support']for r in step_reports)==budget,'route attains exact slot bottleneck')
 visited=[resources[s,mask>>j&1]for mask in route['path']for j,s in enumerate(sources)]
 route_reports.append({'name':route['name'],'sources':sources,'buffers':route['buffers'],'critical_pair_invariant_endpoint_values':critical,'all_binary_states':1<<n,'target_margin_states':balanced,'accepted_states':accepted,'accepted_state_records':state_reports,'rejected_separators_by_generator':dict(separators),'exact_minimax_slot_budget':budget,'exact_minimax_ordinary_clock_budget':clock_minimum,'clock_optimal_path':clock_path,'smaller_budget_reachability':lower,'independent_optimal_path':optimal,'advertised_path':route['path'],'steps':step_reports,'closed_ports_per_state':4*n,'closed_ports_visited_union':len(set.union(*visited)),'ordinary_clock_bottleneck_on_advertised_route':max(r['ordinary_owner_clock_changes']for r in step_reports),'resource_reentry_ledger':reentry_events,'resource_reentry_count':sum(len(r['reentered_normalized_resources'])for r in reentry_events),'satisfies_strict_temporal_nonreturn':not any(r['reentered_normalized_resources']for r in reentry_events)})

image_report=None
if a.images:
 images=json.loads(a.images.read_text());require(images['modulus']==256 and images['input_sha256']==hashes['edges'],'full image witnesses tied to original graph and modulus')
 for extra_image_path in(a.extra_images,a.extra37_images):
  if not extra_image_path:continue
  additional=json.loads(extra_image_path.read_text());require(additional['modulus']==256 and additional['input_sha256']==hashes['edges'],'extended image witnesses tied to original graph and modulus');images['solutions'].extend(additional['solutions'])
 verified=set()
 for solution in images['solutions']:
  z=solution['potential'];charge=dict(solution['charge']);correction=dict(solution['integer_correction_nonzero'])
  require(len(z)==N and z[-1]==0 and all(type(v)is int and 0<=v<256 for v in z),'full grounded image witness vector')
  require(sum(charge.values())==sum(correction.values())==0,'source charge and correction zero sums')
  actual=lap(z)
  for i,v in enumerate(actual):require(v==charge.get(i,0)+256*correction.get(i,0),'every original physical row exact integer image equation')
  verified.add(tuple(sorted(charge.items())))
 # Every unbuffered endpoint charge is a sum of verified compensated step charges.
 base_charges={tuple(r['sources']):tuple(tuple(x)for x in r['steps'][0]['source_charge'])for r in route_reports if r['name'].startswith('no_buffer')}
 required_charges.difference_update(base_charges.values())
 require(required_charges<=verified,'common target and every compensated step have verified image antecedents')
 for r in route_reports:
  if r['name'].startswith('no_buffer'):continue
  original=tuple(r['sources'][:-len(r['buffers'])]);combined=Counter()
  for step in r['steps']:combined.update(dict(step['source_charge']))
  require({s:v for s,v in combined.items()if v}==dict(base_charges[original]),'verified compensated step witnesses add to original endpoint image charge')
 image_report={'verified_vectors':len(images['solutions']),'verified_distinct_charges':len(verified),'required_compensated_and_common_charges':len(required_charges),'physical_rows_per_equation':N}

report={'status':'PASS','checks':checks,'input_sha256':hashes,'scope':'All complete common-clock U/V domains of selected route source sets; full physical nodewise current inclusion. Global all-dyadic conclusion uses independently rerun Session45 filtration. Exact minima belong to each specified menu; full28 critical-pair invariant establishes the larger-universe six lower bound. No full-atlas rewiring or every-ordinary-character claim.','sources':slot_report,'routes':route_reports,'exact_image_equations':image_report}
a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','checks':checks,'routes':[{k:r[k]for k in('name','all_binary_states','accepted_states','exact_minimax_slot_budget','ordinary_clock_bottleneck_on_advertised_route')}for r in route_reports],'images':image_report},sort_keys=True))
