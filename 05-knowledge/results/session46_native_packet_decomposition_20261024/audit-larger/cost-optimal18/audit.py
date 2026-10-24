#!/usr/bin/env python3
"""Independent literal audit of Session46 native circuits; no producer-code imports.
Failures are explicit and remain enabled under python -O.
"""
import argparse,csv,gzip,hashlib,itertools,json,math,struct
from collections import Counter,deque
from pathlib import Path
P=argparse.ArgumentParser();P.add_argument('--repo',type=Path,required=True);P.add_argument('--packet',type=Path,required=True);P.add_argument('--out',type=Path,required=True);P.add_argument('--separators',type=Path);a=P.parse_args()
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
R=a.repo/'05-knowledge/results';B=R/'metagraph_order8_20261023';D=a.packet
files={'nodes':B/'merged_nodes.tsv','owners':B/'owner_map.u16.gz','edges':B/'merged_edges.tsv.gz','electrical':R/'metagraph_order8_electrical_20261023/receipt.json','pool':D/'pool.json','inventory':D/'source_inventory.json','candidates':D/'candidates.json','images':D/'image_witnesses.json','pins':D/'SOURCE_PINS.json','dyadic_audit':R/'session45_injective_multisource_moves_20261023/audit-dyadic/normal.json'}
hashes={k:digest(v)for k,v in files.items()}
pins=json.loads(files['pins'].read_text())
for file,sha in pins['originals'].items():need(digest(a.repo/file)==sha,'original pinned hash: '+file)
need(hashes['inventory']==pins['inherited_inventory_sha256'],'inherited full inventory pin')
N=[{k:int(v)for k,v in row.items()}for row in csv.DictReader(files['nodes'].open(),delimiter='\t')]
need([r['merged_id']for r in N]==list(range(3528)),'complete original node indexing')
own=struct.unpack('<'+'H'*(1<<21),gzip.decompress(files['owners'].read_bytes()))
merge={r[k]:r['merged_id']for r in N for k in ('class_a','class_b')}
edges=[tuple(int(r[k])for k in ('source','target','blue','black'))for r in csv.DictReader(gzip.open(files['edges'],'rt'),delimiter='\t')]
lookup={(i,j):(b,k)for i,j,b,k in edges}
E=json.loads(files['electrical'].read_text())['order8'];Y=E['radical_lifts_mod4'];X=E['higher_layers'][-1]['lifts'][0];U,V=1965,2234
pool=json.loads(files['pool'].read_text());rows=pool['options'];O={r['source']:r for r in rows};cases=json.loads(files['candidates'].read_text())['found']
need(rows==json.loads(files['inventory'].read_text())['records'],'full inventory inherited exactly')
need(len(rows)==len(O)==28 and pool['targets']==[U,V],'28 distinct sources and actual targets')
full=(1<<21)-1;chords=[(i,j)for i in range(8)for j in range(i+2,8)];index={p:i for i,p in enumerate(chords)}
perms=list(itertools.permutations(range(8)))
def matrix(code):
 M=[[0]*8 for _ in range(8)]
 for b,(i,j)in zip(range(27,-1,-1),itertools.combinations(range(8),2)):M[i][j]=(code>>b)&1;M[j][i]=1-M[i][j]
 return M
def trans(M):return [list(z)for z in zip(*M)]
def word(M,p):return sum(M[p[i]][p[j]]<<k for k,(i,j)in enumerate(chords))
def refl(w):return sum(((w>>index[(7-j,7-i)])&1)<<k for k,(i,j)in enumerate(chords))
def tile(w):return {w,full-w,refl(w),full-refl(w)}
def K(M,p):
 kept={frozenset((p[i],p[i+1]))for i in range(7)}
 return [[0 if i==j else M[i][j]if frozenset((i,j))in kept else M[j][i]for j in range(8)]for i in range(8)]
def ham(M,p):return len(set(p))==8 and all(M[p[i]][p[i+1]]for i in range(7))
def lap(v):
 out=[0]*3528
 for i,j,b,k in edges:
  c=(b+k)*(v[i]-v[j]);out[i]+=c;out[j]-=c
 return out
source_reports=[];all_tiles=[];saved={}
for s,r in O.items():
 M=matrix(N[s]['canon_a']);degree=list(map(sum,M));paths=[];aut=[]
 for p in perms:
  if ham(M,p):paths.append(p)
  if all(degree[i]==degree[p[i]]for i in range(8))and all(M[i][j]==M[p[i]][p[j]]for i in range(8)for j in range(i+1,8)):aut.append(p)
 need(aut==[tuple(range(8))],f'complete source automorphism census {s}')
 need(N[s]['class_a']!=N[s]['class_b']and 2*len(paths)==N[s]['fibre'],'complete nonselfconverse path census agrees with original fibre')
 cap=Counter(merge[own[full-word(M,p)]]for p in paths)
 need(all(cap[t]==r['complete_target_capacities'][str(t)]for t in (U,V)),'all-Hamiltonian target capacities')
 p,q=(tuple(r[k]['path'])for k in ('before','after'));fixed=[j for j in range(8)if p[j]==q[j]];free=[j for j in range(8)if p[j]!=q[j]]
 need(len(free)==r['cost']and len(fixed)>=4,'exact changed and fixed clocks')
 domain=[z for z in paths if all(z[j]==p[j]for j in fixed)]
 targets={z:merge[own[full-word(M,z)]]for z in domain}
 need(set(z for z in domain if targets[z]in (U,V))=={p,q},'complete original-target-compatible pinned domain')
 need(len(domain)==(3 if s==1923 else 2),'complete raw pinned domain cardinality')
 if s==1923:need(sorted(targets.values())==[1559,U,V],'unique extra path has target1559')
 start,end,j,c,d=r['pins'];need(j in range(1,6)and len({start,end,c,d})==4,'four distinct shared pin owners')
 need(p[0]==q[0]==start and p[7]==q[7]==end and p[j:j+2]==q[j:j+2]==(c,d)and M[c][d]==1,'actual indexed endpoint and oriented-arc pins')
 roles={start:'first',end:'last',c:'arc0',d:'arc1'};owners={v:('shared',roles[v])if v in roles else('private',s,v)for v in range(8)}
 need(len(set(owners.values()))==8,'literal occurrence map injective')
 domain_report=[]
 for z in domain:
  w=word(M,z);kmat=K(M,z);cmat=trans(M);ckmat=trans(kmat);rz=z[::-1]
  need(ham(kmat,z)and ham(cmat,rz)and ham(ckmat,rz),'all actual K and converse paths Hamiltonian')
  need(K(kmat,z)==M and K(cmat,rz)==ckmat,'actual K involution and converse square')
  words=[word(T,t)for T,t in ((M,z),(kmat,z),(cmat,rz),(ckmat,rz))]
  need(words==[w,full-w,refl(w),full-refl(w)]and len(set(words))==4,'four literal matrix-derived port words form black F/R orbit')
  need(merge[own[w]]==merge[own[refl(w)]]==s,'literal forward/converse source ownership')
  need(merge[own[full-w]]==merge[own[full-refl(w)]]==targets[z],'literal K/converse target ownership')
  colors=lookup[tuple(sorted((s,targets[z])))];need(colors[0]==0 and colors[1]==2*cap[targets[z]],'exact original paired black incidence capacity from full native path census')
  need(all(owners[z[i]]==owners[p[i]]for i in fixed),'every common fixed clock retained in actual owner map')
  all_tiles.append(set(words));domain_report.append({'path':list(z),'target':targets[z],'forward_word':w,'resource_id':min(words),'closed_ports':sorted(words)})
 for name,z in (('before',p),('after',q)):need(z in domain and word(M,z)==r[name]['mask']and targets[z]==r[name]['target'],'supplied endpoint reconstructed independently')
 need((X[s]-X[U])%8==(X[s]-X[V])%8,'chosen primitive per-slot phase')
 label=[sum((Y[i][s]%2)<<i for i in range(4)),X[s]%128];need(pool['labels'][str(s)]==label,'literal electrical source label')
 saved[s]={'paths':(p,q),'fixed':fixed,'cost':len(free),'domain':domain}
 source_reports.append({'source':s,'all_hamiltonian_paths':len(paths),'automorphisms':[list(z)for z in aut],'source_matrix':M,'target_capacities':{str(t):cap[t]for t in (U,V)},'pins':r['pins'],'owner_map':{str(v):list(o)for v,o in owners.items()},'fixed_clocks':fixed,'changed_clocks':free,'complete_pinned_domain':domain_report})
need(len(all_tiles)==57 and len(set.union(*all_tiles))==sum(map(len,all_tiles))==228,'all57 raw pinned alternatives have disjoint closed resources')
# All-dyadic saturation is tied to the independently audited ORIGINAL graph,
# not inferred from coincidental divisibility of the selected debt by 512.
prior=json.loads(files['dyadic_audit'].read_text());F=prior['filtration'];prior_hash={z['path']:z['sha256']for z in prior['source_hashes']}
need(prior['status'].startswith('PASS')and prior['checks']==530454,'prior independent dyadic audit certificate')
for k in ('edges','electrical'):need(prior_hash[str(files[k].relative_to(a.repo))]==hashes[k],'independent audited filtration uses same original bytes')
need(F['valuation_ranks']==[3358,4,0,0,0,0,0,1]and F['paired_charge_type']=='(Z/2)^4 + Z/128','independent original exponent256 and complete paired charge bank')
AX=lap(X);energy=sum(x*y for x,y in zip(X,AX))
need(all(v%256==0 for v in AX)and any(v%512 for v in AX),'raw saved X genuinely harmonic256 but not512')
need(energy==6819178752 and energy//256==26637417 and energy//256%2==1,'direct highest divided energy kills unique surviving dyadic mode')
for y in Y:need(all(v%4==0 for v in lap(y)),'original saved Y actual harmonic4')
need(X[U]==X[V]==146 and all((y[U]-y[V])%2==0 for y in Y),'actual paired targets agree under all marked charge coordinates')
ordinary=[0]*3528;ordinary[U]=4;need(all(v%8==0 for v in lap(ordinary)),'hostile4eU is actual ordinary harmonic8')
need((ordinary[U]-ordinary[V])%8==4,'ordinary per-port phase not preserved')
images=json.loads(files['images'].read_text());need(images['modulus']==256 and images['input_sha256']==hashes['edges'],'image witness original matrix and modulus')
verified={}
for z in images['solutions']:
 potential=z['potential'];charge=dict(z['charge']);corr=dict(z['integer_correction_nonzero'])
 need(len(potential)==3528 and potential[-1]==0 and all(type(v)is int and 0<=v<256 for v in potential),'full grounded integer potential schema')
 need(len(charge)==len(z['charge'])and len(corr)==len(z['integer_correction_nonzero']),'no duplicate image support keys')
 need(all(0<=i<3528 and type(v)is int and v!=0 for i,v in list(charge.items())+list(corr.items())),'image support entries legitimate')
 need(sum(charge.values())==sum(corr.values())==0,'zero-sum image charge and correction')
 actual=lap(potential)
 for i,v in enumerate(actual):need(v==charge.get(i,0)+256*corr.get(i,0),'full original-edge integer image equation')
 verified[tuple(sorted(charge.items()))]=potential
need(((U,2),(V,-2))in verified,'explicit common target charge2v original-L image')
chars256=[[64*v for v in y]for y in Y[:4]]+[X];chars8=[[2*v for v in y]for y in Y[:4]]+[X]
for m,chars in ((8,chars8),(256,chars256)):
 for chi in chars:need(all(v%m==0 for v in lap(chi)),'original-edge separator harmonic equation')
reports=[];ledgers={}
for key,case in sorted(cases.items(),key=lambda x:int(x[0])):
 plus,minus=case['plus'],case['minus'];ss=plus+minus;n=len(ss);signs=[1]*len(plus)+[-1]*len(minus)
 need(n==case['slots']==2*int(key)and len(set(ss))==n and not set(ss)&{U,V},'distinct balanced source slots disjoint from targets')
 need(len(plus)==len(minus)and all(len(saved[s]['domain'])==2 for s in ss),'selected complete raw-pinned binary domains without filtering')
 need(tuple(sorted((s,2*t)for s,t in zip(ss,signs)))in verified,'full source charge2q original-L image')
 old=[O[s]['before'if t==1 else'after']for s,t in zip(ss,signs)];new=[O[s]['after'if t==1 else'before']for s,t in zip(ss,signs)]
 base=Counter(z['target']for z in old);need(base==Counter(z['target']for z in new),'actual target margins')
 oldtiles=[tile(z['mask'])for z in old];newtiles=[tile(z['mask'])for z in new]
 need(len(set.union(*oldtiles))==len(set.union(*newtiles))==4*n and len(set.union(*(oldtiles+newtiles)))==8*n,'closed packet injection at both states and their union')
 cost=sum(saved[s]['cost']for s in ss);need(cost==case['owner_clock_cost'],'ordinary owner-clock distance')
 def action(chi,chosen):
  out=Counter()
  for j in chosen:
   s=ss[j]
   for target,w in ((old[j]['target'],-2),(new[j]['target'],2)):
    v=w*(chi[s]-chi[target]);out[s]+=v;out[target]-=v
  return out
 for chi in chars256:
  row=action(chi,range(n));qdot=sum(t*chi[s]for s,t in zip(ss,signs));vdot=chi[U]-chi[V]
  expected={s:2*t*vdot for s,t in zip(ss,signs)};expected[U]=2*qdot;expected[V]=-2*qdot
  need(all(row[v]==z for v,z in expected.items()),'literal selected-edge current equals physical-row ranktwo formula')
  need(all(v%256==0 for v in row.values()),'both endpoints preserve actual character bank')
 accepted8=[];accepted256=[];balanced=0;counts=Counter();ledger=[];tower=Counter()
 for mask in range(1<<n):
  chosen=[j for j in range(n)if mask>>j&1];margin=base.copy()
  for j in chosen:margin[old[j]['target']]-=1;margin[new[j]['target']]+=1
  if margin!=base:continue
  balanced+=1
  rows8=[action(chi,chosen)for chi in chars8]
  if all(all(v%8==0 for v in row.values())for row in rows8):accepted8.append(mask)
  rows256=[action(chi,chosen)for chi in chars256]
  bad=next(((h,row)for h,row in enumerate(rows256)if any(v%256 for v in row.values())),None)
  if bad is None:accepted256.append(mask)
  else:
   h,row=bad;need(row[U]%256!=0 and all(not(v%256)for s,v in row.items()if s not in (U,V)),'every rejected balanced state has actual owned-target separator and zero source current')
   ledger.append([mask,h,U,row[U]%256]);counts[h]+=1
  # Tower coordinates are complete by the independently audited original filtration.
  parity=[sum(signs[j]*Y[h][ss[j]]for j in chosen)%2 for h in range(4)];xd=sum(signs[j]*X[ss[j]]for j in chosen)
  for k in range(1,9):
   if k==1 or(not any(parity)and xd%(1<<min(k-1,7))==0):tower[1<<k]+=1
 need(balanced==math.comb(n,n//2),'complete balanced menu count')
 need(accepted256==[0,(1<<n)-1],'complete H256 accepted relation exactly two endpoints')
 need(len(ledger)==balanced-2 and sum(counts.values())==balanced-2,'complete proper balanced separator ledger')
 need(len(accepted8)==tower[8],'literal H8 currents equal complete audited tower coordinates')
 # Separate connectivity by low-bank moves is useful, but does not alter the all-bank claim.
 reach=[];h8_path=None
 for budget in range(2,n+1,2):
  parent={0:None};todo=deque([0])
  while todo:
   state=todo.popleft()
   for nxt in accepted8:
    if nxt not in parent and(state^nxt).bit_count()<=budget:parent[nxt]=state;todo.append(nxt)
  reach.append({'budget':budget,'reachable':len(parent)})
  if(1<<n)-1 in parent:
   h8_path=[(1<<n)-1]
   while h8_path[-1]:h8_path.append(parent[h8_path[-1]])
   h8_path.reverse();break
 ledger_data=json.dumps(ledger,separators=(',',':')).encode();ledgers[key]=ledger
 reports.append({'slots':n,'plus':plus,'minus':minus,'all_joint_states':1<<n,'balanced_states':balanced,'H8_accepted':len(accepted8),'H256_accepted':len(accepted256),'proper_balanced_rejected':len(ledger),'separator_counts':dict(counts),'separator_ledger_sha256':hashlib.sha256(ledger_data).hexdigest(),'ordinary_owner_clock_cost':cost,'closed_ports_per_endpoint':4*n,'closed_ports_union':8*n,'integer_X_target_debt':2*sum(t*X[s]for s,t in zip(ss,signs)),'full_dyadic_tower_counts':dict(tower),'H8_minimax_slot_budget':budget,'H8_reachability':reach,'H8_attaining_path':h8_path})
report={'status':'PASS','checks':checks,'input_sha256':hashes,'complete_source_options':source_reports,'raw_pinned_domain_distribution':dict(Counter(len(z['domain'])for z in saved.values())),'pool_closed_ports':228,'verified_full_image_equations':len(verified),'raw_X_harmonic256_not512':True,'original_two_primary_exponent':256,'circuits':reports,'scope':'Complete selected-source fixed-common-clock native binary menu. Full original nodewise harmonic banks preserved at every dyadic modulus by audited original exponent256 and full image equations; chosen X per-slot mod8 retained, not all ordinary per-port phases. No complete-atlas rewiring, common-host amalgam, all-order supply, equality of kernels, or odd-primary preservation asserted.'}
a.out.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
if a.separators:a.separators.write_text(json.dumps({'schema':['selection_mask','character_index','target_vertex','current_mod256'],'characters':['64Y0','64Y1','64Y2','64Y3','X'],'circuits':ledgers},separators=(',',':'))+'\n')
print(json.dumps({'status':'PASS','checks':checks,'circuits':reports},sort_keys=True))
