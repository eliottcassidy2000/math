#!/usr/bin/env python3
"""Independent native570 audit. No producer imports, Git, asserts or external IO."""
from pathlib import Path
from itertools import permutations,product,combinations
from collections import Counter,defaultdict
from array import array
import argparse,csv,gzip,hashlib,json,sys
HERE=Path(__file__).resolve().parent;O=HERE/'inputs/original';P=HERE/'inputs/native'
checks=0
def need(v,s):
 global checks
 checks+=1
 if not v:raise RuntimeError(s)
def get(path):return json.loads(path.read_text())
args=argparse.ArgumentParser();args.add_argument('--out',type=Path,required=True);a=args.parse_args()
for r in get(HERE/'SOURCE_PINS.json')['files']:
 b=(HERE/r['snapshot']).read_bytes();need(len(b)==r['bytes'] and hashlib.sha256(b).hexdigest()==r['sha256'],'snapshot SHA '+r['snapshot'])
N=list(csv.DictReader((O/'metagraph_order8_20261023/merged_nodes.tsv').read_text().splitlines(),delimiter='\t'))
owner=array('H');owner.frombytes(gzip.decompress((O/'metagraph_order8_20261023/owner_map.u16.gz').read_bytes()))
if sys.byteorder!='little':owner.byteswap()
merged={int(r[k]):int(r['merged_id']) for r in N for k in ['class_a','class_b']}
need(len(N)==3528 and len(owner)==2**21,'complete native owner schema')
E=get(O/'metagraph_order8_electrical_20261023/receipt.json')['order8'];X=E['higher_layers'][-1]['lifts'][0];Y=E['radical_lifts_mod4']
C=get(P/'native570.normal.json');S=[1130,1673,1754,2678];T=[570,2234];full=(1<<21)-1
chords=[(i,j) for i in range(8) for j in range(i+2,8)];pairs=list(combinations(range(8),2))
def mat(code):
 code=int(code);B=[[0]*8 for _ in range(8)]
 for k,(i,j) in enumerate(pairs):B[i][j]=(code>>(27-k))&1;B[j][i]=1-B[i][j]
 return B
def transpose(B):return [list(r) for r in zip(*B)]
def word(B,p):return sum(B[p[i]][p[j]]*(1<<k) for k,(i,j) in enumerate(chords))
def reflected(w):return sum(((w>>k)&1)*(1<<chords.index((7-j,7-i))) for k,(i,j) in enumerate(chords))
def orbit(w):return {w,full^w,reflected(w),full^reflected(w)}
def complement(B,p):
 K=transpose(B)
 for i,j in zip(p,p[1:]):K[i][j]=B[i][j];K[j][i]=B[j][i]
 return K
def all_paths(B):
 # Independent iterative stack, not the producer's permutations or recursive DFS.
 out=[];stack=[((),255)]
 while stack:
  p,left=stack.pop()
  if not left:out.append(p);continue
  for v in range(8):
   if left>>v&1 and (not p or B[p[-1]][v]):stack.append((p+(v,),left^(1<<v)))
 return sorted(out)
def isomorphism(B,D):
 groups=defaultdict(list);other=defaultdict(list)
 for v in range(8):groups[sum(B[v])].append(v);other[sum(D[v])].append(v)
 if {k:len(v) for k,v in groups.items()}!={k:len(v) for k,v in other.items()}:return None
 keys=sorted(groups)
 for blocks in product(*(permutations(other[k]) for k in keys)):
  f=[None]*8
  for k,b in zip(keys,blocks):
   for v,w in zip(groups[k],b):f[v]=w
  if all(B[i][j]==D[f[i]][f[j]] for i,j in pairs):return f
 return None
def canon_target(K,t):
 M=mat(N[t]['canon_a'])
 for converse,D in [(False,M),(True,transpose(M))]:
  f=isomorphism(K,D)
  if f is not None:return {'to_canonical_a_or_converse':f,'converse':converse}
 raise RuntimeError('Actual complemented matrix misses canonical target '+str(t))
old=get(O/'session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json')['records']
need(len(old)==28 and len({r['source'] for r in old})==28,'literal old28 source inventory')
used=set();fresh=[];all_raw=[];frames=set();old_sources={r['source'] for r in old}
need(not old_sources.intersection(S+T) and not set(S).intersection(T),'physical source blocks disjoint')
for r in old:
 for side in ['before','after']:
  w=r[side]['mask'];o=orbit(w);need(len(o)==4 and not o&used,'old selected resource is unique');used|=o
need(len(used)==224,'56 old alternatives,224 closed words')
for s,raw_expected,cost in zip(S,[52,31,38,4],[7,6,7,3]):
 r=next(x for x in C['literal_sources'] if x['source']==s);B=mat(N[s]['canon_a']);F=dict(r['fixed_clock_owners']);paths=all_paths(B)
 need(r['matrix']==B,'literal source matrix agrees')
 aut=[f for f in permutations(range(8)) if all(B[i][j]==B[f[i]][f[j]] for i,j in pairs)]
 need(aut==[tuple(range(8))],'rigid literal source')
 raw=[(p,word(B,p),merged[owner[full^word(B,p)]]) for p in paths if all(p[k]==v for k,v in F.items())]
 need(len(raw)==raw_expected==r['raw_domain_size'],'complete raw pinned count')
 need(Counter(t for p,w,t in raw)==Counter({int(k):v for k,v in r['raw_targets'].items()}),'full raw target census')
 chosen=[x for x in raw if x[2] in T]
 need(len(chosen)==2 and {t for p,w,t in chosen}==set(T),'binary ONLY under zero outside-target margins')
 p=next(p for p,w,t in chosen if t==570);q=next(p for p,w,t in chosen if t==2234)
 common=[i for i in range(8) if p[i]==q[i]]
 need(sorted(F)==common and all(F[k]==p[k] for k in common),'all and only equal clocks are pinned')
 need(p[0]!=q[0] and p[-1]!=q[-1],'fresh first and last endpoints really move')
 need(sum(x!=y for x,y in zip(p,q))==cost==r['cost'],'exact native primary indexed owner-clock work')
 need(r['shared_anchor_clock']==common[0],'first listed fixed clock selected as new anchor')
 frame=r['owner_frame'];need(len(set(frame.values()))==8,'injective source-private frame')
 need(frame[str(p[common[0]])]=='new_shared_anchor','actual shared owner remains at its pinned clock')
 private=set(frame.values())-{'new_shared_anchor'};need(not private&frames,'fresh nonanchor names source-private');frames|=private
 out=[]
 for p,w,t in chosen:
  need(merged[owner[w]]==s,'owner-map literal source classification')
  need(all(B[i][j] for i,j in zip(p,p[1:])),'literal forward Hamiltonian path')
  K=complement(B,p);o=orbit(w)
  need(K==next(d['target_matrix'] for d in r['domain'] if d['target']==t),'independent actual target matrix')
  need(word(K,p)==full^w and complement(K,p)==B,'literal complement involution')
  need(complement(transpose(B),p[::-1])==transpose(K),'matrix/path converse square')
  need(word(transpose(B),p[::-1])==reflected(w),'reversed clock converse word')
  need(len(o)==4 and not o&used,'fresh orbit disjoint from full old/fresh reservation');used|=o
  need(sorted(o)==next(d['closed_ports'] for d in r['domain'] if d['target']==t),'exported closed words exact')
  scores=[sum(row) for row in K]
  need(scores==[7-sum(B[i])+int(i==p[0])-int(i==p[-1]) for i in range(8)],'fixed-source endpoint score identity')
  out.append({'path':p,'word':w,'target':t,'canonical_isomorphism':canon_target(K,t),'closed_words':sorted(o)})
 caps=Counter(merged[owner[full^word(B,p)]] for p in paths)
 need({str(t):caps[t] for t in T}==r['complete_target_path_capacities'],'unfiltered full source target capacities')
 fresh.append({'source':s,'fixed_clock_owners':sorted(F.items()),'pinned_paths':len(raw),'extra_raw_paths':len(raw)-2,'raw_targets':dict(sorted(Counter(t for p,w,t in raw).items())),'work':cost,'source_X':X[s],'paired_label':[sum((Y[j][s]%2)<<j for j in range(4)),X[s]%128],'complete_target_capacities':{str(t):caps[t] for t in T},'paths':out})
 all_raw.append(raw)
need([r['paired_label'] for r in fresh]==[[4,35],[15,91],[3,81],[8,45]],'fresh marked labels')
need(sum(r['work'] for r in fresh)==23,'23 owner clocks, not old28 bound22')
need(len(used)==256,'old28+fresh4 reservation64 resources')
raw_joint=0;margin_joint=[]
for z in product(*all_raw):
 raw_joint+=1
 if Counter(t for p,w,t in z)==Counter({570:2,2234:2}):margin_joint.append(tuple(t for p,w,t in z))
need(raw_joint==52*31*38*4 and len(margin_joint)==6,'complete raw Cartesian product leaves exactly six margin2 states')
edges=[];edge={};chars=[X]+[[64*y for y in row] for row in Y[:4]];lapchars=[[0]*3528 for _ in chars]
for r in csv.DictReader(gzip.open(O/'metagraph_order8_20261023/merged_edges.tsv.gz','rt'),delimiter='\t'):
 i,j,b,k=(int(r[c]) for c in ['source','target','blue','black']);edges.append((i,j,b+k));edge[i,j]=(b,k)
 for v,x in zip(lapchars,chars):z=(b+k)*(x[i]-x[j]);v[i]+=z;v[j]-=z
need(len(edges)==438216,'full original edge export')
for j,v in enumerate(lapchars):need(all(x%256==0 for x in v),'full original harmonic equations '+str(j))
for r in fresh:
 for t in T:need(edge[tuple(sorted((r['source'],t)))]==(0,2*r['complete_target_capacities'][str(t)]),'genuine paired-black capacities')
for t in T:need(X[t]==146 and [sum((Y[j][t]%2)<<j for j in range(4)),X[t]%128]==[11,18],'target marked charge and exact X')
def currents(targets,chi):
 d=[0]*3528
 for s,t in zip(S,targets):v=2*(chi[s]-chi[t]);d[s]+=v;d[t]-=v
 return d
start=(570,570,2234,2234);end=(2234,2234,570,570)
base=[currents(start,x) for x in chars];accepted=[];rejected=[]
for ts in sorted(margin_joint):
 witnesses=[]
 for j,x in enumerate(chars):
  d=[v-w for v,w in zip(currents(ts,x),base[j])]
  witnesses.extend({'character':j,'physical_row':i,'residual_mod256':v%256} for i,v in enumerate(d) if v%256)
 if witnesses:rejected.append({'target_choices':ts,'separator':witnesses[0]})
 else:accepted.append(ts)
need(sorted(accepted)==sorted([start,end]) and len(rejected)==4,'primitive exact two-state balanced fibre')
need({sum(int(t==2234)*(1<<i) for i,t in enumerate(ts)) for ts in accepted}==set(C['accepted_masks']),'accepted masks reproduce')
for w in C['proper_balanced_separators']:
 ts=tuple(T[(w['mask']>>i)&1] for i in range(4));j=w['character'];i=w['physical_owner'];z=(currents(ts,chars[j])[i]-base[j][i])%256
 need(z==w['residual_mod256'] and z!=0,'every producer separator independently authentic')
u={1130:1,1673:1,1754:-1,2678:-1};v={570:1,2234:-1};touched=S+T
# Independently build actual incidence Laplacian difference on physical rows.
L0=[[0]*6 for _ in range(6)];L1=[[0]*6 for _ in range(6)]
for ts,L in [(start,L0),(end,L1)]:
 for s,t in zip(S,ts):
  i=touched.index(s);j=touched.index(t);L[i][i]+=2;L[j][j]+=2;L[i][j]-=2;L[j][i]-=2
D=[[b-a for a,b in zip(r0,r1)] for r0,r1 in zip(L0,L1)]
need(D==[[2*(u.get(i,0)*v.get(j,0)+v.get(i,0)*u.get(j,0)) for j in touched] for i in touched],'actual physical rank-two identity including sign')
need(all(sum(r)==0 for r in D),'integer source and target margins preserve degrees')
for chi in chars:
 need(2*sum(u[s]*chi[s] for s in u)%256==0 and 2*sum(v[t]*chi[t] for t in v)%256==0,'full paired code annihilation')
 need(all(sum(D[i][j]*chi[touched[j]] for j in range(6))%256==0 for i in range(6)),'all physical rows vanish modulo256')
debt=[sum(D[i][j]*X[touched[j]] for j in range(6)) for i in range(6)]
need(debt==[0,0,0,0,256,-256],'integer X NOT preserved')
need(sum(X[s] for s in S[:2])==382 and sum(X[s] for s in S[2:])==254,'source integer X sums')
for r in fresh:
 s=r['source'];need((X[s]-X[570])%8==(X[s]-X[2234])%8,'selected ordinary per-slot X phase mod8')
# Pairwise primitive failure is also a cheap explicit hostile control.
for i in range(2):
 for j in range(2,4):
  ts=list(start);ts[i]=2234;ts[j]=570
  need(any(any((a-b)%256 for a,b in zip(currents(ts,x),r)) for x,r in zip(chars,base)),'every proper one-each exchange fails actual harmonic bank')
bank=get(O/'compatible_charged_helpers_20261023/supply/selected_supply.json');helpers={h['source']:h for h in bank['helpers']};main=bank['banks'][bank['main_bank_index']]['source_ids'];helperids=main+sorted(set(helpers)-set(main))
need(len(main)==40 and len(helperids)==53 and not set(helperids)&set(S+T) and not set(helperids)&old_sources,'old/fresh/incoming physical sources and targets disjoint')
for k,s in enumerate(helperids):
 B=mat(N[s]['canon_a']);h=helpers[s]
 for side in ['before','after']:
  r=h['selected'][side];p=r['path'];w=word(B,p);o=orbit(w)
  need(all(B[i][j] for i,j in zip(p,p[1:])) and w==r['mask'],'incoming literal Hamiltonian path')
  need(merged[owner[w]]==s and merged[owner[full^w]]==r['target'],'incoming source/target classification')
  need(len(o)==4 and not o&used,'incoming resource disjointness');used|=o
 if k==39:need(len(used)==576,'old28+fresh4+main40:144 alternatives576 words')
need(len(used)==680,'old28+fresh4+all53:170 alternatives680 words')
# Authenticate obstruction targets directly from original canonical matrices.
scores={t:sorted(map(sum,mat(N[t]['canon_a']))) for t in [570,1965,2234]}
need(scores=={570:[2,2,2,3,4,5,5,5],1965:[2,3,3,3,3,4,5,5],2234:[2,2,3,4,4,4,4,5]},'original target score obstruction')
need(sum(d*d for d in scores[570])==112 and all(sum(d*d for d in scores[t])==106 for t in [1965,2234]),'converse-invariant square-score separation')
result={'status':'PASS: independent native570 finite certificate, with 2-primary-exponent wording correction','checks':checks,'fresh_slots':fresh,'raw_joint_configurations':raw_joint,'margin_two_configurations':len(margin_joint),'accepted_target_choices':accepted,'rejected_target_choices':rejected,'physical_rows':touched,'delta_L':D,'integer_X_current_delta':dict(zip(touched,debt)),'work':23,'target_score_multisets':scores,'resource_counts':{'fresh_reserved':8,'old28_plus_fresh_reserved':64,'with_main40_reserved':144,'with_all53_reserved':170,'fresh_occupied':4,'old28_plus_fresh_occupied':32,'with_main40_occupied':72,'with_all53_occupied':85},'dyadic_scope':'Complete all-dyadic conclusion uses independently audited original 2-primary classification, not the unsupported claim that the full critical-group exponent is256. Original-bank inclusion only.','boundaries':['Old28 slots fixed; fresh endpoints move while common clocks are fixed.','Zero fresh margins outside570/2234 are essential to binary filtering.','Initial and final570 population is2; not entry fromzero.','Exact integer X current is not preserved.','No global minimum work23 or arbitrary72-slot fibre connectivity claim.','No order10 harmonic transfer or Collatz consequence.','One shared new anchor, private other local owners; no common-host amalgam.','Compatible unchanged remainder must avoid every reserved alternative.']}
a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','checks':checks,'output':str(a.out)}))
