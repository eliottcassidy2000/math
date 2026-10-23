#!/usr/bin/env python3
"""Owned multi-cell alignment, centralizers and actual path-arc probes."""
from pathlib import Path
import argparse,collections,gzip,hashlib,itertools,json
from core import paley,double,refine
HERE=Path(__file__).resolve().parent;CHECKS=collections.Counter()
def need(x,msg):
 CHECKS[msg]+=1
 if not x:raise ValueError(msg)
def outmasks(A):return [sum(x<<j for j,x in enumerate(row)) for row in A]
def transpose(A):return [list(x) for x in zip(*A)]
def decode_cells(masks,n=31):return [[v for v in range(n) if m>>v&1] for m in masks]
def automorphisms(A,C,B=None):
 B=A if B is None else B
 return [{u:v for u,v in zip(C,p)} for p in itertools.permutations(C) if all(A[u][v]==B[pu][pv] for u,pu in zip(C,p) for v,pv in zip(C,p))]
def selectors(A,C):
 maps={};modes={}
 for i,X in enumerate(C):
  for j,Y in enumerate(C):
   if i>=j:continue
   q=sum(A[X[0]][v] for v in Y)
   need(all(sum(A[u][v] for v in Y)==q for u in X),'equitable cross block')
   if len(X)==3:
    need(q in(1,2),'triple block normalized unit selector')
    pi={u:next(v for v in Y if A[u][v]==(q==1)) for u in X};mode='delta'
    need(all(A[u][v]==int((pi[u]==v)==(q==1)) for u in X for v in Y),'matching selector reconstructs whole block')
   else:
    need(len(X)==len(Y)==7 and q in(3,4),'seven carrier degree')
    found=[]
    for inward in(0,1):
     pairs=[(u,v) for u in X for v in Y if all((A[u][w] if q==3 else 1-A[u][w])==(A[w][v] if inward else A[v][w]) for w in Y)]
     if len(pairs)==7 and len({v for u,v in pairs})==7:found.append((dict(pairs),inward))
    need(len(found)==1,'one faithful seven row-template transport')
    (pi,mode),=found
    need(all((A[u][w] if q==3 else 1-A[u][w])==(A[w][pi[u]] if mode else A[pi[u]][w]) for u in X for w in Y),'row template reconstructs whole seven block')
   need(set(pi)==set(X) and set(pi.values())==set(Y),'actual selector bijection')
   maps[i,j]=pi;maps[j,i]={v:u for u,v in pi.items()};modes[i,j]=(q,mode)
 return maps,modes

def joint(A,B,cells):
 C=[c for c in cells if len(c)>1];n=len(A)
 if not C:
  need(A==B,'discrete literal alignment');return [list(range(n))],{'root_candidates':1,'internal_survivors':1,'sections':1}
 ma,ta=selectors(A,C);mb,tb=selectors(B,C);need(ta==tb,'retained block template and polarity')
 roots=automorphisms(A,C[0],B);out=[];internal=0
 for g in roots:
  gs=[g]+[{ma[0,i][u]:mb[0,i][g[u]] for u in C[0]} for i in range(1,len(C))]
  if not all(A[u][v]==B[gs[i][u]][gs[i][v]] for i,X in enumerate(C) for u in X for v in X):continue
  internal+=1
  if not all(gs[j][pi[u]]==mb[i,j][gs[i][u]] for (i,j),pi in ma.items() for u in C[i]):continue
  p=list(range(n))
  for g in gs:
   for u,v in g.items():p[u]=v
  need(sorted(p)==list(range(n)),'propagated actual vertex bijection')
  need(all(A[u][v]==B[p[u]][p[v]] for u in range(n) for v in range(n)),'all arcs preserved by common-owner alignment')
  out.append(p)
 sections=0
 for x in C[0]:
  z=[x]+[ma[0,i][x] for i in range(1,len(C))]
  sections+=all(ma[i,j][z[i]]==z[j] for i,j in itertools.combinations(range(len(C)),2))
 return out,{'root_candidates':len(roots),'internal_survivors':internal,'sections':sections}

def partner(A,P):
 n=len(A);S=set(zip(P,P[1:]));B=[[0]*n for _ in range(n)]
 need(sorted(P)==list(range(n)) and all(A[u][v] for u,v in S),'actual owned Hamiltonian path')
 for i,j in itertools.combinations(range(n),2):B[i][j]=A[i][j] if (i,j) in S or(j,i)in S else 1-A[i][j];B[j][i]=1-B[i][j]
 need(all(A[i][j]+B[i][j]-int(i!=j)==int((i,j)in S)-int((j,i)in S) for i in range(n) for j in range(n)),'signed actual path agreement identity')
 return B

def relabel(A,p):
 B=[[0]*len(A) for _ in A]
 for i in range(len(A)):
  for j in range(len(A)):B[p[i]][p[j]]=A[i][j]
 return B

def five_tree():
 p=(1,0,2,3,4);q=tuple(p[(x+1)%5] for x in range(5));C=[list(range(5*i,5*i+5)) for i in range(3)]
 def make(edge2):
  A=[[0]*15 for _ in range(15)]
  for X in C:
   for u,v in itertools.permutations(X,2):A[u][v]=int((v-u)%5 in(1,2))
  for i,j in itertools.combinations(range(3),2):
   for u,v in itertools.product(range(5),repeat=2):
    z=1 if(i,j)==(0,2) else int((p if i==0 else edge2)[u]==v)
    A[5*i+u][5*j+v]=z;A[5*j+v][5*i+u]=1-z
  return A
 A=make(p);B=make(q)
 local=[]
 for old,new in[(p,p),(p,q)]:local.append([(a,b) for a,b in itertools.product(range(5),repeat=2) if all((old[x]+b)%5==new[(x+a)%5] for x in range(5))])
 need(local==[[(0,0)],[(4,0)]],'five-tree individually feasible but conflicting shared middle owner')
 compatible=[]
 for a,b,c in itertools.product(range(5),repeat=3):
  perm=[5*i+(x+[a,b,c][i])%5 for i in range(3) for x in range(5)]
  if all(A[u][v]==B[perm[u]][perm[v]] for u in range(15) for v in range(15)):compatible.append([a,b,c])
 need(not compatible,'all125 joint five-tree rotations rejected')
 need([sum(A[u]) for u in range(15)]==[8]*5+[7]*5+[6]*5,'five-tree cells intrinsically degree distinguished')
 need([sum(B[u]) for u in range(15)]==[8]*5+[7]*5+[6]*5,'five-tree same equitable quotient')
 # Internal R5 automorphisms are precisely the5 rotations.
 for X in C:need(len(automorphisms(A,X))==5,'five-cell internal automorphisms exhausted')
 return {'A_outmasks':outmasks(A),'B_outmasks':outmasks(B),'local_rotation_solutions':local,'global_rotation_solutions':compatible,'support':[[0,1],[1,2]],'scores':[8,7,6]}

def main():
 p=argparse.ArgumentParser();p.add_argument('--inherited',type=Path,default=Path('/tmp/session41-inherited31.json.gz'));p.add_argument('--mod8',type=Path,default=Path('/tmp/session41-mod8.json'));p.add_argument('--write',action='store_true');args=p.parse_args()
 raw=gzip.decompress(args.inherited.read_bytes()).replace(b'\r\n',b'\n');R=json.loads(raw)
 need(hashlib.sha256(raw).hexdigest()=='2ac21f3fe2bbc16b59bf8d9a6c4c62d92460dbfe2aae3271a85afed62e2fdd6b','pinned inherited order31 source and actual ports')
 phase_raw=args.mod8.read_bytes().replace(b'\r\n',b'\n');phase=json.loads(phase_raw)['native']
 need(hashlib.sha256(phase_raw).hexdigest()==R['retained_mod8_receipt_sha256_lf'],'same retained native electrical guard input')
 need(phase['smallest_binary_modulus_rejecting_all_three_switches']==8 and phase['switch_delta_kirchhoff_mod8']=={'1238':2,'2242':6,'2246':6,'2315':2},'native labelled mod8 fingerprint unchanged separate carrier')
 A=double(double(paley(7)));AT=transpose(A);need(outmasks(A)==R['matrix_outmasks'],'literal inherited matrix reconstruction')
 hist=collections.Counter();groups={};carriers={};records=[];section_hist=collections.Counter()
 for row in R['cases']:
  marks=tuple(row['marks']);initial,_,_=refine(A,marks);need(initial==decode_cells(row['initial_cells']),'actual initial endpoint partition')
  cells=decode_cells(row.get('after_overlap_cells',row['initial_cells']));G,info=joint(A,A,cells);GT,_=joint(AT,AT,cells)
  need(sorted(G)==sorted(GT),'joint gauge centralizer preserved by converse')
  need(all(all(g[v]==v for v in marks) for g in G),'every root alignment fixes actual endpoints')
  expected=21 if row['class']=='seven_C21' else 3 if row['class'] in('triple_C3','seven_to_C3') else 1
  need(len(G)==expected,'complete common-owner group agrees with independent stabilizer census')
  hist[row['class'],len(G)]+=1;groups[marks]=G;carriers[marks]=cells
  if len(G)>1:need(all(any(g[v]!=v for v in range(31)) for g in G if g!=list(range(31))),'nonidentity gauges are actual symmetries')
  if any(len(c)>1 for c in cells):section_hist[max(map(len,cells)),info['sections']]+=1
  records.append({'marks':marks,'class':row['class'],'cell_sizes':list(map(len,cells)),**info,'joint_alignments':len(G),'independent_cell_choices':1 if len(cells)==31 else (3**8 if max(map(len,cells))==3 else 21**4)})
 need(hist==collections.Counter({('already_discrete',1):504,('triple_negative_closed',1):126,('triple_C3',3):168,('seven_negative_closed',1):42,('seven_to_C3',3):84,('seven_C21',21):6}),'complete930case group census')
 need(section_hist=={(3,0):420,(7,7):6},'centralizers are not fixed-point sections')
 port_records=[];orbit_mass=0;arc_hist=collections.Counter();occurrence_hist=collections.Counter();window_hostile=None
 for index,row in enumerate(R['actual_ports']['records']):
  P=row['path'];marks=tuple(row['source_endpoints']);G=groups[marks];B=partner(A,P)
  need(outmasks(B)==row['target_outmasks'],'actual target port pinned')
  arc=next(((u,v) for u,v in zip(P,P[1:]) if sum(g[u]==u and g[v]==v for g in G)==1),None) if len(G)>1 else None
  need(len(G)==1 or arc is not None,'one actual path arc kills residual shared alignment')
  if len(G)==21:need(arc[0]%8<7 and arc[1]%8<7 and arc[0]%8!=arc[1]%8,'exceptional arc carries two distinct F7 owners')
  surviving=0;path_orbit=set()
  for g in G:
   Q=[g[v] for v in P];path_orbit.add(tuple(Q));BB=partner(A,Q)
   need(BB==relabel(B,g),'source gauge and target port transform together')
   need(Q[0]==P[0] and Q[-1]==P[-1],'actual endpoints retained in port orbit')
   surviving+=arc is None or(g[arc[0]]==arc[0] and g[arc[1]]==arc[1])
  need(len(path_orbit)==len(G),'group acts freely on full actual path')
  need(surviving==1,'indexed shared actual arc identifies one source-automorphic port alignment')
  occurrence=None;best=None
  if arc is not None:
   edge_sets=[set(zip(Q,Q[1:])) for Q in path_orbit]
   occurrence=sum(arc in EE for EE in edge_sets)
   best=min(sum(e in EE for EE in edge_sets) for e in zip(P,P[1:]))
   occurrence_hist[len(G),occurrence,best]+=1
  if index==697:
   edges=list(zip(P,P[1:]));banks=[(i,j) for i,j in itertools.combinations(range(30),2) if sum(edges[i] in EE and edges[j] in EE for EE in edge_sets)==1]
   windows=[]
   for length in range(1,6):
    counts=[sum(tuple(P[j:j+length+1]) in {Q[t:t+length+1] for t in range(31-length)} for Q in path_orbit) for j in range(31-length)]
    windows.append({'arc_length':length,'minimum_survivors':min(counts),'first_best_start':counts.index(min(counts))})
   need(best==3 and len(banks)==13 and(0,4)in banks,'sharp two unindexed edge bank in actual order31 orbit')
   need([x['minimum_survivors'] for x in windows]==[3,3,3,3,1],'unindexed local windows require five arcs on actual hostile')
   window_hostile={'control_index':index,'path':P,'orbit_size':len(G),'single_edge_minimum_survivors':best,'isolating_two_edge_index_pairs':banks,'explicit_bank':[[23,0],[16,1]],'contiguous_windows':windows}
  need(partner(AT,list(reversed(P)))==transpose(B),'full converse path transport')
  orbit_mass+=len(G);arc_hist[len(G)]+=1
  port_records.append({'index':index,'endpoints':marks,'root_alignment_choices':len(G),'selected_owned_arc':arc,'actual_arc_position':None if arc is None else P.index(arc[0]),'remaining_indexed_choices':surviving,'unindexed_selected_arc_survivors':occurrence,'best_unindexed_single_arc_survivors':best})
 need(orbit_mass==1517 and arc_hist=={1:686,3:242,21:5},'actual933port aligned orbit census')
 # A different triple system still needs the tree chirality guard.
 def two_triples(perm):
  M=[[0]*6 for _ in range(6)]
  for off in(0,3):
   for x in range(3):M[off+x][off+(x+1)%3]=1
  for x,y in itertools.product(range(3),repeat=2):M[x][3+y]=int(perm[x]==y);M[3+y][x]=1-M[x][3+y]
  return M
 triple_bad,triple_info=joint(two_triples((0,1,2)),two_triples((1,0,2)),[[0,1,2],[3,4,5]])
 need(not triple_bad and triple_info['root_candidates']==3 and triple_info['internal_survivors']==0,'different triple systems require tree chirality compatibility')
 # Native order31 same-local-orbits hostile, with an actual path unchanged.
 marks=(0,7);cells=carriers[marks];C=[c for c in cells if len(c)>1];maps,_=selectors(A,C);P=next(x['path'] for x in R['actual_ports']['records'] if tuple(x['source_endpoints'])==marks)
 ci={v:i for i,X in enumerate(C) for v in X};used={tuple(sorted((ci[u],ci[v]))) for u,v in zip(P,P[1:]) if u in ci and v in ci and ci[u]!=ci[v]}
 i,j=next((i,j) for i,j in itertools.combinations(range(1,8),2) if(i,j)not in used);X,Y=C[i],C[j];rot={u:next(v for v in Y if A[u][v]) for u in Y};q=sum(A[X[0]][v] for v in Y);AA=[r[:] for r in A];pi={u:rot[maps[i,j][u]] for u in X}
 for u,v in itertools.product(X,Y):AA[u][v]=int((pi[u]==v)==(q==1));AA[v][u]=1-AA[u][v]
 need(refine(A,marks)[0]==refine(AA,marks)[0],'hostile retains actual endpoint equitable cells')
 need(all(sum(r)==15 for r in AA),'hostile remains regular31')
 Gbad,_=joint(A,AA,cells);need(not Gbad,'local block orbit equality fails joint holonomy intertwining')
 BM=partner(AA,P);need(partner(transpose(AA),list(reversed(P)))==transpose(BM),'hostile has an actual converse-compatible Hamiltonian port')
 spectrum=lambda M:dict(sorted(collections.Counter(sum(M[u][w]*M[v][w] for w in range(31)) for u,v in itertools.combinations(range(31),2)).items()))
 need(spectrum(A)=={7:465} and spectrum(AA)=={6:36,7:393,8:36},'ordinary and converse-merged nonisomorphism hostile spectrum')
 hostile={'marks':marks,'changed_cell_edge':[i,j],'cell_X':X,'cell_Y':Y,'path':P,'A_outmasks':outmasks(A),'twisted_outmasks':outmasks(AA),'original_common_out_spectrum':spectrum(A),'twisted_common_out_spectrum':spectrum(AA),'root_intertwiners':len(Gbad),'target_outmasks':outmasks(BM)}
 result={'status':'PASS','scope':'All930 actual endpoint carriers; shared internal-chart plus cycle-intertwiner criterion;933 inherited actual ports and1517 source-automorphic alignments. Not all Hamiltonian paths/order31 wiring and not arbitrary target alignments.','input_sha256':{'order31_receipt_lf':hashlib.sha256(raw).hexdigest(),'mod8_receipt_lf':hashlib.sha256(phase_raw).hexdigest()},'group_census':[[*k,v] for k,v in sorted(hist.items())],'section_census':[[*k,v] for k,v in sorted(section_hist.items())],'endpoint_carriers':records,'actual_port_controls':port_records,'port_orbit_mass':orbit_mass,'unindexed_arc_occurrence_histogram':[[*k,v] for k,v in sorted(occurrence_hist.items())],'unindexed_window_hostile':window_hostile,'native31_twisted_chord_hostile':hostile,'five_tree_hostile':five_tree(),'triple_tree_orientation_hostile':triple_info,'checks':dict(sorted(CHECKS.items()))}
 rawout=(json.dumps(result,sort_keys=True,indent=2)+'\n').encode();dest=HERE/'receipt.json'
 if args.write:dest.write_bytes(rawout)
 else:need(dest.read_bytes()==rawout,'frozen exact receipt reproduction')
 print(json.dumps({'status':'PASS','checks':sum(CHECKS.values()),'endpoint_cases':len(records),'actual_ports':len(port_records),'port_orbit_mass':orbit_mass,'one_indexed_arc_suffices':True,'sha256':hashlib.sha256(rawout).hexdigest()},indent=2))
if __name__=='__main__':main()
