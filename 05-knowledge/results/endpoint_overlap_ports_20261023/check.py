"""Actual endpoint overlap/holonomy closure and literal Hamiltonian-port controls."""
from pathlib import Path
from collections import Counter
from itertools import combinations,permutations
from hashlib import sha256
import argparse,gzip,json
from core import paley,double,refine

HERE=Path(__file__).resolve().parent
CHECKS=0
def need(ok,message):
 global CHECKS
 CHECKS+=1
 if not ok:raise ValueError(message)
def seed(A,c,selected):
 records=[(c[v],v in selected)for v in range(len(A))];bank={r:i for i,r in enumerate(sorted(set(records)))}
 return refine(A,colors=[bank[r]for r in records])
def shape(cells):return sorted(map(len,cells))
def partition(cells):return sorted(tuple(C)for C in cells)
def masks(cells):return [sum(1<<v for v in C)for C in cells]
def automorphism(A,p,marks):
 need(sorted(p)==list(range(len(A))) and all(p[v]==v for v in marks),'actual endpoint-preserving permutation')
 need(all(A[i][j]==A[p[i]][p[j]]for i in range(len(A))for j in range(len(A))),'literal all-arc automorphism')

def overlap_seven(A,cells):
 big=[C for C in cells if len(C)>1];owners=set();records=[]
 need(len(big)==4 and all(len(C)==7 for C in big),'four seven-cell carrier')
 for i,C in enumerate(big):
  need(all(sum(A[u][v]for v in C)==3 for u in C),'regular internal seven-cell')
  for j,D in enumerate(big):
   if i==j:continue
   q=sum(A[C[0]][w]for w in D)
   need(q in(3,4) and all(sum(A[u][w]for w in D)==q for u in C),'regular cross-seven block')
   rows={u:tuple(A[u][w]if q==3 else 1-A[u][w]for w in D)for u in C}
   need(all(sum(r)==3 for r in rows.values()),'normalized three-subsets')
   row_design=all(sum(x*y for x,y in zip(rows[u],rows[v]))==1 for u,v in combinations(C,2))
   for inward in (False,True):
    templates={v:tuple(A[w][v]if inward else A[v][w]for w in D)for v in D}
    need(all(sum(r)==3 for r in templates.values()),'internal three-subset template size')
    template_design=all(sum(x*y for x,y in zip(templates[u],templates[v]))==1 for u,v in combinations(D,2))
    matches=[(u,v)for u in C for v in D if sum(x*y for x,y in zip(rows[u],templates[v]))==3]
    need(all(rows[u]==templates[v]for u,v in matches),'owned overlap three is equality of literal subsets')
    if len(matches)==1:owners.update(matches[0])
    records.append([i,j,int(inward),q,matches,row_design,template_design])
 return owners,records

def triple_transport(A,cells):
 C=[c for c in cells if len(c)>1];maps={}
 need(len(C)==8 and all(len(c)==3 for c in C),'eight cyclic-cell carrier')
 for D in C:need(all(sum(A[u][v]for v in D)==1 for u in D),'actual oriented cyclic triple')
 for i in range(len(C)):
  for j in range(i+1,len(C)):
   q=sum(A[C[i][0]][v]for v in C[j])
   need(q in(1,2) and all(sum(A[u][v]for v in C[j])==q for u in C[i]),'nonuniform regular triple block')
   pi={u:next(v for v in C[j]if A[u][v]==(q==1))for u in C[i]}
   need(len(set(pi.values()))==3,'literal selected matching bijective')
   maps[i,j]=pi;maps[j,i]={v:u for u,v in pi.items()}
 negatives=[];profile=Counter();common=set(C[0])
 for i in range(1,len(C)):
  for j in range(i+1,len(C)):
   H={x:maps[j,0][maps[i,j][maps[0,i][x]]]for x in C[0]};fixed=[x for x in H if H[x]==x]
   need(len(fixed)in(0,1,3),'affine three-point holonomy type')
   common.intersection_update(fixed)
   profile[len(fixed)]+=1
   if len(fixed)==1:negatives.append([i,j,fixed[0],[[v,H[v]]for v in C[0]]])
 direct=[];boundary=[]
 for x in C[0]:
  assigned=[x]+[maps[0,i][x]for i in range(1,len(C))]
  if all(maps[i,j][assigned[i]]==assigned[j]for i,j in combinations(range(len(C)),2)):
   direct.append(x);boundary.append(assigned[:3])
 need(sorted(common)==direct,'joint based-loop fixed set equals literal all-edge parallel sections')
 need(len(common)in(0,1,3),'global section cardinality on connected affine-three carrier')
 return C,maps,negatives,profile,{'root_fixed_vertices':direct,'selected_boundary_cells':C[:3],'selected_boundary_relation':boundary}

def close(A,marks):
 cells,h,c=refine(A,marks);initial=shape(cells);data={'marks':list(marks),'initial_cells':masks(cells)}
 if len(cells)==31:data['class']='already_discrete';return cells,c,data
 if max(initial)==7:
  owners,records=overlap_seven(A,cells);data['overlap_owners']=sorted(owners);data['overlap_records']=records
  need(len(owners)in(0,4),'observed overlap owner cardinality')
  if owners:need(all(len(set(C)&owners)==1 for C in cells if len(C)==7),'one authenticated owner in each seven-cell')
  cells,h,c=seed(A,c,owners);data['after_overlap_cells']=masks(cells)
  if max(map(len,cells))==7:
   data['class']='seven_C21';return cells,c,data
 C,transport,negative,profile,sections=triple_transport(A,cells);data['holonomy_fixed_count_histogram']=dict(sorted(profile.items()));data['parallel_sections']=sections
 if negative:
  data['selected_negative_holonomy']=negative[0]
  cells,h,c=seed(A,c,{negative[0][2]})
  need(len(cells)==31,'intrinsic negative-holonomy owner closes every vertex')
  data['class']='seven_negative_closed'if max(initial)==7 else 'triple_negative_closed'
 else:
  data['class']='seven_to_C3'if max(initial)==7 else 'triple_C3'
  rotation={u:next(v for v in C[0]if A[u][v])for u in C[0]};perm=list(range(31))
  for i,D in enumerate(C):
   for u in C[0]:
    source=u if i==0 else transport[0,i][u]
    target=rotation[u]if i==0 else transport[0,i][rotation[u]]
    perm[source]=target
  automorphism(A,perm,marks);need(any(perm[v]!=v for v in range(31)),'genuine residual C3, not only refinement failure')
  data['residual_C3_generator']=perm
 return cells,c,data

def insertion_path(A,order):
 P=[]
 for v in order:
  where=next((i for i,w in enumerate(P)if A[v][w]),len(P));P.insert(where,v)
 return P
def fixed_endpoint_path(A,a,b):
 out=[sum(v<<j for j,v in enumerate(row))for row in A];calls=0
 def visit(P,remaining):
  nonlocal calls
  calls+=1
  if calls>200000:raise ValueError('bounded actual path search exhausted')
  if not remaining:return P+[b]if A[P[-1]][b]else None
  choices=[j for j in range(31)if remaining>>j&1 and A[P[-1]][j]]
  choices.sort(key=lambda j:((out[j]&remaining).bit_count(),j))
  for j in choices:
   found=visit(P+[j],remaining^(1<<j))
   if found:return found
  return None
 P=visit([a],((1<<31)-1)^(1<<a)^(1<<b));need(P is not None,'literal prescribed endpoint path found');return P
def partner(A,P):
 n=len(A);forced=set(zip(P,P[1:]));B=[[0]*n for _ in range(n)]
 for i,j in combinations(range(n),2):
  B[i][j]=A[i][j]if(i,j)in forced or(j,i)in forced else 1-A[i][j];B[j][i]=1-B[i][j]
 need(all(B[u][v]for u,v in forced),'target retains actual forced path')
 need(all(A[i][j]+B[i][j]-(i!=j)==int((i,j)in forced)-int((j,i)in forced)for i in range(n)for j in range(n)),'literal aligned source-target signed path')
 return B

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
 A=double(double(paley(7)));rev=[list(row)for row in zip(*A)];cases=[];hist=Counter();by_marks={}
 inherited_matrix=HERE.parent/'metagraph_minimal_signatures_20261023/receipt.json';matrix_raw=inherited_matrix.read_bytes().replace(b'\r\n',b'\n')
 source_rows=next(r['rows']for r in json.loads(matrix_raw)['all_round_all_edge_moment_hostile']if r['name']=='double_double_Paley7')
 need(source_rows==[sum(x<<j for j,x in enumerate(row))for row in A],'exact inherited literal DDP7 adjacency')
 need(all(sum(row)==15 for row in A),'actual DRT order31 scores')
 need(all(sum(A[u][w]*A[v][w]for w in range(31))==7 for u,v in combinations(range(31),2)),'actual DRT common-neighbor law')
 for a in range(31):
  for b in range(31):
   if a==b:continue
   cells,c,record=close(A,(a,b));hist[record['class']]+=1;by_marks[a,b]=(cells,c,record)
   rcells,rc,rrecord=close(rev,(b,a))
   need(record['class']==rrecord['class'],'converse with reversed actual endpoint order retains closure type')
   need(partition(cells)==partition(rcells),'converse closure retains literal unresolved cell memberships')
   if 'C3'in record['class']:
    unresolved=[v for C in cells if len(C)>1 for v in C]
    for v in unresolved:need(len(seed(A,c,{v})[0])==31,'every actual owner in a residual triple closes all cells')
    record['minimal_extra_actual_owners']=1
   elif record['class']=='seven_C21':record['minimal_extra_actual_owners']=2
   else:record['minimal_extra_actual_owners']=0
   record['final_intrinsic_cells']=masks(cells);cases.append(record)
 need(hist==Counter({'already_discrete':504,'triple_C3':168,'triple_negative_closed':126,'seven_to_C3':84,'seven_negative_closed':42,'seven_C21':6}),'complete 930 endpoint census')
 special=[tuple(r['marks'])for r in cases if r['class']=='seven_C21']
 need(set(special)=={(a,b)for a in(7,15,23)for b in(7,15,23)if a!=b},'exact six residual center pairs')
 group=[]
 for mult in(1,2,4):
  for shift in range(7):
   p=[v if v in(7,15,23)else 8*(v//8)+(mult*(v%8)+shift)%7 for v in range(31)];automorphism(A,p,(7,15,23));group.append(p)
 need(len({tuple(p)for p in group})==21,'21 distinct aligned seven-cell symmetries')
 root7=next(C for C in by_marks[special[0]][0]if len(C)==7)
 autos=sum(all(A[root7[i]][root7[j]]==A[p[i]][p[j]]for i in range(7)for j in range(7))for p in permutations(root7))
 need(autos==21,'internal seven-cell has exactly 21 automorphisms')
 for marks in special:
  cells,c,_=by_marks[marks]
  big=[C for C in cells if len(C)>1]
  for C,D in combinations(big,2):need(len({tuple(A[u][v]for u in C)for v in D})==7,'distinct columns force aligned neighboring cell action')
  for v in [v for C in big for v in C]:
   cs,ccdata,cc=seed(A,c,{v});need(shape(cs)==[1]*7+[3]*8,'one genuine owner leaves unavoidable triple symmetry')
   need(sum(p[v]==v for p in group)==3,'one owner has nontrivial order-three stabilizer')
   for D in cs:
    if len(D)>1:
     for w in D:need(len(seed(A,cc,{w})[0])==31,'every second unresolved actual owner completes seven case')
 paths=[insertion_path(A,[(mult*i+shift)%31 for i in range(31)])for mult in range(1,31)for shift in range(31)]
 paths += [fixed_endpoint_path(A,*marks)for marks in((0,1),(0,7),(7,15))]
 distinct_paths=len({tuple(P)for P in paths});need(distinct_paths==932,'distinct literal Hamiltonian paths among 933 controls')
 port_hist=Counter();endpoint_pairs=set();port_records=[]
 for P in paths:
  need(sorted(P)==list(range(31)) and all(A[u][v]for u,v in zip(P,P[1:])),'literal Hamiltonian path with actual vertices')
  marks=(P[0],P[-1]);endpoint_pairs.add(marks);cells,c,data=by_marks[marks];pivots=[]
  while len(cells)<31:
   unresolved={v for C in cells if len(C)>1 for v in C};v=next(v for v in P[1:-1]if v in unresolved);pivots.append(v)
   cells,h,c=seed(A,c,{v});need(len(pivots)<=2,'actual path-owned pivot budget')
  need(len(pivots)==data['minimal_extra_actual_owners'],'path pivots attain proved minimum')
  B=partner(A,P);RB=partner(rev,list(reversed(P)))
  need(RB==[list(row)for row in zip(*B)],'converse target with reversed literal path')
  port_hist[len(pivots)]+=1
  port_records.append({'path':P,'pivots':pivots,'source_endpoints':list(marks),'target_outmasks':[sum(x<<j for j,x in enumerate(row))for row in B]})
 need(set(port_hist)=={0,1,2},'actual controls exercise every pivot budget')
 section_hist=Counter(len(r['parallel_sections']['root_fixed_vertices'])for r in cases if 'parallel_sections'in r)
 need(section_hist==Counter({0:420}),'all 420 authentic eight-triple systems have no global parallel section')
 section_by_class={kind:dict(sorted(Counter(len(r['parallel_sections']['root_fixed_vertices'])for r in cases if r['class']==kind).items()))for kind in hist if any(r['class']==kind and 'parallel_sections'in r for r in cases)}
 inherited=HERE.parent/'metagraph_hall_bott_20261023/receipt.json';phase_raw=inherited.read_bytes().replace(b'\r\n',b'\n');phase=json.loads(phase_raw)['native']
 need(phase['smallest_binary_modulus_rejecting_all_three_switches']==8 and phase['switch_delta_kirchhoff_mod8']=={'1238':2,'2242':6,'2246':6,'2315':2},'retained prior native mod-eight guard, separate carrier')
 result={'status':'FINITE-EXACT all930 ordered marks plus authenticated actual Hamiltonian ports; scoped structural proofs in README',
  'matrix_outmasks':[sum(x<<j for j,x in enumerate(row))for row in A],
  'endpoint_census':dict(sorted(hist.items())),'minimum_additional_owner_census':{'0':672,'1':252,'2':6},
  'parallel_section_census':dict(sorted(section_hist.items())),'parallel_section_census_by_closure_class':section_by_class,
  'actual_ports':{'count':len(paths),'distinct_literal_paths':distinct_paths,'distinct_endpoint_pairs':len(endpoint_pairs),'pivot_histogram':dict(sorted(port_hist.items())),'records':port_records},
  'inherited_matrix_receipt_sha256_lf':sha256(matrix_raw).hexdigest(),
  'retained_mod8_receipt_sha256_lf':sha256(phase_raw).hexdigest(),'cases':cases,'checks':CHECKS}
 raw=(json.dumps(result,sort_keys=True,indent=2)+'\n').encode()
 if args.write:(HERE/'receipt.json.gz').write_bytes(gzip.compress(raw,mtime=0))
 else:need(gzip.decompress((HERE/'receipt.json.gz').read_bytes()).replace(b'\r\n',b'\n')==raw,'frozen endpoint closure receipt')
 print(json.dumps({'checks':CHECKS,'census':dict(hist),'sections':dict(section_hist),'sections_by_class':section_by_class,'actual_ports':len(paths),'distinct_paths':distinct_paths,'distinct_endpoints':len(endpoint_pairs),'actual_pivot_histogram':dict(port_hist),'sha256':sha256(raw).hexdigest()}))
if __name__=='__main__':main()
