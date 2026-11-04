from pathlib import Path
import hashlib,json,subprocess,sys,csv
D=Path(__file__).resolve().parent

def need(x,m):
 if not x:raise RuntimeError(m)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build(q,a,c):
 rows=[0]*29
 for i in range(27):
  bi,pi=divmod(i,3)
  for j in range(27):
   bj,pj=divmod(j,3)
   if (bi==bj and (pj-pi)%3==1) or (bi!=bj and q[bi]>>bj&1):rows[i]|=1<<j
  if a>>bi&1:rows[i]|=1<<27
  else:rows[27]|=1<<i
  if c>>bi&1:rows[i]|=1<<28
  else:rows[28]|=1<<i
 rows[27]|=1<<28
 return rows

def main():
 for pin in json.loads((D/'DEPENDENCIES.json').read_text()):need(digest(D/pin['path'])==pin['sha256'],'reviewed dependency pin '+pin['path'])
 cs=json.loads((D/'input/candidates.json').read_text())['candidates'];ws=json.loads((D/'input/first_witnesses.json').read_text());qs=json.loads((D/'input/queries.json').read_text())
 need(qs==[[0,1],[1,0],[1,2],[1,3],[3,1],[3,3],[3,5],[0,7],[0,18],[40,256],[4,3],[130,128]],'fixed twelve-query battery')
 need(len(cs)==len(ws)==2348 and len({c['node'] for c in cs})==2348,'complete distinct endpoint-zero cases')
 inp=str(len(cs))+'\n'+''.join(' '.join(map(str,[c['node'],9]+qs[w['query_index']]+c['rows']))+'\n' for c,w in zip(cs,ws))
 need(inp==(D/'all_recount.in').read_text(),'complete replay input alignment')
 if '--recount' in sys.argv:
  subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'orbit_physical.cpp'),'-o',str(D/'orbit_physical')],check=True)
  with (D/'all_recount.out').open('wb') as out,(D/'all_recount.log').open('wb') as err:subprocess.run([str(D/'orbit_physical')],input=inp.encode(),stdout=out,stderr=err,check=True)
 output=[list(map(int,l.split())) for l in (D/'all_recount.out').read_text().splitlines()]
 need(len(output)==2348,'all physical orbit counts completed')
 hist=[0]*12;phys=[]
 for c,w,r in zip(cs,ws,output):
  need(c['node']==w['node']==r[0],'same census node')
  q=c['rows'];need(len(q)==9 and all(not(x>>i&1) and all(((x>>j)&1)+((q[j]>>i)&1)==1 for j in range(9) if i!=j) for i,x in enumerate(q)),'quotient tournament')
  need(not any((s+e)%2 for s,e in zip(c['starts'],c['ends'])),'endpoint-zero branch')
  need(0<=w['query_index']<12,'fixed battery witness')
  need(r[1:]==[w['source_residue'],w['target_residue']],'independent physical counts match run counts')
  need((r[1]-r[2])%8==4 and r[1]%2==r[2]%2==1,'nonzero actual Hamilton gap')
  hist[w['query_index']]+=1;a,b=qs[w['query_index']];lift=lambda m:sum(7<<(3*i) for i in range(9) if m>>i&1)
  left=build(q,a,b);right=build([511^v^(1<<i) for i,v in enumerate(q)],a,b)
  for rows in (left,right):
   need(all(not(x>>i&1) and all(((x>>j)&1)+((rows[j]>>i)&1)==1 for j in range(29) if i!=j) for i,x in enumerate(rows)),'full physical tournament')
   need(rows[27]>>28&1 and not(rows[28]>>27&1),'fixed x to y')
   need(sum(((rows[i]>>27)&1)<<i for i in range(27))==lift(a),'actual a mask')
   need(sum(((rows[i]>>28)&1)<<i for i in range(27))==lift(b),'actual c mask')
  for i in range(27):
   for j in range(i+1,27):need(((left[i]>>j)&1)^((right[i]>>j)&1)==(i//3!=j//3),'exact quotient converse with same interiors')
  phys.append(dict(node=c['node'],query_index=w['query_index'],block_query=[a,b],physical_a=lift(a),physical_c=lift(b),source_rows=left,target_rows=right,source_residue=r[1],target_residue=r[2]))
 (D/'physical_records.json').write_text(json.dumps(dict(child_sets=[list(range(3*i,3*i+3)) for i in range(9)],port_names={'x':27,'y':28},records=phys),separators=(',',':'))+'\n')
 small=json.loads((D/'controls.normal.json').read_text());need(small==json.loads((D/'controls.optimized.json').read_text()) and small['status']=='PASS' and small['query_pairs']==580,'normal/optimized independent literal controls')
 for f,h in small['sha256'].items():need(digest(D/f)==h,'small literal control pin '+f)
 census_dir=next(p for p in (D.parent/'session189-c3-census-audit',D.parent/'c3_census_audit') if p.is_dir());cr=json.loads((D/'input/census_receipt.json').read_text());need(digest(census_dir/'receipt.json')==digest(D/'input/census_receipt.json'),'live census receipt equals frozen approval')
 need(cr['status']=='PASS','separate independent census audit');last=cr['summary'][-1]
 need(last['n']==9 and last['merged']==97144 and last['ordinary']==191536 and last['labeled_mass']==2**36 and last['prime_merged']==69445 and last['p_nonzero']==67097 and last['p_zero']==2348,'complete census domain')
 need(cr['sha256']['session189-c3-frontier/candidates.json']==digest(D/'input/candidates.json'),'census proof applies to frozen candidate data')
 with (census_dir/'candidates.tsv').open() as f:
  checked=[dict(node=int(x['node']),code=int(x['code']),rows=list(map(int,x['rows'].split(','))),starts=list(map(int,x['starts'].split(','))),ends=list(map(int,x['ends'].split(',')))) for x in csv.DictReader(f,delimiter='\t')]
 need(checked==cs,'separately canonicalized and literally counted candidates')
 producer=next(p for p in (D.parent/'session189-c3-frontier',D.parent/'c3_frontier') if p.is_dir())
 need((D/'PRODUCER_PINS.json').exists(),'frozen producer pins required')
 if (D/'PRODUCER_PINS.json').exists():
  for name,h in json.loads((D/'PRODUCER_PINS.json').read_text()).items():need(digest(producer/name)==h,'frozen producer file '+name)
 battery=json.loads((producer/'battery.json').read_text());need(battery==[dict(index=i,a=a,c=c,physical_a=sum(7<<(3*j) for j in range(9) if a>>j&1),physical_c=sum(7<<(3*j) for j in range(9) if c>>j&1)) for i,(a,c) in enumerate(qs)],'producer physical battery')
 pw=json.loads((producer/'physical_witnesses.json').read_text());need(len(pw)==len(phys),'producer physical witness count')
 for a,b in zip(phys,pw):
  need(a['node']==b['node'] and a['query_index']==b['query_index'] and a['block_query']==[b['coarse_a'],b['coarse_c']],'producer witness physical alignment')
  for k in ('physical_a','physical_c','source_residue','target_residue'):need(a[k]==b[k],'producer witness field '+k)
  need(a['source_rows']==b['augmented_source_rows'] and a['target_rows']==b['augmented_target_rows'],'producer physical full augmented rows')
 def table(p):
  with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
 bank=table(producer/'residues.tsv');need(len(bank)==2348,'producer complete twelve-query bank rows');nonzero=0
 for a,b in zip(phys,bank):
  need(a['node']==int(b['node']),'bank node alignment');j=a['query_index'];need([int(b['S'+str(j)]),int(b['T'+str(j)])]==[a['source_residue'],a['target_residue']],'independent count equals selected bank entries')
  pairs=[(int(b['S'+str(i)]),int(b['T'+str(i)])) for i in range(12)];need(all(0<=x<8 and 0<=y<8 and x%2==y%2==1 and (x-y)%8 in (0,4) for x,y in pairs),'complete bank residue ranges')
  hits=[i for i,(x,y) in enumerate(pairs) if x!=y];need(hits[0]==j,'bank agrees with retained first witness');nonzero+=len(hits)
 ca=[r for r in table(census_dir/'all_rows.tsv') if r['n']=='9'];pf=table(producer/'filter.tsv');need(len(ca)==len(pf)==97144,'entire prime-filter domain');anchors=[]
 for a,b in zip(ca,pf):
  need(a['node']==b['node'] and a['source_code']==b['code'],'independent all-prime filter physical code');need(bool(int(a['prime']))==(int(b['proper_module'])==0),'independent all-prime classification')
  if int(a['prime']):
   need(all(a[k]==b[k] for k in ('p_mask','starts','ends')),'all prime endpoint coordinates independently agree')
   p=int(a['p_mask'])
   if p:
    j=(p&-p).bit_length()-1;anchors.append(dict(node=int(a['node']),p_mask=p,anchor=j,coarse_a=1<<j,coarse_c=0,physical_a=7<<(3*j),physical_c=0,contrast_mod4=2))
 endpoint=[{k:int(v) for k,v in r.items()} for r in table(producer/'endpoint_witnesses.tsv')];need(endpoint==anchors and len(anchors)==67097,'all endpoint witnesses from independent exact p data')
 names=['REVIEW.md','PRODUCER_PINS.json','input/producer_proof.md','orbit_physical.cpp','literal_physical.cpp','controls.py','verify.py','all_recount.in','all_recount.out','physical_records.json','input/candidates.json','input/queries.json','input/first_witnesses.json','controls.normal.json','controls.optimized.json','DEPENDENCIES.json']
 receipt=dict(status='PASS',approved_original_module_minimum_old_order=28,approved_family_battery_scope='root-containing families on a fixed Q[C3^9] only',scope='All 2348 endpoint-zero prime-nine converse-merged representatives; exactly C3 children; no ambient lift claimed',independent_physical_count_pairs=2348,independent_prime_endpoint_coordinates=1250010,producer_endpoint_witnesses_checked=67097,producer_full_residue_entries_sanity_checked=56352,producer_nonzero_query_pairs=nonzero,physical_old_order=27,augmented_order=29,subset_orbit_states_per_side=1048576,first_witness_histogram=hist,all_residue_differences_mod8=4,literal_controls=small['query_pairs'],census=dict(merged=97144,ordinary=191536,prime_merged=69445,p_nonzero=67097,p_zero=2348,labeled_mass=2**36),sha256={n:digest(D/n) for n in names},census_receipt_sha256=digest(census_dir/'receipt.json'))
 print(json.dumps(receipt,sort_keys=True,separators=(',',':')))
if __name__=='__main__':main()
