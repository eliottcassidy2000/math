#!/usr/bin/env python3
"""Portable synthesis consistency checks, not a rerun of all research engines."""
import hashlib,json,re,tempfile,shutil,subprocess,sys
from pathlib import Path
from decimal import Decimal
R=Path(__file__).resolve().parent
def require(value,message):
 if not value:raise ValueError(message)
def load(p):return json.loads((R/p).read_text())
def checks():
 pins=load('SOURCE_PINS.json')
 for x in pins:
  b=(R/x['snapshot']).read_bytes();require(hashlib.sha256(b).hexdigest()==x['sha256'] and len(b)==x['bytes'],'source snapshot mismatch')
 # Authenticate reviewed prose and selected receipts against the frozen lane manifests.
 manifest_maps={}
 for p in (R/'inputs/manifests').glob('*/*.json'):
  raw=json.loads(p.read_text());entries=raw.get('files',raw)
  manifest_maps[p.parent.name]=({x['path']:x['sha256'] for x in entries} if isinstance(entries,list) else entries)
 authenticated=0
 for x in pins:
  parts=Path(x['source']).parts
  if not parts[0].startswith('session52-') or parts[0]=='session52-root':continue
  lane=parts[0][len('session52-'):];rel='/'.join(parts[1:])
  if lane in manifest_maps and rel in manifest_maps[lane]:
   require(manifest_maps[lane][rel]==x['sha256'],'source differs from lane manifest: '+x['source']);authenticated+=1
 links=[]
 for name in ['README.md','SYNTHESIS.md']:
  p=R/'inputs/publish'/name
  for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
   if '://' not in link:
    require((p.parent/link.split('#')[0]).is_file(),'missing publish-layout root link: '+link);links.append([name,link])
 d=load('inputs/receipts/session52-balanced-metagraph/normal.json');row=next(r for r in d['counts'] if r['n']==12)
 require(row['balanced_ports']==10188802555033,'balanced ports')
 require(row['internal_ports']==3213825974560==2*row['internal_lines'],'internal pairing')
 require(row['internal_lines']==1606912987280==row['internal_blue_lines']+row['internal_black_lines'],'color total')
 require(row['internal_blue_lines']==3401824 and row['internal_black_lines']==1606909585456,'line colors')
 require(sum(row['endpoint_ports'].values())==row['balanced_ports'],'endpoint total')
 b=load('inputs/receipts/session52-audit-balanced-blue-wiring/independent-normal.json')
 require([b[k] for k in ['SC_classes','internal_lines','loops','boundary_edges','boundary_targets','topology_certificate_edges']]==[229,11452,24,18795,786,455],'blue census')
 require(b['ports']==2*b['internal_lines']+b['boundary_edges']==41699,'blue port mass')
 require(b['internal_lines']-b['SC_classes']+1==11224,'cycle rank')
 t=load('inputs/receipts/session52-balanced-blue-wiring/topology.json')
 require(t['positive_core_edges']==90 and t['core_to_zero_edges']==1857 and t['zero_to_zero_edges']==9505,'core edge split')
 require(90+1857+9505==11452 and len(t['positive_core_nodes'])==20 and t['zero_error_nodes']==209,'core totals')
 require(t['positive_core_component_sizes']==[20],'core connectivity')
 dist=load('inputs/receipts/session52-audit-balanced-blue-wiring/distances-normal.json')
 require(dist['unordered_pair_distance_histogram']==[[1,9219],[2,16886],[3,1]] and dist['diameter']==3,'distance histogram')
 require(9219+16886+1==229*228//2==dist['pair_total'],'pair total')
 require(dist['unique_distance_three_pair']==[63,86] and dist['exception_simple_degrees']==[17,17],'exception pair')
 ns=[set(row[1]) for row in dist['exception_neighbor_lists']]
 require(len(ns[0])==len(ns[1])==17 and not ns[0]&ns[1],'disjoint neighbor sets')
 require([[r[k] for k in ['atlas_index','automorphism_order','internal_degree','blue_presentation_words','boundary_edges']] for r in dist['exception_node_records']]==[[63,5,17,37,20],[86,5,17,35,18]],'exception node records')
 require(37==17+20 and 35==17+18,'exception boundary accounting')
 require(11*903753248+152690==9941438418,'backsolved consistency arithmetic')
 require(4129*28406001782370300553==117288381359406970983337,'first extension factor')
 bits=''.join(str(int(i in [1,2,3,6,7,8,14,16,17])) for i in range(1,18))
 require(bits=='11100111000001011','saturation bits')
 require([len(x.group()) for x in re.finditer(r'0+|1+',bits)]==[3,2,3,5,1,1,2],'run lengths')
 m=load('inputs/receipts/session52-mrb-supersingular/results.json')
 require([(r['p'],r['ogg'],r['supersingular'],r['trace']) for r in m['CM_curve']['rows'] if r['p'] in [17,37]]==[(17,True,False,-6),(37,False,True,0)],'CM counterexamples')
 require(m['CM_curve']['j']==8000 and m['CM_curve']['discriminant']==512,'CM model')
 f=m['periodic_filter'];lo,hi=map(Decimal,f['unit_MRB'].strip('[]').split(','))
 require(abs(lo-Decimal('0.537584423937395341'))<Decimal('5e-19') and abs(hi-Decimal('0.537584423937395341'))<Decimal('5e-19'),'reported decimal')
 require(f['signed_fundamental_period']==30 and f['unsigned_divergence_leading_coefficient']=='4/15','filter period and divergence')
 ex=load('inputs/receipts/session52-blue-existence-probe/normal.json')
 require([[r[k] for k in ['n','anti_invariant_matrices','balanced_matrices','minimum_symmetric_downhill_paths','every_paired_arc_realized']] for r in ex['rows']]==[[2,2,2,1,True],[4,16,8,2,True],[6,512,112,5,True]],'existence finite controls')
 return {'status':'PASS','source_snapshots':len(pins),'manifest_authenticated_sources':authenticated,'root_relative_links':len(links),'links':links,'scope':'Synthesis source consistency, exact arithmetic, and root-link closure; not independent reexecution of all producer engines.'}
def manifest():
 m=load('MANIFEST.json');expected=set(m['files'])
 actual={str(p.relative_to(R)) for p in R.rglob('*') if p.is_file() and p != R/'MANIFEST.json'}
 require(actual==expected,'outer file-set closure')
 for p,h in m['files'].items():require(hashlib.sha256((R/p).read_bytes()).hexdigest()==h,'outer hash '+p)
if __name__=='__main__':
 manifest();result=checks()
 if '--replay' in sys.argv:
  with tempfile.TemporaryDirectory() as td:
   dst=Path(td)/'package';dst.mkdir()
   for p in ['MANIFEST.json',*load('MANIFEST.json')['files']]:
    dest=dst/p;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(R/p,dest)
   logs=[]
   for flags in [[],['-O']]:logs.append(subprocess.check_output([sys.executable,'-I','-B',*flags,str(dst/'check.py')],cwd=td))
   require(logs[0]==logs[1],'normal/optimized replay mismatch')
   require(json.loads(logs[0])==result,'relocated result mismatch')
 print(json.dumps(result,sort_keys=True))
