"""Compile and replay each declared local search, then independently recount every positive."""
from pathlib import Path
from math import comb,factorial
import hashlib,importlib.util,json,shutil,subprocess,tempfile
HERE=Path(__file__).resolve().parent
p=HERE.parent.parent/'changing_core_returns_20261113/preparation/check.py';spec=importlib.util.spec_from_file_location('search_replay_independent',p);e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
A=(518,388,456,1011,487,455,387,257,1,502)
small=json.loads((HERE/'star_successor_probe.json').read_text(encoding='utf-8'))['cases'][1]['source_rows']
alts=json.loads((HERE/'alternate_1691_inputs.json').read_text(encoding='utf-8'))
helper,_=e.birth(A,3,True);larger=json.loads((HERE/'alternate303_ready3311.json').read_text(encoding='utf-8'))['rows']
cases=[
 ('a2_parent_patch.cpp',(),A,'a2_parent_patch.raw.jsonl',1691),
 ('a2_selective_parent_patch.cpp',(),A,'a2_selective_parent_patch.raw.jsonl',1691),
 ('two_star_search.cpp',(11,3311,0),small,'two_star_3311.raw.jsonl',3311),
 ('two_star_search.cpp',(11,3311,0),alts[0]['rows'],'alternate298_two_star.raw.jsonl',3311),
 ('two_star_search.cpp',(10,1691,183),A,'a2_alt298_two_star.raw.jsonl',1691),
 ('two_star_search.cpp',(12,3311,0),alts[2]['rows'],'alternate303_two_star.raw.jsonl',3311),
 ('two_star_search.cpp',(13,6551,0),larger,'alternate6551_two_star.raw.jsonl',6551),
 ('two_star_search.cpp',(11,1691,351),helper,'helper_alt303_two_star.raw.jsonl',1691)]
compiler=shutil.which('g++');e.need(compiler is not None,'g++ available');records=[];positive_controls=0
with tempfile.TemporaryDirectory(prefix='chan_pak_star_replay_') as directory:
 built={}
 for filename,config,T,rawfile,target in cases:
  key=(filename,config)
  if key not in built:
   exe=Path(directory)/('probe'+str(len(built))+'.exe');args=[compiler,'-O3','-std=c++17']
   if config:args+=['-DORDER='+str(config[0]),'-DTARGET='+str(config[1]),'-DREQUIRED_H='+str(config[2])]
   args+=[str(HERE/filename),'-o',str(exe)];subprocess.run(args,check=True,capture_output=True,text=True);built[key]=exe
  e.need(e.tournament(tuple(T)),'valid source input');n=len(T)
  e.need(factorial(n)//2**(n//2)<2**32,'32bit path count bound via disjoint adjacent pairs')
  out=subprocess.run([str(built[key])],input=' '.join(map(str,T))+'\n',check=True,capture_output=True,text=True).stdout
  actual=[json.loads(line) for line in out.splitlines()];expected=[json.loads(line) for line in (HERE/rawfile).read_text(encoding='utf-8-sig').splitlines()]
  e.need(actual==expected,'exact local search replay '+rawfile);hits=[]
  for row in actual:
   patch=row.get('changed_arcs',row.get('patch'))
   if patch is None:continue
   V=list(T)
   for a,b in patch:V[a]^=1<<b;V[b]^=1<<a
   V=tuple(V);root=row.get('root',9);z=e.state(V,root);expected_state=row.get('parent_state',row.get('state'))
   e.need(z==tuple(expected_state),'independent positive full parent state')
   if 'rows' in row:e.need(V==tuple(row['rows']),'literal endpoint rows')
   W,r=e.birth(V,root,True);zz=e.state(W,r);e.need(zz[0]==target and zz[2]-zz[1]==zz[3]-zz[0]==1,'independent literal unit-ready birth')
   hits.append(dict(root=root,parent_state=z,output_state=zz,output_strong=e.strong(W)));positive_controls+=1;e.tables.cache_clear()
  metadata=[row for row in actual if 'checked' in row or 'inner_stars_checked' in row]
  if not config:e.need([row['checked'] for row in metadata]==[comb(45,w) for w in range(5)],'complete width0..4 universe')
  records.append(dict(raw=rawfile,input_order=n,target=target,source_sha256=hashlib.sha256((HERE/filename).read_bytes()).hexdigest(),raw_sha256=hashlib.sha256((HERE/rawfile).read_bytes()).hexdigest(),metadata=metadata,positive_recounts=hits))
result=dict(status='PASS',searches=len(records),independent_positive_recounts=positive_controls,records=records,scope='Every declared C++ search exactly replayed; all positive profiles and unit-ready births independently recounted in Python; first-hit searches not treated as complete classifications')
(HERE/('search_replay.json' if __debug__ else 'search_replay.optimized.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(dict(status='PASS',searches=len(records),positive_recounts=positive_controls)))
