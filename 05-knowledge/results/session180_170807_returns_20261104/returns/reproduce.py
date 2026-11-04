"""Rebuild exhaustive suffix/prefix menus, compare receipts, and verify matching."""
from pathlib import Path
import json, subprocess, tempfile, hashlib, sys
D=Path(__file__).resolve().parent
def require(ok,msg):
 if not ok: raise RuntimeError(msg)
def rows(path):return [json.loads(s) for s in Path(path).read_text().splitlines() if s.strip()]
def main():
 subprocess.run([sys.executable,*(['-O'] if sys.flags.optimize else []),str(D/'build.py'),'--check'],check=True,stdout=subprocess.DEVNULL)
 source=json.loads((D/'source_zero_ready.json').read_text())[1]
 inp=str(len(source['suffix_rows']))+'\n'+' '.join(map(str,source['suffix_rows']))+'\n'
 with tempfile.TemporaryDirectory(prefix='marked-return-') as t:
  t=Path(t)
  for name in ['search','prefix']:
   subprocess.run(['g++','-std=c++17','-O3',str(D/(name+'.cpp')),'-o',str(t/name)],check=True)
  r=subprocess.run([str(t/'search'),'2'],input=inp,text=True,cwd=t,capture_output=True,check=True)
  summary=json.loads(r.stdout)
  require(summary==json.loads((D/'search_summary.json').read_text()),'suffix menu summary')
  hits=rows(t/'hits.jsonl');require(hits==rows(D/'hits.jsonl'),'complete lattice menu')
  r=subprocess.run([str(t/'prefix')],text=True,capture_output=True,check=True)
  got={tuple(map(int,s.split()[1:3])) for s in r.stdout.splitlines()}
  saved={tuple(map(int,s.split()[1:3])) for s in (D/'prefixes.txt').read_text().splitlines()}
  require(got==saved,'complete prefix profiles')
  lattice=[x for x in hits if x.get('lattice')]
  require(len(lattice)==53 and max(x['P2A'] for x in lattice)==198,'all-order prefix bound')
  matches=[x for x in lattice if (x['hA'],x['P2A']) in got]
  require(len(matches)==1,'unique realization')
  x=matches[0]
  require((x['flip_a'],x['flip_b'],x['mark'],x['hA'],x['P2A'],x['target'])==(31,40,2,1,16,512423),'realized profile')
  for target in [512423,512425]:
   require(not any((target-2195*x)>0 and (target-2195*x)%42153==0 for x in range(1,target//2195+1)),'fixed suffix obstruction')
  report={'status':'PASS','summary':summary,'lattice_candidates':len(lattice),'maximum_candidate_P2':198,'unique_match':x,'all_order_prefix_bound':7,'source_pin':'4ad839b20dd280cdcfbadf86cb62b2be3c9615ad','sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [D/'search.cpp',D/'prefix.cpp',D/'source_interface.json',D/'source_zero_ready.json']}}
  (D/'reproduction.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
