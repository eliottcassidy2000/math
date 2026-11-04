"""Source-pinned restoration integration with complete named physical controls."""
from pathlib import Path
import hashlib,json,shutil,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent
SOURCES=['controls.py','build_banks.py','physical.cpp']
ARTIFACTS=['controls.json','named_banks.json']
def need(c,m):
 if not c:raise RuntimeError(m)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dependency(pin):return ROOT/pin['path']

def main():
 deps=json.loads((ROOT/'DEPENDENCIES.json').read_text())
 for pin in deps:need(sha(dependency(pin))==pin['sha256'],'reviewed dependency changed: '+pin['path'])
 with tempfile.TemporaryDirectory(prefix='prime8-frontier-integration-') as td:
  D=Path(td)
  for f in SOURCES:shutil.copyfile(ROOT/f,D/f)
  py=[sys.executable]+(['-O'] if sys.flags.optimize else [])
  subprocess.run(py+[str(D/'controls.py')],cwd=D,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
  subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'physical.cpp'),'-o',str(D/'physical')],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  with open(D/'models.in','rb') as a,open(D/'named_values.txt','wb') as b:subprocess.run([str(D/'physical')],stdin=a,stdout=b,stderr=subprocess.PIPE,check=True)
  subprocess.run(py+[str(D/'build_banks.py')],cwd=D,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
  C=json.loads((D/'controls.json').read_text());B=json.loads((D/'named_banks.json').read_text());ncoords=0;physical=[]
  for case in C:
   a,c=case['query'];vals=[]
   for role in ['source','intermediate','target']:
    bank=next(r for r in B if r['case']==case['name'] and r['role']==role);ms=bank['physical_masks'];p={v:i for i,v in enumerate(ms)}
    value=sum(bank['coefficients'][i][j] for i,u in enumerate(ms) if u&a==u for j,v in enumerate(ms) if v&c==v);
    if a in p and c in p:need(value==bank['values'][p[a]][p[c]],'actual basis query matches complete coefficient polynomial')
    vals.append(value);ncoords+=len(ms)**2
   need(vals==case['actual_response_values'],'independent C++ physical response equals Python literal DP');need(vals[1]==vals[2] and (vals[0]-vals[1])%4==2,'exact restoration and nonzero lifted contrast')
   physical.append(dict(name=case['name'],old_vertices=len(case['source']),highest_nodes=case['highest_changed_nodes'],actual_masks=case['query'],response_values=vals,induced_subsets_with_complete_equal_profiles=len(case['hereditary_profiles_source'])))
  for f in ARTIFACTS:need((D/f).read_bytes()==(ROOT/f).read_bytes(),'complete artifact byte mismatch: '+f)
  allsources=['PROOF.md','verify.py','DEPENDENCIES.json']+SOURCES+[p['path'] for p in deps]
  return {'status':'PASS','claim':'Every highest changed quotient of a nontrivial original-module full-response collision is prime of order at least nine.','source_pins':{f:sha(ROOT/f) for f in allsources},'reviewed_dependency_pins':{p['path']:p['sha256'] for p in deps},'complete_named_coordinates':ncoords,'controls':physical,'artifacts':{f:{'bytes':(ROOT/f).stat().st_size,'sha256':sha(ROOT/f)} for f in ARTIFACTS},'scope':'Source-relative original-module XOR, fixed physical names, exact equal-profile restoration of disjoint highest nodes. No unrestricted reconstruction or later-created-module sequences.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
