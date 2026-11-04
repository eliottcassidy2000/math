"""Independent formal and literal audit of the cyclic-child kernel."""
from pathlib import Path
import hashlib,importlib.util,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;SOURCE=ROOT.parent/'cyclic_child_disjoint'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 names=['PROOF.md','verify.py','kernel.py','count.cpp','README.md'];pins={n:sha(SOURCE/n) for n in names}
 need(pins['PROOF.md']=='2c625e1a7c7cff166813d7f7243220b160fe04c3f32418e75bb07f267d9fb1be','frozen proof')
 spec=importlib.util.spec_from_file_location('formal_audit',ROOT/'formal.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 with tempfile.TemporaryDirectory() as td:
  binary=Path(td)/'literal';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'literal.cpp'),'-o',str(binary)],check=True);literal=json.loads(subprocess.check_output([str(binary)]))
 need(literal['physical_query_count']==192,'complete literal query count')
 need([r['local_contrast'] for r in literal['records']]==[60,0,-109584,0],'all local named quartics')
 need(all(r['local_contrast']==r['root_contrast'] for r in literal['records']),'balanced exterior restoration')
 need(pins=={n:sha(SOURCE/n) for n in names},'unchanged producer')
 return {'status':'PASS','scope':'Sections1--3 only: exact cyclic-child disjoint kernel, arbitrary-module same-port zero, and endpoint-balanced exterior restoration. Sections4--6 reviewed separately.','formal_identity':m.RESULT,'literal_controls':literal,'producer_pins':pins,'audit_pins':{p.name:sha(p) for p in [ROOT/'REVIEW.md',ROOT/'audit.py',ROOT/'formal.py',ROOT/'literal.cpp']}}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
