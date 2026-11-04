"""Exact coordinate and physical controls for the small-child ambient corollary."""
from pathlib import Path
import hashlib,importlib.util,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('v',ROOT.parent/'covered_response'/'verify.py');v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
E=v.add(v.X,v.Y,v.Z,v.prod(v.X,v.Y),v.prod(v.X,v.Z),v.prod(v.Y,v.Z),v.scale(v.P,3))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 dependencies=[ROOT.parent/'covered_ambient'/'PROOF.md',ROOT.parent/'profile_cone_bounded'/'PROOF.md',ROOT.parent/'profile_cone_bounded'/'verify.py',ROOT.parent/'covered_response'/'PROOF.md',ROOT.parent/'covered_response'/'verify.py',ROOT.parent/'covered_positivity'/'PROOF.md'];pins={str(p.relative_to(ROOT.parent)):sha(p) for p in dependencies}
 cases=[('TT'+str(n),v.tt(n)) for n in range(1,8)]+[('C3',[2,4,1]),('strong4',[8,1,3,6])];records=[]
 for name,rows in cases:
  n=len(rows);q=v.profile(rows);p=v.ONE;gaps=[]
  for h in range(1,3*n+3):
   num={e:c for e,c in v.mul(v.N,p).items() if max(e)<=n};z=h*v.functional(num,h+1,q);v.need(z%3==0,'integer ordered-cover contrast');gaps.append(z//3);p={e:c for e,c in v.mul(p,E).items() if max(e)<=n}
  qq=q+[0];P=sum(j*(qq[j+1]+2*qq[j])*(qq[j]**2-qq[j-1]*qq[j+1]) for j in range(1,n+1));v.need(gaps[0]==-P<0 and all(x<=0 for x in gaps),'full coordinate domination')
  records.append({'name':name,'child_rows':rows,'q':q,'all_augmented_cover_gaps':gaps,'P':P})
 with tempfile.TemporaryDirectory() as td:
  binary=Path(td)/'literal';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'literal.cpp'),'-o',str(binary)],check=True);literal=json.loads(subprocess.check_output([str(binary)]))
 v.need(literal['physical_queries']==72 and len(literal['records'])==36,'complete physical battery')
 v.need(pins=={str(p.relative_to(ROOT.parent)):sha(p) for p in dependencies},'unchanged dependencies')
 return {'status':'PASS','scope':'Equal complete child profiles of common actual order1..7 in a cyclic quotient, arbitrary fixed ambient and profile-preserving target interiors; one lifted physical query has a strict negative source-minus-target integer contrast.','coordinate_controls':records,'literal_controls':literal,'dependency_pins':pins,'source_pins':{p.name:sha(p) for p in [ROOT/'PROOF.md',ROOT/'verify.py',ROOT/'literal.cpp']}}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
