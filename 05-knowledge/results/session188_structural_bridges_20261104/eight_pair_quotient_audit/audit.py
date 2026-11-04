"""No producer imports: complete seven-table, alternative eight-extensions, literal ports."""
from pathlib import Path
import hashlib,itertools,json,math,struct,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;SOURCE=ROOT.parent/'eight_pair_quotient'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def digest(b):return hashlib.sha256(b).hexdigest()
def decode(code):
 a=[0]*8;k=0
 for u in range(8):
  for v in range(u+1,8):
   if code>>k&1:a[u]|=1<<v
   else:a[v]|=1<<u
   k+=1
 return a
def relabel(a,p):
 return sum(((a[p[u]]>>p[v])&1)<<k for k,(u,v) in enumerate(itertools.combinations(range(8),2)))
def main():
 sources=['PROOF.md','DATA.md','verify.py','census.cpp','check_endpoint.cpp','literal.cpp','orbits.py'];pins={n:digest((SOURCE/n).read_bytes()) for n in sources}
 with tempfile.TemporaryDirectory() as td:
  folder=Path(td);binary=folder/'audit'
  subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'independent.cpp'),'-o',str(binary)],check=True)
  run=subprocess.run([str(binary)],cwd=folder,capture_output=True,check=True);counts=json.loads(run.stdout)
  seven=(folder/'independent_seven.bin').read_bytes();eight=(folder/'independent_eight.bin').read_bytes()
 need(seven==(SOURCE/'t7.bin').read_bytes(),'all 2097152 seven-table bytes')
 need(eight==(SOURCE/'blind8.bin').read_bytes(),'all 30720 alternative-vertex residual codes')
 need(counts['seven_zero']==50688 and counts['extension_candidates']==6488064 and counts['eight_residual']==30720 and counts['eight_prime_residual']==26880 and counts['endpoint_xor_255_count']==30720,'complete exact counts')
 need(counts['small_physical_queries']==2196 and counts['eighteen_source_mod8']==5 and counts['eighteen_target_mod8']==7,'independent physical counts')
 bank=set(struct.unpack('<30720I',eight));remaining=set(bank);orbits=[];perms=list(itertools.permutations(range(8)))
 while remaining:
  code=min(remaining);a=decode(code);orbit={relabel(a,p) for p in perms};need(orbit<=remaining,'disjoint orbit cover');remaining-=orbit
  op=[r^(255^(1<<i)) for i,r in enumerate(a)];opcodes={relabel(op,p) for p in perms}
  orbits.append({'code':code,'orbit_size':len(orbit),'automorphisms':math.factorial(8)//len(orbit),'converse_minimum':min(opcodes),'rows':a})
 expected=json.loads((SOURCE/'orbits.json').read_text())
 need(len(orbits)==4 and [r['code'] for r in orbits]==[600626,600653,602967,1472106],'complete four orbits')
 for r,s in zip(orbits,expected):
  need(all(r[k]==s[k] for k in ['code','orbit_size','automorphisms','rows']),'independent orbit detail')
  need(r['converse_minimum']==expected[s['converse_orbit']]['code'],'converse orbit pairing')
 need(pins=={n:digest((SOURCE/n).read_bytes()) for n in sources},'source files unchanged')
 return {'status':'PASS','scope':'Independent finite-exact eight-label TT2 quotient-converse theorem, proof review including all-order physical and ambient bridges. No unrestricted even-order or arbitrary tournament reconstruction claim.','complete_census':counts,'complete_banks_sha256':{'t7.bin':digest(seven),'blind8.bin':digest(eight)},'independent_orbits':orbits,'producer_pins':pins,'audit_source_pins':{p.name:digest(p.read_bytes()) for p in [ROOT/'REVIEW.md',ROOT/'audit.py',ROOT/'independent.cpp']}}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
