"""Synthesis link and reviewed-source inventory, not a mathematical reproof."""
from pathlib import Path
import hashlib,json,re
D=Path(__file__).resolve().parent
s=(D/'NOTE.md').read_text()
links=re.findall(r'\]\(([^)]+)\)',s)
local=[x for x in links if not x.startswith(('http://','https://'))]
for x in local:
 if not (D/x).is_file(): raise RuntimeError('missing relative target '+x)
required={
 '../spectrum/PROOF.md':['{empty, 2, 4, 22, 42, 422, 4222}','112 distinct actual candidate words'],
 '../spectrum_four/PROOF.md':['{empty, 2, 4, 22, 42, 422, 4222}','order at least five'],
 '../collatz/RESULT.md':['first descent at step 2410','not context-free'],
 '../oddperfect/PROOF.md':['1 + max(|L_r|, |U_r|)','does **not** demonstrate an arithmetic improvement'],
 '../nearby/NOTE.md':['Omega(n)<=2K+1','Omega(n)<=2K+2','alpha=1/2'],
}
for p,phrases in required.items():
 text=(D/p).read_text()
 for phrase in phrases:
  if phrase not in text:raise RuntimeError('source claim changed: '+p+' '+phrase)
files=[D/'NOTE.md']+[D/x for x in local]
r={'status':'PASS','scope':'Relative links, source fingerprints and selected source phrases only; mathematical review is separate.','relative_links':local,'sha256':{str(p.relative_to(D.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
print(json.dumps(r,indent=2,sort_keys=True))
