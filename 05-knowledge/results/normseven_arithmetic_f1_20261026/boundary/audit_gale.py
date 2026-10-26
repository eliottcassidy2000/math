"""Independent arithmetic controls for the sibling Gale proof; no producer imports."""
import sys as _output_sys
_output_sys.stdout.reconfigure(newline="\n")
from fractions import Fraction as F
from itertools import product
from math import gcd
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import tempfile

gates=0
def require(ok,label):
    global gates
    gates+=1
    if not ok: raise RuntimeError(label)

rows=[]
for p in (5,7,11,13,17,19):
    all_counts={}
    for delta in range(1,p):
        # Literal triples (x,z,w), without Legendre symbols or histograms.
        affine=sum((w*w-delta*(1-x*x)*(1-z*z)*(1+x+z))%p==0
                   for x,z,w in product(range(p),repeat=3))
        resolved=affine+(p+1)+17*p
        all_counts[delta]=resolved
        points=[(1,v,w) for v,w in product(range(p),repeat=2)]
        points += [(0,1,w) for w in range(p)]+[(0,0,1)]
        conic=sum((w*w-delta*u*v)%p==0 for u,v,w in points)
        require(conic==p+1, ("twisted node conic",p,delta))
    require(sum(all_counts.values())==(p-1)*(p*p+18*p+1),("mean",p))
    rows.append({"p":p,"plus":all_counts[1],"minus":all_counts[p-1],
                 "all_twists":all_counts})
require([(r['plus'],r['minus']) for r in rows if r['p'] in (7,13)] == [(180,172),(408,408)],"witness")

# Independent exact specialization control for inherited nontorsion.
def valuation2(x):
    n,d=x.numerator,x.denominator
    out=0
    while n%2==0:n//=2;out+=1
    while d%2==0:d//=2;out-=1
    return out
u=F(-3)
vals=[]
for n in range(1,9):
    u=(u*u+36)**2/(4*u*(u*u-36))
    vals.append(valuation2(u))
    require(vals[-1]==-2*n,"nontorsion valuation witness")

published=Path(__file__).parent.parent/'gale'/'normal.json'
with tempfile.TemporaryDirectory(prefix='audit-gale-') as tmp:
    reproduced=[]
    for flags,name in (([],"normal.json"),(["-O"],"optimized.json")):
        target=Path(tmp)/name
        subprocess.run([sys.executable,"-B",*flags,str(published.with_name('check.py')),
                        "--output",str(target)],check=True,stdout=subprocess.DEVNULL)
        reproduced.append(target.read_bytes())
    require(reproduced[0]==published.read_bytes(),"producer replay exact")
    require(reproduced[0]==reproduced[1],"producer optimized replay exact")
print(json.dumps({"status":"FINITE-EXACT independent controls; proof audit in AUDIT_GALE.md",
                  "gates":gates,"prime_field_rows":rows,"specialization_v2":vals,
                  "producer_receipt_sha256":hashlib.sha256(published.read_bytes()).hexdigest()},
                 indent=2,sort_keys=True))
