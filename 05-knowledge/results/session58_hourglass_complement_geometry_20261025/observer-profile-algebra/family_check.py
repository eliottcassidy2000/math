from fractions import Fraction as F
from math import comb
import json
from pathlib import Path

def require(v,m):
 if not v:raise RuntimeError(m)
rows=[]
for r in range(1,20):
 # 2 integral d(1-d)(1-2d)^r, exact binomial expansion.
 value=2*sum((F(comb(r,k)*(-2)**k,k+2)-F(comb(r,k)*(-2)**k,k+3)) for k in range(r+1))
 require(value==0 if r%2 else value>0,'integral parity')
 rows.append({'r':r,'integral':str(value),'valid_Royle_even_family':r>=3 and r%2==1})
Path(__file__).with_name('family_certificate.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
