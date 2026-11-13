"""Independent exact five-period/parity controls; standard library only."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json

def require(ok,why):
    if not ok: raise RuntimeError(why)

vals=sorted({F(a,b) for a in range(1,9) for b in range(1,9)})
for pair in product(vals,repeat=2):
    seen=[];p=pair
    for j in range(5):
        require(p not in seen,'premature positive rational period')
        seen.append(p);p=(p[1],(1+p[1])/p[0])
    require(p==pair,'period five return')
for a,b in product(range(1,32,2),repeat=2):
    c=F(1+b,a)
    require(c.numerator%2==0 and c.denominator%2==1,'odd full count obstruction')
result={'status':'PASS','positive_rational_pairs':len(vals)**2,'numerator_and_denominator_bounds':[1,8],'odd_input_pairs':16**2,'period_five_example':[1,1,2,3,2]}
Path(__file__).with_name('exchange_receipts.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
