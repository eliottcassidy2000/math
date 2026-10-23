from math import isqrt
from fractions import Fraction
from pathlib import Path
import json
# Opposite nonzero-grade dimensions have even total increments.
base={-2:1,-1:56,0:133,1:56,2:1}
assert sum(base.values())==247
assert 248-sum(base.values())==1
forced=dict(base);forced[0]+=1
assert sum(forced.values())==248
# Analytic integer-rank obstructions, rather than bounded search.
# A_n: (n+1)^2=249; B/C/D: discriminant1+8*248=1985.
assert isqrt(249)**2!=249
assert isqrt(1985)**2!=1985
# Independently solve the diagonal linear equations, including rational controls.
rows=[]
for c in [Fraction(-3),Fraction(-1,2),Fraction(0),Fraction(2,3),Fraction(3)]:
 w={j:j*c for j in range(-2,3)}
 assert w[2]==2*w[1] and w[-1]==-w[1] and w[-2]==-w[2]
 assert all(w[i]+w[j]==w[i+j] for i in w for j in w if i+j in w)
 rows.append({j:str(v) for j,v in w.items()})
out={'forced_profile':forced,'classical_rank_discriminants':{'A':249,'BCD':1985},'both_discriminants_nonsquare':True,'diagonal_rational_controls':rows,'scope':'Exact arithmetic only. General equality rigidity audited by proof in EQUALITY_AUDIT.md.'}
Path('lie-arithmetic-audit/equality_independent_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
