"""Exact last-atom gate; source data from pinned session 155 atom atlas."""
import json
import check as c
from pathlib import Path
bank=json.loads((Path(__file__).parent/'bank65.json').read_text())
c.check(len(bank)==65,'65 labels')
rays=((2,1,1),(1,0,1),(0,1,1)); groups={'low':[],'high':[],'unresolved':[]}
for label, data in sorted(bank.items()):
 c.M[label]=tuple(map(tuple,data['matrix'])); c.G[label]=tuple(data['source_rows']); P=c.M[label]
 v0=c.observed((label,));v1=c.observed(('a',label));v2=c.observed(('a','a',label))
 rows=(v0,tuple(a-b for a,b in zip(v2,v1)),tuple(a-b for a,b in zip(v1,v0)))
 c.check(rows==P,'independent matrix '+label)
 values=[c.vm(c.vm(r,P),c.B) for r in rays]
 c.check(all(a>0 and b>0 for a,b in values),'positive output bulk '+label)
 low=all(11*b<17*a for a,b in values);high=all(5*b>8*a for a,b in values)
 c.check(not (low and high),'disjoint gates')
 groups['low' if low else 'high' if high else 'unresolved'].append(label)
c.check(tuple(len(groups[x]) for x in ('low','high','unresolved'))==(9,44,12),'three class counts')
for k,v in groups.items():print(k,len(v),', '.join(v))
print('65 independently calibrated atom matrices; physical cone rays; PASS')
# Actual source swap outside the seven-macro target family.
V=c.observed(('q5','a'));w=c.observed(('q1',))
c.check(V==(13,5,9) and w==(3,3,3),'actual distinct sources')
c.check(c.vm(w,c.M['b'])==tuple(3*V[i]-6*c.t[i] for i in range(3)),'nontrivial source replacement')
v=c.vm(V,c.D);lhs=c.vm(c.vm(w,c.M['b']),c.D)
c.check(lhs==tuple(3*v[i]+24*c.t[i] for i in range(3)),'mixed physical return')
print('new source swap: q1 b = 3 (q5 a) - 6 t')
print('actual target (q5 a) D =',v,'actual predecessor q1 =',w,'carry 24')
print('PASS')
