#!/usr/bin/env python3
"""Independent substitution-only replay of the mod144 finite-clock hostile."""
import json
w='C'
while len(w)<40000:w=''.join({'C':'CA','A':'B','B':'CA'}[c] for c in w)
pos=[i+1 for i,c in enumerate(w) if c=='A'];seen=set();first=None
for k,n in enumerate(pos,1):
    seen.add(n%144)
    if len(seen)==144:first=k;break
count=len({n%144 for n in pos[:5000]})
if count!=95 or first!=11001:raise ValueError('clock mismatch')
print(json.dumps({'status':'PASS','method':'substitution only, no radical floors','first_5000_residue_count':count,'modulus':144,'first_full_coverage_occurrence':first},indent=2))
