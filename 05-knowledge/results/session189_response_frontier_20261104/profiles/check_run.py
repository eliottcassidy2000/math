from analyze_models import *
S=json.load(open(D/'affine_banks.json'));vals=iter(map(int,(D/'run_values.txt').read_text().split()));bs=[a for a in range(128) if a.bit_count()<=2];pos={a:i for i,a in enumerate(bs)};checks=0
for case in S:
 for B in case['full_run_coefficients']:
  channels=[B['K0']]
  for i in range(7):channels += [B['K2'][i],B['K3'][i]]
  channels+=B['K22']
  for a in bs:
   for c in bs:
    for ch in channels:
     want=sum(ch[pos[s]*29+pos[t]] for s in subs(a) for t in subs(c));actual=next(vals)
     if actual!=want:raise RuntimeError(('run count mismatch',case['Q'],a,c,actual,want))
     checks+=1
try:next(vals);raise RuntimeError('trailing')
except StopIteration:pass
print('complete independent run coefficient checks',checks)
C=json.load(open(D/'quotients.json'));cc=[list(map(int,z.split())) for z in (D/'census_values.txt').read_text().splitlines()]
if len(C)!=197 or len(cc)!=197:raise RuntimeError('prime count')
for q,v in zip(C,cc):
 if q['S']!=v[1:8] or q['E']!=v[8:]:raise RuntimeError('permutation endpoint mismatch')
print('complete prime representative endpoints checked',197*14,'labeled prime coverage',sum(v[0] for v in cc))
