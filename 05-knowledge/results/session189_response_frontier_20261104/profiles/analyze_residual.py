from analyze_models import *
A=json.load(open(D/'anchor_exploration.json'));M=json.load(open(D/'named_banks.json'));X=readmodels(json.load(open(D/'extra_models.json')),(D/'extra_named_values.txt').read_text());json.dump(X,open(D/'extra_named_banks.json','w'),separators=(',',':'));M+=X
lookup={(tuple(m['Q']),tuple(m['sizes']),m['side']):flatten(m['block_values']) for m in M}
bs=M[0]['block_masks'];out=[]
for a in A:
 for V in a['blind']:
  if V==127:continue
  Q=tuple(a['rows']);v=[]
  for side in range(2):
   f0=lookup[Q,(1,)*7,side];vals=f0.copy()
   for j in range(7):
    if V>>j&1:
     ss=[1]*7;ss[j]=2;fj=lookup[Q,tuple(ss),side];vals=[a+b-c for a,b,c in zip(vals,fj,f0)]
   v.append(vals)
  diff=[(a-b)%4 for a,b in zip(*v)];nz=next((k for k,x in enumerate(diff) if x),None)
  out.append(dict(Q=Q,V=V,nonzero_mod4=sum(bool(x) for x in diff),query=None if nz is None else [bs[nz//29],bs[nz%29]],contrast=diff))
json.dump(out,open(D/'residual_mod4.json','w'),separators=(',',':'));print(json.dumps([{k:v for k,v in r.items() if k!='contrast'} for r in out],indent=2))
