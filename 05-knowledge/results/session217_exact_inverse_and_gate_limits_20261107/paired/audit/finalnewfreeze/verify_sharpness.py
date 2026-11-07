"""Exact finite controls of the C3 two-counter sharpness argument."""
import json
previous_max=-1;layers=0;pairs=0
for n in range(2,161):
 values={}
 for q in range(1,n):
  z=n-q;x=2**(n+3)-2**(q+2)-2**(z+1)+3
  # Stable tuple h=3,a=rho=2,b=4,d=1, j=0, via unexpanded sum.
  u=2**q;v=2**z
  h=3+2*(u-1)+2*(v-1)+4*(u-1)*(v-1)
  a=u*(2+4*(v-1))
  if x!=h+a:raise RuntimeError('triangle formula')
  if not 6*2**n-1<=x<8*2**n+3:raise RuntimeError('layer interval')
  values[x]=values.get(x,0)+1;pairs+=1
 if max(values.values())>2 or 2*len(values)<n-1:raise RuntimeError('binary multiplicity')
 if min(values)<=previous_max:raise RuntimeError('layer separation')
 previous_max=max(values);layers+=1
print(json.dumps({'status':'PASS','layers':layers,'exponent_pairs':pairs,'max_layer':160,'arithmetic':'arbitrary-precision integers','scope':'finite controls, not a replacement for the uniqueness-of-binary-powers proof'},indent=2))
