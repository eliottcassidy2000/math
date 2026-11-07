#!/usr/bin/env python3
import pathlib,json
# Reuse only the independent audit's literal graph/permutation definitions.
s=(pathlib.Path(__file__).parent/'independent.py').read_text();exec(s[:s.index('graphchecks=')])
literal=0
r=(2,4,1)
if profile(r)!=(3,1,1,3,3,1,4):raise RuntimeError('triangle base')
for n in range(2,6):
 for q in range(1,n):
  z=n-q;g=r
  for c in 'Q'*q+'Z'*z:g=add(g,c)
  h,start,end,L,R,d,b=profile(g)
  x=2**(n+3)-2**(q+2)-2**(z+1)+3
  if start!=d or d!=1 or h+L-start!=x:raise RuntimeError('literal compiler base')
  literal+=1
all_values=set();pairs=0;boundchecks=0
for n in range(2,81):
 layer={}
 for q in range(1,n):
  z=n-q;u=2**q;v=2**z
  x=3+2*(u-1)+2*(v-1)+4*(u-1)*(v-1)+u*(2+4*(v-1))
  layer.setdefault(x,[]).append((q,z));pairs+=1
 if len(layer)!=(n+1)//2:raise RuntimeError('exact layer count')
 if min(layer)!=6*2**n-1 or max(layer)>=8*2**n+3:raise RuntimeError('interval')
 if all_values.intersection(layer):raise RuntimeError('layer collision')
 for x,ps in layer.items():
  if len(ps)>2 or x%2!=1:raise RuntimeError('parity or multiplicity')
  if len(ps)==2:
   (q,z),(qq,zz)=ps
   if (qq,zz)!=(z-1,q+1):raise RuntimeError('collision involution')
 all_values.update(layer)
 if len(all_values)!=(n+1)**2//4-1:raise RuntimeError('exact cumulative count')
 N=24*2**n+11
 K=((N-11)//24).bit_length()-1
 if K!=n or max(3*x+2 for x in all_values)>=N or 4*len(all_values)<K*(K-1):raise RuntimeError('N bound')
 boundchecks+=1
print(json.dumps(dict(status='PASS',literal_permutation_compiler_base_checks=literal,exponent_pairs=pairs,integer_N_threshold_checks=boundchecks,exact_layer_count='ceil(n/2)',exact_cumulative_count='floor((K+1)^2/4)-1',max_layer=80),indent=2))
