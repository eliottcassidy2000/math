from fractions import Fraction as F
import json
from check_one_reset import coeff,params,reference_rate,words,orbit,modpair,v5
checks={}
def gate(x,label):
 checks[label]=checks.get(label,0)+1
 if not x:raise RuntimeError(label)
def modulus(a,L,e):
 t,r=params(coeff((1,)*L+(a,)));h=v5(2**a-2)
 if t>1:
  if e<=h:return 5**max(0,min(r,e)-1)
  return 2*5**max(0,e-h-1,min(r,e)-1)
 return 2*5**max(0,min(r,e)-1,e+r-h-1)
def formula(a,L):
 t,r=params(coeff((1,)*L+(a,)));z=F(4*r+1,4*t);e=1
 while True:
  g=modulus(a,L,e)
  if g>L:break
  z+=F(L//g,t*5**max(0,e-r));e+=1
 return z
def phase_good(a,L,j):
 e=1
 while True:
  g=modulus(a,L,e)
  if g>L:return True
  if j%g>L%g:return False
  e+=1
def max_entry(a,L,u):
 w=(1,)*L+(a,)
 for j in range(L+1):
  full=(1,)*j+u
  if len(full)%len(w)==0 and full==w*(len(full)//len(w)):return phase_good(a,L,j)
 return False
def radices(a,L):
 t,r=params(coeff((1,)*L+(a,)));h=v5(2**a-2)
 if h==0 or t==1:return [2]
 if r<=h:return [5]*(r-1)+[2]
 return [5]*(h-1)+[10]
def digit_good(a,L,j):
 initial=radices(a,L);i=0
 while L or j:
  base=initial[i] if i<len(initial) else 5
  if j%base>L%base:return False
  L//=base;j//=base;i+=1
 return True
def phase_count(a,L):
 initial=radices(a,L);i=0;p=1
 while L:
  base=initial[i] if i<len(initial) else 5
  p*=L%base+1;L//=base;i+=1
 return p
rows=[]
for a in range(2,22):
 for L in range(16):
  w=(1,)*L+(a,);mx=formula(a,L)
  gate(mx==reference_rate(w,()),'exact_max_formula')
  for u in words(7):
   z=reference_rate(w,u)
   gate(z<=mx,'entry_upper_bound')
   gate((z==mx)==max_entry(a,L,u),'complete_equality')
  count=0
  for j in range(L+1):
   ok=phase_good(a,L,j);count+=ok
   gate(ok==digit_good(a,L,j),'mixed_radix_equivalence')
  gate(count==phase_count(a,L),'phase_count')
  if a in (2,3,5,17) and L in (0,2,8,10):rows.append({'a':a,'L':L,'mean':str(mx),'clock_mean':str(mx/(L+a)),'initial_radices':radices(a,L),'maximizing_phases':count})
# Independent cyclic subgroup membership, including h2,h3,h4 reset factors.
for a in list(range(2,22))+[41,101,501]:
 for L in list(range(12))+[24,124]:
  w=(1,)*L+(a,);previous=1
  for e in range(1,5):
   q=5**e;M=orbit(modpair(w,q),q);g=modulus(a,L,e)
   gate(g%previous==0,'nested_moduli');previous=g
   x,y=1,0;ma,ba=modpair((1,),q)
   for j in range(2*5**(e-1)):
    gate(((x,y) in M)==(j%g==0),'literal_intersection_membership')
    x,y=ma*x%q,(ma*y+ba)%q
for a in [2,3,5,9,17,21,101]:
 for L in range(200):
  count=0
  for j in range(L+1):
   ok=phase_good(a,L,j);count+=ok
   gate(ok==digit_good(a,L,j),'extended_radix_equivalence')
  gate(count==phase_count(a,L),'extended_phase_count')
A=orbit(modpair((1,),5),5);M=orbit(modpair((1,1,5),5),5)
gate(A==M=={(1,0),(4,3)},'nonunit_intersection_hostile')
gate(reference_rate((1,1,5),())==F(69,40),'nonunit_mean_hostile')
gate(formula(5,8)==F(45,8) and phase_count(5,8)==9,'radix10_example')
print(json.dumps({'status':'PASS','gates':checks,'total_gates':sum(checks.values()),'reset_letters_small_entries':'2..21','L_small_entries':'0..15','entry_clock_budget':7,'literal_depths':'1..4; resets2..21,41,101,501; L0..11,24,124','examples':rows},indent=2))
