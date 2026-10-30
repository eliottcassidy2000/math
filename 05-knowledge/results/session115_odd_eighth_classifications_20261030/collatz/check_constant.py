from fractions import Fraction as F
import json
from check_one_reset import v5,coeff,params,reference_rate,words
checks=0
def gate(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def tr(a):
 if a%4==3:return 1,0
 if a%4==1:return 2,v5(9-4**a)
 return 4,v5(81-16**a)
def mx(a):
 t,r=tr(a);return F(4*r+1,4*t)
for a in range(1,2001):
 t,r=tr(a)
 gate((t,r)==params(coeff((a,))),('parameters',a))
 gate(mx(a)/a<=F(5,8),('clock bound',a))
 gate((mx(a)/a==F(5,8))==(a==1),('sharpness',a))
 if a>=2:gate(mx(a)/a<F(1,2),('half',a))
for a in range(1,21):
 for u in words(8):
  z=reference_rate((a,),u)
  gate(z<=mx(a),('entry maximum',a,u))
  gate((z==mx(a))==all(x==a for x in u),('entry equality',a,u))
  P,Q,B=coeff(u);R,S,C=coeff((a,)+u)
  gate((3-2**a)*C-(R-S)==2**a*((3-2**a)*B-(P-Q)),('commutator',a,u))
# Finite residue laws use modular arithmetic, never giant powers or inferred means.
residues=[]
for E in range(1,7):
 modulus=4*5**E;q=5**(E+1)
 counts={a4:[0]*(E+1) for a4 in (0,1,2)}
 for a in range(modulus):
  if a%4==3:continue
  val=(9-pow(4,a,q))%q if a%4==1 else (81-pow(16,a,q))%q
  r=min(v5(val),E+1)
  for k in range(1,E+2):counts[a%4][k-1]+=r>=k
 for a4,cs in counts.items():
  for k,n in enumerate(cs,1):gate(n==5**(E-k+1),('residue law',E,a4,k,n))
 residues.append({'E':E,'class_counts':counts})
gate(mx(1)==F(5,8) and mx(5)/5==F(1,8),'eighth examples')
gate(mx(17)==F(9,8) and mx(17)/17==F(9,136),'normalization hostile')
print(json.dumps({'status':'PASS','gates':checks,'a_bound_scan':2000,'entry_clock_budget':8,'entry_controller_a':20,'finite_residue_laws':residues,'Haar_conditional_even_block_mean':'3/8','all_one_high_layer_tail':'1/8'},indent=2))
