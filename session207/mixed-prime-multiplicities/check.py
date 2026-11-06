import json,math,pathlib
from fractions import Fraction as F
P=pathlib.Path(__file__).resolve().parent
for b in range(21):
 for c in range(21):
  k=b+c
  if not k:continue
  denom=math.factorial(b)*math.factorial(c);t=3*b+4*c
  top_second=F(k*(k+1)//2-k*t,denom)
  pairs=F(0)
  if b>=2:pairs+=F(2,math.factorial(b-2)*math.factorial(c))
  if b and c:pairs+=F(4,math.factorial(b-1)*math.factorial(c-1))
  if c>=2:pairs+=F(3,math.factorial(b)*math.factorial(c-2))
  closed=-F(b*b+4*b*c+c*c+3*b+5*c,2*denom)
  assert top_second+pairs==closed
  assert pairs*denom==2*b*(b-1)+4*b*c+3*c*(c-1)
print('PASS exact Fraction coefficient/defect identities for 440 positive exponent pairs b,c<=20.')
j=json.loads((P/'inherited_controls.json').read_text())
expected={'3':(1,-2),'5':(1,-3),'9':(F(1,2),-F(5,2)),'15':(1,-7),'25':(F(1,2),-F(7,2)),'27':(F(1,6),-F(3,2))}
for h,coeffs in expected.items():
 terms=j[h]['core_word_terms'];degree=max(k for k,t,w in terms)
 lead=second=F(0)
 for k,t,w in terms:
  if k==degree:
   lead+=F(w,math.factorial(k));second+=F(w*(k*(k+1)//2-k*t),math.factorial(k))
  elif k==degree-1:second+=F(w,math.factorial(k))
 assert (lead,second)==coeffs,(h,lead,second)
 print('PASS inherited full core-word control H='+h,'top coefficients',str(lead),str(second))
assert j['3']['strong_atoms']=='z**3' and j['5']['strong_atoms']=='z**4'
assert j['9']['strong_atoms']=='2*z**5'
assert j['15']['strong_atoms']=='2*z**6 + 2*z**5'
assert j['25']['strong_atoms']=='z**7 + 2*z**6'
print('PASS inherited complete strong inventory input checks.')
