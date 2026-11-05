from pathlib import Path
from itertools import product
import json,re
P=Path(__file__).resolve().parent

def need(x,s):
 if not x:raise RuntimeError(s)
def f(n):
 prefix=0;value=1
 for digit in bin(n)[2:]:
  if digit=='0':value+=prefix+(prefix&1)
  prefix=2*prefix+int(digit)
 return value
def physical_g(L):
 m=(L-1).bit_length();word=L-1;count=0;allv=(1<<m)-1
 for A in range(1<<m):
  B=allv^A
  before=not A or bool(word>>(A.bit_length()-1)&1)
  after=not B or not bool(word>>((B&-B).bit_length()-1)&1)
  count+=before and after
 return count
def carry(L,z,y,h,e,R,padding=0):
 n=max(L.bit_length(),z.bit_length(),y.bit_length(),1)+padding;c=0
 for i in range(n):
  value=c+((y>>i)&1)-h*((z>>i)&1)-R*((L>>i)&1)
  if i==0:value+=e*(L&1)
  if value&1:return False
  c=value//2
  need(abs(c)<=1+h+R+e,'finite carry bound')
 return c==0
physical=0
for L in range(1,2049):
 need(physical_g(L)==f(L-1),'literal one-apex before/after partitions');physical+=1
forms=0;equal=0
for a in range(3,258):
 for b in range(3,258):
  n=((1<<a)-1)<<b;y=f(n);need(y==(1<<a)*(1<<b)-(1<<a)-(1<<b)+3,'balanced-run formula')
  s,t=min(a,b),max(a,b)
  expected='1'*(s-1)+'0'*(s-1)+'11'if a==b else'1'*(s-1)+'0'+'1'*(t-s)+'0'*(s-2)+'11'
  need(bin(y)[2:]==expected,'exact bit positions')
  accepted=re.fullmatch(r'1{2,}0{2,}11',bin(y)[2:])is not None
  need(accepted==(a==b),'regular track filter');forms+=1;equal+=accepted
relations=0;forced_padding=0
for h,e,R in ((1,1,1),(3,1,3),(5,1,7),(13,2,18),(2,0,3),(17,35,101),(100,0,100)):
 for L in range(1,2049):
  z=f(L-1);y=h*z+R*L-e*(L%2);need(y>=0,'nonnegative sample function')
  for pad in (0,1,7):
   need(carry(L,z,y,h,e,R,pad),'exact fixed-coefficient relation')
   need(not carry(L,z,y+1,h,e,R,pad),'false relation rejected')
   need(not carry(L,z+1,y,h,e,R,pad),'false recovered value rejected');relations+=3
  width=max(L.bit_length(),z.bit_length(),y.bit_length())
  pair=list(zip(format(L,'0%db'%width),format(z,'0%db'%width)))
  if pair[0]==('0','0'):forced_padding+=1
  while len(pair)>1 and pair[0]==('0','0'):pair.pop(0)
  need(int(''.join(x for x,y in pair),2)==L and int(''.join(y for x,y in pair),2)==z,'projection zero-column normalization')
  need(pair[0]!=('0','0'),'canonical projected encoding')
addition=0
for n in range(8192):
 L=n+1;c=1
 for i in range(max(n.bit_length(),L.bit_length())+3):
  q=((n>>i)&1)+c
  need((q&1)==((L>>i)&1),'addition relation');c=q//2
 need(c==0,'terminal carry');addition+=1
result=dict(status='PASS',literal_physical_g_values=physical,big_integer_bit_forms=forms,balanced_filter_hits=equal,fixed_coefficient_carry_controls=relations,projection_cases_requiring_zero_trim=forced_padding,addition_controls=addition,scope='All-order nonregularity is the balanced-language proof; these are separate physical, arithmetic and normalization controls.')
print(json.dumps(result,indent=2,sort_keys=True))
