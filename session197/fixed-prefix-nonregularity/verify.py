import json,re

def need(b,s):
 if not b:raise RuntimeError(s)
def g(L):
 need(L>=1,'positive gadget index');v=1
 while L>1:
  if L%2==0:L//=2
  else:
   L=(L+1)//2;v+=L-L%2
 return v
def f(n):return g(n+1)
def bits(n):return bin(n)[2:]
def main():
 recurrence=0
 for n in range(65536):
  need(f(2*n+1)==f(n),'odd f recurrence')
  need(f(2*n)==f(n)+n+n%2,'even f recurrence');recurrence+=2
 forms=filters=0;examples=[]
 for a in range(1,129):
  for b in range(1,129):
   n=((1<<a)-1)<<b;y=f(n);want=((1<<a)-1)*((1<<b)-1)+2
   need(y==want,'two run exact formula');forms+=1
   if min(a,b)<3:continue
   s,t=min(a,b),max(a,b)
   expected='1'*(a-1)+'0'*(a-1)+'11' if a==b else '1'*(s-1)+'0'+'1'*(t-s)+'0'*(s-2)+'11'
   need(bits(y)==expected and len(bits(n))==len(bits(y))==a+b,'canonical binary form and no padding')
   need(bool(re.fullmatch(r'1{2,}0{2,}11',bits(y)))==(a==b),'regular filter exactly balanced');filters+=1
   if (a,b) in ((3,3),(3,4),(4,3),(4,4),(7,11)):examples.append({'a':a,'b':b,'input':bits(n),'output':bits(y)})
 recovery=0
 for h,e,R in ((1,1,1),(3,1,3),(5,1,7),(13,2,18)):
  for L in range(1,4097):
   y=h*g(L)+R*L-e*(L%2)
   need(y>=0 and y+e*(L%2)-R*L==h*g(L),'fixed coefficient recovery')
   need((y+e*(L%2)-R*L)//h==g(L),'unique recovered gadget value');recovery+=1
 # Projection can retain zeros forced by a longer auxiliary track.
 padding=0;hostile=None
 for L in range(1,257):
  z=g(L);y=13*z+18*L-2*(L%2)
  width=max(L.bit_length(),z.bit_length(),y.bit_length())
  pair=list(zip(format(L,f'0{width}b'),format(z,f'0{width}b')))
  trim=0
  while len(pair)>1 and pair[0]==('0','0'):
   pair=pair[1:];trim+=1
  canonical_width=max(L.bit_length(),z.bit_length())
  need(pair==list(zip(format(L,f'0{canonical_width}b'),format(z,f'0{canonical_width}b'))),'left-quotient padding normalization')
  if trim and hostile is None:hostile={'L':L,'g':z,'auxiliary_y':y,'forced_zero_columns':trim}
  padding+=1
 need(hostile is not None,'a projected pair requiring trimming must be tested')
 # A concrete finite set of the all-order distinguishing suffixes in the proof.
 distinctions=0
 for i in range(3,65):
  for j in range(3,65):
   if i==j:continue
   good='1'*i+'0'*i;bad='1'*j+'0'*i
   need(good.count('1')==good.count('0') and bad.count('1')!=bad.count('0'),'distinguishing suffix');distinctions+=1
 print(json.dumps({'status':'PASS','recurrence_checks':recurrence,'two_run_formulas':forms,'binary_filter_checks':filters,'fixed_prefix_recoveries':recovery,'padding_controls':padding,'forced_padding_hostile':hostile,'distinguishing_suffix_controls':distinctions,'examples':examples,'scope':'Nonregular complete synchronous function graphs; no conclusion about regular or cofinite target ranges.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
