import json,math
C={}
def check(v,k):
 if not v:raise RuntimeError(k)
 C[k]=C.get(k,0)+1
def U(n):
 x=3*n+1;a=0
 while x%2==0:x//=2;a+=1
 return x,a
def P(n):
 a=2 if n%8==1 else 1;x=(3*n+1)//2**a;r=0
 while x%8==5:x=(x-1)//4;r+=1
 return x,a+2*r,r
def s(j):return(32*64**j-5)//9
def phi(n):return sum(math.gcd(n,k)==1 for k in range(1,n+1))
rows=[]
for j in range(31):
 n=s(j);check((32*64**j-5)%9==0 and n%8==3,'source')
 check(P(n)==(1,6*j+3,3*j+1),'primitive_clock')
 x,a=U(n);y,b=U(x);check((a,b,y)==(1,6*j+4,1),'odd_word')
 z=n;t=0;first=None
 while z!=1:
  z=(3*z+1)//2 if z%2 else z//2;t+=1
  if first is None and z<n:first=t
 check((t,first)==(6*j+5,4),'shortcut_clock_and_first_exit')
 z=n;t=0;first=None
 while z!=1:
  z=3*z+1 if z%2 else z//2;t+=1
  if first is None and z<n:first=t
 check((t,first)==(6*j+7,6),'ordinary_clock_and_first_exit')
 check((s(j+6)-n)%840==0 and P(s(j+6))[1]-P(n)[1]==36,'same_guard_full_clock')
 check(P(s(j+6))[1]%12==P(n)[1]%12,'same_modular_clock')
 rows.append({'j':j,'source':n,'P_clock':6*j+3,'shortcut_ROOT':6*j+5,'ordinary_ROOT':6*j+7,'first_shortcut_exit':4})
for M in [1,2,3,5,7,8,9,12,15,16,35,105,210,840,1024]:
 a=0;o=M
 while o%2==0:o//=2;a+=1
 L=phi(9*o);j=max(0,(a-5+5)//6)
 check(pow(64,L,9*o)==1,'period')
 check((s(j+L)-s(j))%M==0,'general_modulus')
 check(P(s(j+L))[1]-P(s(j))[1]==6*L,'general_clock_difference')
check(s(6)==244335917283,'printed_source')
print(json.dumps({'status':'PASS','checks':C,'total_checks':sum(C.values()),'rows':rows,'scope':'All-j result proved algebraically; finite replay is control, not global coverage.'},indent=2,sort_keys=True))
