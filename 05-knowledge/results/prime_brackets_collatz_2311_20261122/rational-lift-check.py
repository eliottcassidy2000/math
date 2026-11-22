from math import gcd

def val(a,p):
 if not a:return float("inf")
 j=0
 while not a%p:a//=p;j+=1
 return j

def data(w):
 n=k=Z=0
 for q in w:
  if q=="1":Z=3*Z+2**n;k+=1
  n+=1
 return n,k,Z

def rational(x):
 return [-x[0],-x[1]] if x[1]<0 else x

def actual(word):
 n,k,B=data(word);D=(1<<n)-3**k
 x,D=rational([B,D]);start=x
 for b in word:
  assert x%2==int(b)
  x=(3*x+D)//2 if b=="1" else x//2
 assert x==start
 return start,D

def parity_prefix(x,b):
 out=[]
 for i in range(b):
  v=x%2;out.append(str(v));x=(3*x+1)//2 if v else x//2
 return "".join(out)

def synth(sign,b,r2,s,r3,znum,zexp,u,r):
 A,B=("001010","001100") if sign==1 else ("0111","1011")
 n=len(A);k=A.count("1");d=max(0,r+u)
 prefix=parity_prefix(r2,b)
 eb=0
 while (z:=(-pow(2,-eb,3**s)))%3**s!=r3%3**s:eb+=1
 m0=1
 while m0*k<prefix.count("1") or m0*(n-k)<prefix.count("0"):m0+=1
 first=prefix+"1"*(m0*k-prefix.count("1"))+"0"*(m0*(n-k)-prefix.count("0"))
 end=1
 while end*(n-k)-eb<2 or end*k-s<=len(first)+s+2*n:end+=1
 final="0"+"1"*(end*k-s)+"0"*(end*(n-k)-eb-1)+"1"*s+"0"*eb
 pairs=[]
 for j in range(d):
  p,q,v=data(A);assert data(B)[:2]==(p,q)
  assert val(v-data(B)[2],11)==j
  pairs.extend([(A,B)]*10);A,B=A+B,B+A
 core="001010" if sign==1 else "0111"
 clocks=m0+end+10*(2**d-1);L=max(clocks+1,11**(r-1))
 L+=(-L)%11**(r-1)
 if val(L,11)>r-1:L+=11**(r-1)
 words=[a for a,_ in pairs]
 extra=core*(L-clocks);clockN=n*L;clockK=k*L
 clockD=(1<<clockN)-3**clockK
 assert val(clockD,11)==r
 denom=11**zexp;assert clockD%denom==0
 target=(clockD//denom)*znum
 for j in range(d):
  ok=False
  for digit in range(11):
   words[10*j:10*j+10]=[pairs[10*j][1]]*digit+[pairs[10*j][0]]*(10-digit)
   if (data(first+"".join(words)+extra+final)[2]-target)%11**(j+1)==0:ok=True;break
  assert ok
 word=first+"".join(words)+extra+final
 assert len(word)==clockN and word.count("1")==clockK
 runs=(word+"0").split("0")
 maxrun=end*k-s
 assert runs.count("1"*maxrun)==1 and len(max(runs,key=len))==maxrun
 assert word[0]=="0" or len(runs[0])+len(runs[-2])<maxrun
 U,V=actual(word)
 assert (U>0)==(sign>0) and gcd(V,6)==1
 assert (U-r2*V)%2**b==0 and (U-r3*V)%3**s==0
 assert val(U*denom-znum*V,11)-val(V*denom,11)>=u
 return clockN,val(V//gcd(U,V),11)

for sign in (1,-1):
 for q in [(5,13,3,10,7,0,1,1),(4,12,2,5,1,1,1,1),(5,0,2,7,23,2,0,2),(3,2,1,1,0,0,0,2)]:
  print(sign,synth(sign,*q))
for n in range(1,260):
 for k in range(1,n+1):
  v=val(2**n-3**k,11)
  assert (v>0)==(n%10==8*k%10)
  if v:
   assert (v>1)==(n%11==0)
   if v>1:assert (v>2)==((n//11-3*k)%11==0)
print("8 signed primitive simultaneous lifts; exhaustive short clock valuations pass")
