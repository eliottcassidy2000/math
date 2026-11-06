import json, hashlib
from pathlib import Path
from functools import lru_cache
ROOT=Path(__file__).resolve().parents[1]
MP=ROOT/'dependencies/MENUS.json'
menu=json.loads(MP.read_text())['5']
def need(b,msg):
 if not b: raise RuntimeError(msg)
@lru_cache(None)
def g(n):
 if n==1: return 1
 if n%2==0:return g(n//2)
 h=(n+1)//2
 return g(h)+h-h%2
forms=[]
for e in menu:
 for tr in [False,True]:
  M=e['matrix'];M=list(map(list,zip(*M))) if tr else M
  forms.append((e['word'],tr,M))
coeff=[]
for w,tr,M in forms:
 if not w:continue
 l,c,d,b=M[0][0],M[0][2],M[2][0],M[2][2]
 e,f,h,i=-M[0][1],-M[1][0],-M[1][2],-M[2][1]
 j=M[1][1]
 need(l%2==1 and 1<=l<=5 and 0<=c<=7 and 0<=d<=7 and 0<=b<=8 and min(e,f,h,i)>=0 and -1<=j<=1,('coefficients',w,tr))
 need(all(M[a][3]==M[3][a]==0 for a in range(4)),('support',w,tr))
 need(l+3*c-e>0,('odd fixed source positivity',w,tr))
 coeff.append(dict(word=w,transpose=tr,lambda_=l,c=c,d=d,b=b,e=e,f=f,h=h,i=i,j=j,odd_coefficient_lower_bound=l+3*c-e))
# Source tests include odd-core lemma and Fermat-shaped sources at large heights.
for n in range(1,65537):
 a=n//(n&-n)
 need(g(n)==g(a) and 1<=g(n)<=n,('source',n))
 if a>=3:need(2*g(a)>=a+3,('odd bound',a))
for j in range(1,101):need(g((1<<j)+1)==(1<<j)+1,('Fermat shaped',j))
# Complete rectangular domain, plus sparse dyadic towers, actual coordinates only.
sources=sorted(set(range(2,513))|{a<<r for a in range(1,32,2) for r in range(1,21)})
hits=[];zeros={'even_even':0,'even_odd':0};pairs=0
for u in sources:
 A=g(u);p=u%2;r=(u&-u).bit_length()-1
 for w,tr,M in forms:
  row=[sum(t*M[ii][j] for ii,t in enumerate((A,p,u,1))) for j in range(4)]
  for v in sources:
   B=g(v);q=v%2;s=(v&-v).bit_length()-1;G=max(A,B)
   n=row[0]*B+row[1]*q+row[2]*v+row[3];pairs+=1
   need(2*n>=u*v,('physical bound',w,u,v))
   if w=='' and r>=max(1,s):
    Z=A+B-v+2-q
    if Z==0:
     zeros['even_odd' if q else 'even_even']+=1
     need((n+1)&n!=0,('zero constant Mersenne',u,v))
   if (n+1)&n:continue
   need(n<=18630*G**8,('main bound',w,u,v))
   if r<s:continue
   need(n>=u and n+1>(1<<r),('valuation premise',w,u,v))
   if w:
    if r==0:need(n<=66*G**2,'odd odd')
    elif s==0:
     C=A*(M[0][0]*B+M[0][2]*v+M[0][1])
     need(0<C<=19*G**2 and (C+1)%(1<<r)==0 and (1<<r)<=20*G**2 and n<=940*G**4,'mixed')
    else:
     C=A*(M[0][0]*B+M[0][2]*v)
     need((M[0][0]*A*B+1)%(1<<s)==0 and (1<<s)<=6*G**2 and (C+1)%(1<<r)==0 and (1<<r)<=90*G**4 and n<=18629*G**8,'even even')
   elif r:
    Z=A+B-v+2-q
    need(Z!=0 and Z%(1<<r)==0,'empty divisibility')
    need(n<=(23*G**3 if q else 195*G**5),'empty bound')
   if len(hits)<20:hits.append([w,tr,u,v,n])
result={'status':'PASS','matrix_sha256':hashlib.sha256(MP.read_bytes()).hexdigest(),'sources':len(sources),'pair_forms':pairs,'zero_constant_nonmersenne_controls':zeros,'sample_mersenne_hits':hits,'coefficients':coeff,'scope':'Finite-exact controls; general theorem audited separately in AUDIT.md. Both orientations included, not a coverage proof.'}
print(json.dumps(result,indent=2))
