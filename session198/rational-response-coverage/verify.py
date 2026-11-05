"""Exact periodic-word and compressed physical recipe checks."""
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
import json,math,hashlib,sys
HERE=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def move(s,c):
 D,E,U,d=s
 return (D,0,2*U,d)if c=='O'else(D-E+U,d,2*U-d,d)

def apply(s,w):
 for c in w:s=move(s,c)
 return s

@lru_cache(None)
def g(L):
 if L==1:return 1
 if L%2==0:return g(L//2)
 h=(L+1)//2;return g(h)+h-h%2

def eval_bits(w):
 n=0;D=1
 for bit in w:
  if bit=='0':D+=n+n%2;n*=2
  else:n=2*n+1
 return n,D

def modular_word(M,u,z,extra=0):
 q=M;k=0
 while q%2==0:q//=2;k+=1
 a=2
 while pow(2,a,q)!=1%q:a+=1
 delta=((1<<a)-1)*((1<<k)-1)+1
 need(math.gcd(delta,M)==2,'joint compiler loop content')
 r=(u-1)%M;m=M//2
 j=0 if m==1 else (((z-g(r+1))//2)*pow(delta//2,-1,m))%m
 j+=(1+extra)*m
 return ('1'*a+'0'*k)*j+format(r,f'0{(M-1).bit_length()}b')

def egcd(a,b):
 if b==0:return a,1,0
 d,x,y=egcd(b,a%b);return d,y,x-(a//b)*y

def clock(a,b):
 need(a>0 and b>0 and a%2==0 and b%2==1 and math.gcd(a,b)==1 and b<2*a<2*b,'admissible slope')
 k=1
 while (2**k-1)%b:k+=1
 B=2**k;A=a*(B-1)//b
 eps=[(A>>i)&1 for i in range(k)];W=''.join('Z'if x else'O'for x in eps)
 R=sum(eps[i]*2**(k-1-i)for i in range(k))
 C=sum(eps[j]*eps[i]*2**(i-j-1)for i in range(k)for j in range(i))+sum(eps[i-1]*eps[i]for i in range(1,k))
 K=a*R-b*C
 return dict(a=a,b=b,k=k,B=B,A=A,R=R,C=C,K=K,word=W)

def construct(a,b,target):
 """Return a compressed actual source: V(L), then O, then W^repetitions.
 The interface dependency defines the literal V(L). No graph-size promise.
 """
 need(target%2==1,'odd target')
 c=clock(a,b);K=c['K'];need(K>0 and K%2==0,'positive even decrement')
 gcd,alpha,beta=egcd(2*a,b);need(gcd==1,'Bezout input')
 z=(target-2*a+b)//2
 D0=(1-2*beta*z)%K;U0=(1+2*alpha*z)%K
 need(D0%2==U0%2==1 and (2*a*U0-b*D0-target)%K==0,'odd seed congruence')
 u=(U0+1)//2;hres=(D0-u+u%2)%K;extra=0
 while True:
  seedword=modular_word(K,u,hres,extra)
  nn,gg=eval_bits(seedword);L=nn+1;D=gg+L-L%2;U=2*L-1
  need(L>=2 and D%K==D0 and U%K==U0,'physical modular seed')
  response=2*a*U-b*D
  if response>=target+K:break
  extra+=1
 repetitions=(response-target)//K
 need(repetitions>=1 and response-K*repetitions==target,'exact response correction')
 return {'a':a,'b':b,'target':target,'clock':c,'seed_L':L,'seed_bits':seedword,'seed_profile':(D,1,U,1),'prefix':'O','repetitions':repetitions,'final_s':1,'final_E':1,'final_d':1}

def literal_clock(W):
 T=[2,0];root=0
 for c in W:
  n=len(T);T.append(1 if c=='Z'else 0)
  for i in range(n):
   if not(c=='Z'and i==root):T[i]|=1<<n
 n=len(T);F=(1<<n)-1
 end=[[0]*n for _ in range(F+1)];start=[[0]*n for _ in range(F+1)];hs=[1]+[0]*F
 for S in range(1,F+1):
  for v in range(n):
   if not S>>v&1:continue
   Q=S^(1<<v)
   if Q==0:end[S][v]=start[S][v]=1
   else:
    end[S][v]=sum(end[Q][u]for u in range(n)if Q>>u&1 and T[u]>>v&1)
    start[S][v]=sum(start[Q][u]for u in range(n)if Q>>u&1 and T[v]>>u&1)
  hs[S]=sum(end[S])
 U=sum(hs[S]*start[F^S][root]for S in range(F+1))
 need(start[F][root]==1,'literal start interface')
 need(all(T[i]>>j&1 for i in range(1,n)for j in range(i+1,n)),'literal transitive deletion')
 return hs[F],end[F][root],U,hs[F^(1<<root)]

def main():
 for n in range(2,10001):
  ss=(g(n),n%2,n,1)
  need(move(ss,'O')==(g(2*n),0,2*n,1),'arithmetic orbit O')
  need(move(ss,'Z')==(g(2*n-1),1,2*n-1,1),'arithmetic orbit Z')
  if n%2:
   bad=3*n-4*g(n)
   need(bad<=0 or bad%4==1,'even-denominator hostile control')
 basis=[tuple(int(i==j)for i in range(4))for j in range(4)];clocks=[];recipes=[];literal=[];materialized=0
 for b in range(3,32,2):
  for a in range(2,b,2):
   if 2*a<=b or math.gcd(a,b)>1:continue
   c=clock(a,b);B,A,R,C,K,W=(c[x]for x in('B','A','R','C','K','word'))
   need(W[0]=='O'and W[-1]=='Z'and B//2<=A<B-1,'physical word boundaries')
   for s in basis:
    D,E,U,d=s
    need(apply(s,W)==(D+A*U-C*d,d,B*U-R*d,d),'word transfer basis')
    DD,EE,UU,dd=apply(s,W)
    need(a*UU-b*DD==a*U-b*D-K*d,'response basis')
   eps=[int(x=='Z')for x in W];u=Fraction(R,B-1);states=[];terms=[]
   for i,e in enumerate(eps):
    need(0<u<1 and e==(2*u).numerator//(2*u).denominator,'periodic digit orbit')
    states.append(u);terms.append(e*(u-eps[i-1]));u=2*u-e
   need(u==states[0]and sum(terms)==Fraction(K,b),'periodic increment')
   runs=[]
   for i,e in enumerate(eps):
    if e and not eps[i-1]:
     j=1
     while eps[(i+j)%len(eps)]:j+=1
     v=states[(i+j)%len(eps)];part=v+(1-v)/2**j
     need(part>0 and part==sum(terms[(i+h)%len(eps)]for h in range(j)),'positive one-run contribution')
     runs.append(part)
   need(sum(runs)==Fraction(K,b)and K>0 and K%2==0,'positive even clock')
   clocks.append(c)
   if len(W)<=6:
    actual=literal_clock(W);expected=apply((1,0,2,1),W)
    need(actual==expected and actual[1]==actual[3]==1,'literal clock response')
    literal.append({'a':a,'b':b,'word':W,'profile':actual})
   for target in range(-9,10,2):
    q=construct(a,b,target);D,E,U,d=q['seed_profile'];h=q['repetitions']
    need(D>=q['seed_L'] and U>=D and 2*a*U-b*D>0,'positive start from strong seed')
    if h*len(W)<=20000:
     # Evaluate the full enormous positive profile in these manageable cases,
     # independently of the signed-response jump used by construct().
     ss=(D,0,2*U,1)
     for _ in range(h):ss=apply(ss,W)
     need(ss[1]==ss[3]==1 and a*ss[2]-b*ss[0]==target,'materialized full response')
     materialized+=1
    recipes.append({'a':a,'b':b,'target':target,'seed_bits':len(q['seed_bits']),'seed_sha256':hashlib.sha256(q['seed_bits'].encode()).hexdigest(),'repetition_bits':h.bit_length(),'decrement':K})
 sys.path.insert(0,str(HERE/'dependencies'))
 import bilinear
 from core import inventories,cycles
 t=(1,0,2,1)
 for ss in basis:
  D,E,U,d=ss
  jq=(D-E+U,D-E,D+U,D)
  parent=bilinear.value('N',t,move(ss,'O'));child=bilinear.value('M',t,jq)
  need(child-2*parent==7*D-6*U-2*E+4*d,'C2 offset response basis')
 def graph(word):
  T=[2,0];root=0
  for c in word:
   n=len(T);T.append(0)
   for i in range(n):
    if c in('Z','J')and i==root:T[n]|=1<<i
    else:T[i]|=1<<n
   if c=='J':root=n
  return T,root
 c2=[]
 for word,nn,mm,delta in [('ZOZ',157,313,-1),('OZZ',187,375,1)]:
  T,r=graph('');V,c=graph(word+'O');W,d=graph(word+'J')
  pp=bilinear.construct('N',T,r,V,c)[0];cc=bilinear.construct('M',T,r,W,d)[0]
  measured=(inventories(pp)['H'],cycles(pp),inventories(cc)['H'],cycles(cc))
  need(measured==(nn,2,mm,2)and mm==2*nn+delta and (nn-19+2*delta)%34==0,'literal C2 offset pair')
  c2.append({'source_word':word,'parent_order':len(pp),'child_order':len(cc),'counts':measured,'offset':delta})
 out={'status':'PASS','C2_offset_controls':c2,'maximum_denominator':31,'arithmetic_orbit_targets':9999,'even_denominator_hostile_control':'All positive 3n-4g(n) tested for odd n<=10000 are1 modulo4; proof is all-order.','clocks':clocks,'compressed_exact_recipes':recipes,'materialized_positive_profiles':materialized,'literal_clocks':literal,'dependency_hashes':{str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted((HERE/'dependencies').glob('*'))if p.is_file()},'scope':'Exact all-odd signed response images for all reduced admissible rational slopes; no positive scalar coverage.'}
 text=json.dumps(out,sort_keys=True,indent=2)+'\n';(HERE/'VALIDATION.json').write_text(text);print(text,end='')

if __name__=='__main__':main()
