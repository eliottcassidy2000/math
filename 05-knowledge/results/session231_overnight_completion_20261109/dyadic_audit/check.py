"""Independent audit. Imports no producer code. All checks survive python -O."""
from itertools import permutations
from fractions import Fraction
import json

def require(x, label):
 if not x: raise ValueError(label)

def valid(a):
 n=len(a)
 require(all(not(a[i]>>i&1) and a[i] < 1<<n for i in range(n)), 'loops/range')
 require(all(((a[i]>>j&1)+(a[j]>>i&1))==1 for i in range(n) for j in range(i)), 'tournament')

def path_tables(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    r=s^(1<<v)
    d[s][v]=sum(d[r][u] for u in range(n) if a[u]>>v&1) if r else 1
 return d

def count(a):
 valid(a);return sum(path_tables(a)[-1])

def state(a,cut):
 valid(a);n=len(a);full=(1<<n)-1;h=[];s=[];e=[]
 for mask in range(1<<n):
  good=[p for p in permutations([v for v in range(n) if mask>>v&1]) if all(a[u]>>v&1 for u,v in zip(p,p[1:]))]
  h.append(len(good));s.append(sum(bool(p) and bool(cut>>p[0]&1) for p in good));e.append(sum(bool(p) and bool(cut>>p[-1]&1) for p in good))
 require(h[-1]==count(a),'enumeration/dp')
 return [h[-1],sum(h[x]*h[full^x] for x in range(full+1)),sum(h[x]*s[full^x] for x in range(full+1)),sum(e[x]*h[full^x] for x in range(full+1)),sum(e[x]*s[full^x] for x in range(full+1)),s[-1],e[-1]]

def build(core,cut,word,tail=''):
 n=len(word)+len(core)+len(tail);a=[0]*n;marks=[]
 for i in range(n):
  for j in range(i+1,n):
   if len(word)<=i<j<len(word)+len(core):
    u,v=(i,j) if core[i-len(word)]>>(j-len(word))&1 else (j,i)
   else:u,v=i,j
   a[u]|=1<<v
 for i,b in enumerate(word):
  if b=='1':marks.append(i)
 marks.extend(len(word)+i for i in range(len(core)) if cut>>i&1)
 marks.extend(len(word)+len(core)+i for i,b in enumerate(tail) if b=='1')
 return a,sum(1<<i for i in marks)

def root(a,cut):
 n=len(a);return [z|(0 if cut>>i&1 else 1<<n) for i,z in enumerate(a)]+[cut]

def close(a,cut):
 n=len(a);return root([z|(1<<n) for z in a]+[0],cut)

def strong(a):
 for v in range(len(a)):
  reached={v};todo=[v]
  while todo:
   x=todo.pop()
   for y in range(len(a)):
    if a[x]>>y&1 and y not in reached:reached.add(y);todo.append(y)
  require(len(reached)==len(a),'strong')

def energy(word):
 return sum(2**(j-i-1)+(j==i+1) for i in range(len(word)) for j in range(i+1,len(word)) if word[i]==word[j]=='1')

def composition(x,y):
 h,q,a,g,t,s,e=x;H,Q,A,G,T,S,E=y
 return [h*H,q*Q,a*Q+h*A,G*q+H*g,a*G+t*H+h*T+e*S,s*H,h*E]

cores=[([4,65,26,99,43,135,181,31],189),([68,9,50,21,35,139,62,95],190)]
expected=[[323,5096,3672,3672,2610,252,254],[323,5096,3672,3672,2610,254,252]]
records=[]
for ix,(core,cut) in enumerate(cores):
 st=state(core,cut);require(st==expected[ix],'core state')
 a,marks=build(core,cut,'11001');F=count(root(a,marks));g=sum(z for i,z in enumerate(path_tables(a)[-1]) if marks>>i&1)
 q=count(close(a,marks))-F-g
 require([F,q,g]==[31957,131072,254-2*ix],'source counts')
 for b in (0,1):
  a,marks=build(core,cut,'11001','0'*b);graph=close(a,marks);strong(graph);H=count(graph)
  require(H==2**(18+b)-98861-2*ix and len(graph)==15+b,'initial graph')
  records.append({'core':ix,'vertices':len(graph),'count':H})
 # Arbitrary binary contexts, including the absent-prefix boundary.
 for pre,suf in [('', ''),('1',''),('0','1'),('1','1'),('10','01'),('11','10')]:
  a,marks=build(core,cut,pre,suf);P=int(pre or '0',2);R=int(suf or '0',2);b=len(suf)
  formula=(5096*P+3672)*(2**(b+1)-int(suf[::-1] or '0',2))+323*(2*R+1-energy(suf))-3672*P-323*energy(pre)-2610-(P%2)*(252+2*ix)-(254-2*ix)*int(suf[0] if suf else '0')
  require(count(close(a,marks))==formula,'arbitrary context formula')
require(next(j for j in range(1,638) if pow(2,j,637)==1)==84,'exact order')
require([j for j in range(3,185) if pow(2,j,5096)==3672]==[17,101],'all exponent residues')
M=2**84;B=(M-1)*459//637;w=bin(B)[2:];R=int(w[::-1],2);fb=energy(w)
require(len(w)==84 and w[0]==w[-1]=='1','block boundaries')
P=25;f=14
for t in range(9):
 m=17+84*t;word=bin(P)[2:];c=3672*P+323*f+2539
 require(5096*P+3672==2**m and P%2==1,'dyadic entry')
 require(len(word)==m-12 and energy(word)==f,'length/energy')
 require(2**(m-1)<c<c+2<2**m,'cost interval')
 prefix=[1,2,int(word[0]),int(word[0]),0,int(word[0]),int(word[0])]
 for z in word[1:]:
  z=int(z);prefix=composition(prefix,[1,2,z,z,0,z,z])
 for ix,st in enumerate(expected):
  total=composition(prefix,st)
  require(total[2]==2**m and total[4]-total[0]==c+2*ix,'all-word composition')
 P=M*P+B;f=f+R*((P-B)//M)+fb+1
alpha=Fraction(459,637)+Fraction(323*R,5096*(M-1))
require(Fraction(1,2)<alpha<Fraction(4,5),'limit coefficient')
print(json.dumps({'status':'PASS','literal_graphs':records,'order_637':84,'block':B,'reverse_block':R,'block_energy':fb,'limit_c_over_2m':str(alpha),'limit_decimal':float(alpha),'all_parameter_checks':'t=0..8; analytic proof in AUDIT.md','arbitrary_context_checks':12},sort_keys=True,indent=2))
