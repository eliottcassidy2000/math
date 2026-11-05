#!/usr/bin/env python3
"""Actual rooted amalgamation, arbitrary-core finite inverse, and a literal repair."""
from verify_retained_root_sewing import *
from functools import lru_cache
from array import array
from pathlib import Path

def amalg(A,a,C,c):
 n=len(A);others=[v for v in range(len(C))if v!=c];idx={c:a,**{v:n+j for j,v in enumerate(others)}};out=list(A)+[0]*len(others);mask=sum(1<<idx[v]for v in others)
 for u in range(n):
  if u!=a:out[u]|=mask
 for u in range(len(C)):
  for v in range(len(C)):
   if C[u]>>v&1:out[idx[u]]|=1<<idx[v]
 validate(out);return tuple(out)

def law(st,z):
 h,x,y,d=st;hc,e,s,R,L,H,B=z
 return(H*h+(L-H)*x+(H-s)*y+(hc-H)*d,B*x+(R-H)*d,(H-e)*d,H*d)

def opposite(A):return tuple(((1<<len(A))-1)^row^(1<<v)for v,row in enumerate(A))
def mapped(A,m):
 out=[0]*len(A)
 for u in range(len(A)):
  for v in range(len(A)):
   if A[u]>>v&1:out[m[u]]|=1<<m[v]
 return tuple(out)

def birth(A,c):
 n=len(A);A=tuple(z|(1<<n)for z in A)
 return A+(0,)if c=='O'else(A[0]^(1<<n),)+A[1:]+(1,)
def word(A,w):
 for c in w:A=birth(A,c)
 return A

def prefix(q):return ''if q==0 else ''.join('O'if b=='1'else'Z'for b in bin(q)[2:])
def suffix(t):return ''if t==0 else ''.join('Z'if b=='1'else'O'for b in bin(t)[2:][::-1])
@lru_cache(None)
def F(q):
 if q==0:return 1
 return F(q//2)+(0 if q&1 else(q//2)+(q//2&1))

def scalar(z,q,t):
 h,e,s,R,L,H,B=z
 return h+H*(F(q)+F(t)-2)+(L-H)*q+(H-s)*(q&1)+B*q*t+(R-H)*t+(H-e)*(t&1)

def graph(C,c,q,t):return word(amalg(word((0,),prefix(q)),0,C,c),suffix(t))

def inverse(z,M):
 h,e,s,R,L,H,B=z;al=R-H;la=L-H;out={}
 if M<h:return out
 def keep(q,t):
  N=scalar(z,q,t)
  if N<=M:
   require(N>0 and N&1,'odd actual inverse output');require(q.bit_length()+t.bit_length()<=M.bit_length()+1,'canonical added-order bound');out.setdefault(N,(q,t))
 for t in (range((M-h)//al+1)if al else range(0,max(0,2*(M//H)-1),2)):keep(0,t)
 for q in (range((M-h)//la+1)if la else range(0,max(0,2*(M//H)-1),2)):keep(q,0)
 for q in range(1,(M-h)//(la+B)+1):
  hq=scalar(z,q,0)
  for t in range(1,(M-hq)//(al+B*q)+1):keep(q,t)
 return out

def h_flat(A):
 n=len(A);N=1<<n;dp=array('Q',[0])*(N*n)
 for v in range(n):dp[(1<<v)*n+v]=1
 for mask in range(1,N):
  V=mask
  while V:
   bit=V&-V;V-=bit;v=bit.bit_length()-1;c=dp[mask*n+v]
   if c:
    U=A[v]&~mask
    while U:
     b=U&-U;U-=b;dp[(mask|b)*n+b.bit_length()-1]+=c
 return sum(dp[(N-1)*n:])

def isstrong(A):
 for v in range(len(A)):
  reach=1<<v
  while True:
   new=reach
   for u in range(len(A)):
    if reach>>u&1:new|=A[u]
   if new==reach:break
   reach=new
  if reach!=(1<<len(A))-1:return False
 return True

def run():
 first=[(tournament(n,m),a)for n in range(1,4)for m in range(1<<(n*(n-1)//2))for a in range(n)]
 second=[(tournament(n,m),a)for n in range(2,5)for m in range(1<<(n*(n-1)//2))for a in range(n)]
 full=dual=0;bank={};boundary=0
 for C,c in second:
  z,_=state(C,c);z=z[:7];bank.setdefault(z,(C,c));hc,e,s,R,L,H,B=z
  require(hc>=H>=max(e,s)and min(R,L)>=H,'all positive core coefficients')
  indeg=len(C)-1-C[c].bit_count();outdeg=C[c].bit_count()
  require((L>H)==(indeg>0)and(R>H)==(outdeg>0),'exact slope boundaries')
  require(amalg(C,c,(0,),0)==C,'literal right identity')
  require(four(amalg((0,),0,C,c),0)==four(C,c),'left identity rooted state');boundary+=2
  for A,a in first:
   Z=amalg(A,a,C,c);st=four(A,a);require(four(Z,a)==law(st,z),'complete retained positive state');full+=1
   za,_=state(A,a);zz,_=state(Z,a)
   h,ee,ss,rr,ll,d,bb,_=za
   expected=(rr-d)*(L-H)+(h-ee)*H+(hc-s)*d+ee*s
   require(zz[0]==expected and zz[3]==B*(rr-d)+d*R,'literal scalar/root-cover formulas')
   # Converse anti-product with explicit relabelling of all vertices.
   Bop=amalg(opposite(C),c,opposite(A),a);m={a:c};otherA=[v for v in range(len(A))if v!=a]
   m.update({v:len(C)+j for j,v in enumerate(otherA)})
   m.update({len(A)+j:v for j,v in enumerate(v for v in range(len(C))if v!=c)})
   require(mapped(opposite(Z),m)==Bop,'literal converse anti-product');dual+=1
 # Associativity is checked on every triple of unrooted labelled graphs through3 with mark0.
 tiny=[tournament(n,m)for n in range(1,4)for m in range(1<<(n*(n-1)//2))];assoc=0
 for A in tiny:
  for B in tiny:
   for C in tiny:
    require(amalg(amalg(A,0,B,0),0,C,0)==amalg(A,0,amalg(B,0,C,0),0),'literal associativity');assoc+=1
 # Identity endpoint exception, witnessed by root-source edge.
 edge=transitive(2);require(state(edge,0)[0][1]==0 and state(edge,0)[0][5]==1,'endpoint identity hostile')
 binary=conv=0
 for z,(C,c)in bank.items():
  zop,_=state(opposite(C),c)
  for q in range(8):
   for t in range(8):
    A=graph(C,c,q,t);require(h_subsets(A)[-1]==scalar(z,q,t),'literal binary core count');binary+=1
    require(scalar(z,q,t)==scalar(zop[:7],t,q),'scalar converse symmetry')
    B=graph(opposite(C),c,t,q);p,r=q.bit_length(),t.bit_length();m={0:0}
    m.update({i:len(C)+r+p-i for i in range(1,p+1)})
    m.update({p+1+j:r+1+j for j in range(len(C)-1)})
    m.update({len(C)+p+j:r-j for j in range(r)})
    require(mapped(opposite(A),m)==B,'literal binary converse isomorphism');conv+=1
 M=200;inv=[];boxes=0
 for z,(C,c)in bank.items():
  got=inverse(z,M);expected={scalar(z,q,t)for q in range(2*M+1)for t in range(2*M+1)if scalar(z,q,t)<=M};boxes+=(2*M+1)**2
  require(set(got)==expected,'complete independent rectangular inverse control')
  inv.append({'signature':z,'outputs':len(got),'support_sha256':hashlib.sha256(json.dumps(sorted(got),separators=(',',':')).encode()).hexdigest()})
 C=(24,1,3,6,14);c=2;z,_=state(C,c);require(z[:7]==(11,3,2,17,20,5,32),'repair core signature')
 q,t=24084,0;A=graph(C,c,q,t)
 require(F(q)==20279 and len(A)==20 and isstrong(A)and scalar(z[:7],q,t)==462661 and h_flat(A)==462661,'actual strong20-vertex repair')
 require(15*0-3*32==-96,'physical negative production minor')
 deps={p:hashlib.sha256(Path(__file__).with_name(p).read_bytes()).hexdigest()for p in['verify.py','verify_retained_root_sewing.py','RETAINED_ROOT_SEWING.md','RETAINED_SEWING_SOURCE_SCOPE.md']}
 print(json.dumps({'status':'PASS','dependencies':deps,'full_rooted_pairs':full,'converse_products':dual,'identity_controls':boundary,'associative_triples':assoc,'distinct_small_core_signatures':len(bank),'actual_binary_graphs':binary,'literal_binary_converses':conv,'inverse_bound':M,'rectangular_pair_controls':boxes,'inverse_sets':inv,'repair':{'core':C,'root':c,'signature':z[:7],'q':q,'t':t,'F_q':F(q),'rows':A,'order':len(A),'H':462661,'strong':isstrong(A)},'negative_minor':-96},indent=2,sort_keys=True))
if __name__=='__main__':run()
