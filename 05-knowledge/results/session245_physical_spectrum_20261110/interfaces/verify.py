"""Self-contained exact literal reversal-bridge and cut-interface certificate."""
from itertools import permutations
from array import array
import json,sys

def need(x,s):
 if not x: raise RuntimeError(s)
def rev(a): return [((1<<len(a))-1)^x^(1<<i) for i,x in enumerate(a)]
def join(a,b):return [x|(((1<<len(b))-1)<<len(a)) for x in a]+[x<<len(a) for x in b]
def word(a,r,w):
 for c in w:
  n=len(a)
  if c=='O':a=join(a,[0])
  else:a=[x if i==r else x|(1<<n) for i,x in enumerate(a)]+[1<<r]
  if c=='J':r=n
 return a,r

def delete(a,r):
 ids=[i for i in range(len(a)) if i!=r]
 return [sum((a[i]>>j&1)<<k for k,j in enumerate(ids)) for i in ids],sum((a[r]>>j&1)<<k for k,j in enumerate(ids))
def tables(a):
 n=len(a);F=1<<n;incoming=[sum(1<<i for i in range(n) if a[i]>>j&1) for j in range(n)]
 e=[array('Q',[0])*F for _ in a];h=array('Q',[0])*F;h[0]=1
 for s in range(1,F):
  bits=s
  while bits:
   b=bits&-bits;bits-=b;j=b.bit_length()-1;t=s^b
   if not t:z=1
   else:
    z=0;pred=t&incoming[j]
    while pred:
     c=pred&-pred;pred-=c;z+=e[c.bit_length()-1][t]
   e[j][s]=z;h[s]+=z
 return h,e

def interface(a,c):
 h,e=tables(a);_,s=tables(rev(a));f=(1<<len(a))-1
 E=[sum(e[j][m] for j in range(len(a)) if c>>j&1) for m in range(f+1)]
 S=[sum(s[j][m] for j in range(len(a)) if c>>j&1) for m in range(f+1)]
 return (h[f],sum(h[m]*h[f^m] for m in range(f+1)),sum(h[m]*S[f^m] for m in range(f+1)),sum(E[m]*h[f^m] for m in range(f+1)),sum(E[m]*S[f^m] for m in range(f+1)),S[f],E[f])
def mul(x,y):
 h,q,A,G,T,S,E=x;j,r,B,F,U,V,W=y
 return h*j,q*r,A*r+h*B,F*q+j*G,A*F+T*j+h*U+E*V,S*j,h*W

def dual(x):h,q,A,G,T,S,E=x;return h,q,G,A,T,E,S

def closure(a,c):
 b=join(a,[0]);n=len(b)
 return [x|((1<<n) if not(c>>i&1) else 0) for i,x in enumerate(b)]+[c]

def strong(a):
 for r in range(len(a)):
  s=1<<r
  while True:
   t=s
   for i,x in enumerate(a):
    if s>>i&1:t|=x
   if t==s:break
   s=t
  if s!=(1<<len(a))-1:return False
 return True

def main():
 if '--fail' in sys.argv:need(False,'intentional failing sentinel')
 out=[]
 for i,w in enumerate(('JJJJ','JOZJ')):
  a,c=delete(*word([8,1,3,6],0,w));b,d=delete(*word([24,1,3,6,14],2,w));b=rev(b)
  za=interface(a,c);zb=interface(b,d);z=mul(mul(za,(1,4,0,0,0,0,0)),zb)
  g=join(join(a,[2,0]),b);cut=c|(d<<(len(a)+2));closed=closure(g,cut)
  need(len(closed)==19,'order19');need(strong(closed),'strong closure')
  H,e=tables(closed);count=H[-1];E=e[-1][-1];del H,e
  need(count==z[0]+2*z[2]-z[4],'independent literal scalar')
  need(E==z[0],'independent literal root ending')
  out.append(dict(word=w,interface=z,rows=closed,count=count,root_end=E,root_response=z[0]+2*z[2],deletion=z[0]))
 need([x['count'] for x in out]==[589901,589089],'advertised counts')
 need(out[0]['root_response']==out[1]['root_response']==603229,'same response')
 coprime=[]
 for w in ('JJJJ','JOZJ'):
  a,c=delete(*word([44,21,48,54,33,2],2,w));b,d=delete(*word([58,20,25,34,40,6],5,w));b=rev(b)
  z=mul(mul(interface(a,c),(1,4,0,0,0,0,0)),interface(b,d))
  g=join(join(a,[2,0]),b);cut=c|(d<<(len(a)+2));closed=closure(g,cut)
  need(strong(closed),'coprime strong physical graph')
  coprime.append(dict(word=w,interface=z,count=z[0]+2*z[2]-z[4],rows=closed,order=len(closed)))
 need([z['count'] for z in coprime]==[15669451,15654639],'coprime bridge counts')
 need(coprime[0]['interface'][2]==coprime[1]['interface'][2],'coprime cancellation')
 a=[18,8,3,5,14];small=[]
 for cut in (3,6):
  z=interface(a,cut);g=closure(a,cut);h,e=tables(g)
  need(h[-1]==z[0]+2*z[2]-z[4],'small closure')
  need(strong(g),'small strong')
  small.append(dict(cut=cut,interface=z,count=h[-1],rows=g))
 need([x['count'] for x in small]==[71,69],'small counts')
 print(json.dumps(dict(status='PASS',reversal_bridge=out,coprime_bridge=coprime,parent_cut_bridge=small,scope='Two physical pairs and proved symbolic laws; no anchored interval or universal coverage.'),indent=2))
if __name__=='__main__':main()
