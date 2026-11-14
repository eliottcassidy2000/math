from counter import *
import numpy as np
def build(M):
 T=np.zeros((4,36,36),dtype=np.int64)
 for d in range(4):
  for i,j in E[d]:T[d,i,j]+=1
 P={}
 for d in '02':
  powers=[np.eye(36,dtype=np.int64)]
  nxt=T[int(d)]
  seen={powers[0].tobytes():0}
  while True:
   a=powers[-1]@nxt%M;b=a.tobytes()
   if b in seen:
    P[d]=(powers,seen[b],len(powers)-seen[b]);break
   seen[b]=len(powers);powers.append(a)
 return T,P
def raw(T,P,w,M,e=1):
 v=np.zeros(36,dtype=np.int64);v[0]=v[18]=1
 for x in w:
  if isinstance(x,tuple):
   d,n=x;u,j,p=P[d];v=v@u[n if n<j else j+(n-j)%p]%M
  else:v=v@T[int(x)]%M
 return sum(v[18*r+6*a+2+c]+(v[18*r+6*a+c] if e!=r else 0) for r,a,c in product(range(2),range(2),range(2)))%M
def listres(M):
 T,P=build(M)
 print('params',M,[(d,j,p) for d,(u,j,p) in P.items()])
 def reps(d): u,t,p=P[d];return range(t+p)
 def high(d): u,t,p=P[d];return range(t,t+p)
 # All templates at run-level; extend editing inside constant subwords and each boundary.
 F=('1','1','0','2','3','1','1')
 def exp(F,k,l):return tuple(('0',k) if x=='0' else ('2',l) if x=='2' else x for x in F)
 R={}
 def rec(group,w,e=1): R.setdefault(group,set()).add(raw(T,P,w,M,e))
 for k,l in product(range(P['0'][1]+1,sum(P['0'][1:])+1),range(P['2'][1]+1,sum(P['2'][1:])+1)):
  w=exp(F,k,l);rec('hub flip',w,0)
  for i in [0,1,4,5,6]:
   rec('shell delete',w[:i]+w[i+1:])
   for x in '0123':rec('shell replace',w[:i]+(x,)+w[i+1:])
  for i in range(8):
   for x in '0123':rec('shell boundary insert',w[:i]+(x,)+w[i:])
  for i in range(1,7):
   a,b=F[i-1:i+1]
   w2=list(w[:i-1])
   if a in '02':w2.append((a,k-1 if a=='0' else l-1))
   w2.extend((b,a))
   if b in '02':w2.append((b,k-1 if b=='0' else l-1))
   w2.extend(w[i+1:])
   rec('swap',w2)
 for i,d in [(2,'0'),(3,'2')]:
  f='2' if d=='0' else '0'
  t=P[d][1]
  for a,b,l in product(reps(d),reps(d),high(f)):
   if max(a,b)<t:continue
   w=exp(F,l,l)
   for x in ('',*'0123'):
    mid=((d,a),)+tuple(x)+((d,b),)
    rec('run split edit '+d,w[:i]+mid+w[i+1:])
 return R
if __name__=='__main__':
 M=288
 T,P=build(M)
 assert [(j,p) for _,j,p in P.values()]==[(5,24)]*2
 R=listres(M)
 U=set.union(*R.values())
 assert all(x%2==1 for x in U)
 assert 3 not in U and len(U)==125
 print({s:len(rs) for s,rs in R.items()})
 print("joint obstruction: 3 is forbidden in",len(U),"possible residues mod 288")
 # Controls: neither of the two factors alone rules out the desired residue.
 w="11"+"0"*10+"2"*10+"311"
 h=H(w)
 assert h==1+36*2**10*3**10
 rep=lambda i,d:w[:i]+d+w[i+1:]
 d32=H(rep(14,"1"))-h
 d9=H(rep(0,"2"))-h
 assert (d32%32,d32%9)==(2,5) and (d9%32,d9%9)==(30,2)
 print("individual moduli fail:",d32,d9)
