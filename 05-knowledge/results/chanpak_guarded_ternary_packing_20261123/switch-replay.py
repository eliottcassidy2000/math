from collections import defaultdict
from itertools import combinations,permutations
from functools import cache
A8=(16,1,3,7,14)
A10=(20,1,2,7,14)
def sub(A,I):
 return tuple(sum((A[i]>>k&1)<<j for j,k in enumerate(I)) for i in I)
@cache
def h(A):
 n=len(A)
 if not n:return 1
 d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,1<<n):
  for i,c in enumerate(d[s]):
   if not c:continue
   q=A[i]&~s
   while q:
    b=q&-q;q-=b;d[s|b][b.bit_length()-1]+=c
 return sum(d[-1])
def signature(A,S):
 U=[0,0];V=[0,0];W=[0,0];K=[0,0]
 for k in range(1,len(A)+1):
  for P in permutations(range(len(A)),k):
   if not all(A[i]>>j&1 for i,j in zip(P,P[1:])):continue
   a,b=((S>>P[0])&1),1-((S>>P[-1])&1)
   w=h(sub(A,tuple(i for i in range(len(A)) if i not in P)))
   U[k%2]+=a*w;V[k%2]+=b*w;W[k%2]+=a*b*w;K[k%2]+=w
 return (W[0],U[1],U[0],V[1],V[0]),tuple(K)

def apex(A,S):
 n=len(A)
 return tuple(x+(0 if S>>i&1 else 1<<n) for i,x in enumerate(A))+(S,)

def tt(n):
 return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def join(A,B):
 n,m=len(A),len(B)
 return tuple(a|(((1<<m)-1)<<n) for a in A)+tuple(b<<n for b in B)

def cyclic(A,k):
 n=0
 for I in combinations(range(len(A)),k):
  for W in permutations(I[1:]):
   P=(I[0],)+W
   n+=all(A[i]>>j&1 for i,j in zip(P,P[1:]+P[:1]))
 return n

def autom(A):
 return sum(all(A[i]>>j&1==A[P[i]]>>P[j]&1
                for i in range(len(A)) for j in range(len(A))) for P in permutations(range(len(A))))

def probes(A,S):
 z=h(apex(A,S))
 l=[h(apex(join(tt(j),A),(S<<j)+(1<<j)-1)) for j in (1,2)]
 r=[h(apex(join(A,tt(j)),S)) for j in (1,2)]
 U=(r[0]-z)//2;V=(l[0]-z)//2
 return (z-9)//2,U,(r[1]-z)//2-2*U,V,(l[1]-z)//2-2*V



def flip(A,e):
 a,b=e;B=list(A);B[a]^=1<<b;B[b]^=1<<a;return tuple(B)
def five(A,S):
 D=h(A);(W,Uo,Ue,Vo,Ve),(Ke,Ko)=signature(A,S)
 Q=D+Ke+Ko
 return D,2*(Q+Uo+Ue),2*(Q+Vo+Ve),4*Q,D+Q+2*(W+Uo+Vo)
def extend(A,S,left,right):
 L=len(left); B=join(join(tt(L),A),tt(len(right)))
 bit=lambda w:sum((x=='O')<<i for i,x in enumerate(w))
 return apex(B,bit(left)+(S<<L)+(bit(right)<<(L+len(A))))
def tail(R,w):
 D,A,B,Q,H=R
 prev='I'
 for s in w:
  if s=='O': A,Q=2*A+D,2*Q
  else: A,B,Q,H=2*A,B+Q,2*Q,H+A+D*(prev=='O')
  prev=s
 return D,A,B,Q,H

if __name__=='__main__':
 A=(58,4,9,18,38,14);B=flip(A,(2,5));S=28
 u,v=five(A,S),five(B,S)
 assert u==(29,898,964,1280,751) and v==(29,898,964,1280,753)
 assert h(extend(A,S,'O','I'))==751 and h(extend(B,S,'O','I'))==753
 assert h(extend(A,S,'O','IOII'))==6255 and h(extend(B,S,'O','IOII'))==6257
 X=(82,36,9,19,38,77,30); Y=(250,84,25,50,32,70,156,62)
 assert h(X)==139 and h(Y)==251
 assert h(join(join((2,4,1),(2,4,9,3)),X))==2085
 assert h(extend(A,S,'O','IOIOI'))==9963 and h(extend(B,S,'O','IOIOI'))==9965
 assert h(extend(A,S,'O','IIOI'))==5299 and h(extend(B,S,'O','IIOI'))==5301
 Z={19:(48,1,3,39,46,6),31:(52,33,34,7,46,8),41:(48,33,35,7,14,24)}
 assert all(h(C)==n for n,C in Z.items())
 assert 3**4*h(Z[41])==3321 and 3*h(Z[19])*h(Z[31])==1767
 for m in range(36):
  for w,c,a in [('OIOI',-4695,14658),('IOI',-1943,7242)]:
   g0,g1=map(lambda x:tail(x,w+'I'*m),(u,v))
   assert g0[:4]==g1[:4] and g0[-1]==a*2**m+c and g1[-1]==g0[-1]+2
  r0,r1=map(lambda x:tail(x,'OII'+'I'*m),(u,v))
  assert r0[-1]==7300*2**m-1045 and r1[-1]==r0[-1]+2
  r0,r1=map(lambda x:tail(x,'I'*(2*m)),(u,v))
  assert r0[-1]==898*4**m-147 and r1[-1]==r0[-1]+2
 print('guard positive/negative seed, renewable rails, and quotient controls PASS')

if __name__=='__main__':
 A=(58,4,9,18,38,14); B=flip(A,(2,5)); S=28
 P=(214,68,72,51,38,71,152,62)
 assert h(extend(P,71,'O','I'))==12693
 for left,right,hs in [('OIOIO','IOI',(38077,38079)),('OIIOOIO','I',(38079,38081))]:
  assert tuple(h(extend(C,S,left,right)) for C in (A,B))==hs
 print('guarded parent 12693, two center hosts 38079, both signed exits pass')

if __name__=='__main__':
 assert five((2,0),3)==(1,14,8,16,9)
 assert five((2,4,1),7)==(3,42,24,48,27)
 assert h(extend((2,0),3,'O','I'))==9
 assert h(extend((2,4,1),7,'O','I'))==27
 print('literal uniform-contact neutral digit scales all five reads')
