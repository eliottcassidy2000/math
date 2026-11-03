"""All-word finite-monoid certificate and mixed formal counterexample."""
from collections import Counter,deque

def need(ok,msg):
 if not ok:raise AssertionError(msg)
def mm(A,B,mod=None):
 C=tuple(tuple(sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))) for i in range(len(A)))
 return C if mod is None else tuple(tuple(v%mod for v in row) for row in C)
def vm(v,A):return mm((v,),A)[0]
def reduce(A,p):return tuple(tuple(x%p for x in row) for row in A)
def span(A,v,p):return any(tuple((i*A[0][k]+j*A[1][k])%p for k in range(2))==v for i in range(p) for j in range(p))
M={'a':((1,0,1),(2,0,1),(0,1,0)), 'b':((1,2,2),(4,3,4),(4,2,3)), 'q2':((5,3,4),(8,7,8),(6,6,7)), 'r21':((13,8,10),(20,14,17),(16,11,14)), 'r24':((5,16,16),(18,28,30),(18,22,24)), 'q4':((3,6,6),(6,8,9),(6,5,6)), 'r19':((11,7,10),(28,12,19),(26,11,16)), 'r32':((3,12,12),(12,22,23),(12,19,20)), 'r48':((5,16,16),(10,20,21),(10,13,14)), 'r18':((13,12,14),(28,19,24),(24,16,19))}
I3=((1,0,0),(0,1,0),(0,0,1));I2=((1,0),(0,1));t=(2,-1,0);B=((0,1),(0,2),(1,-1))
words={'L2':('r24','a','a','a'), 'L6':('r24','q4'), 'L45':('r19','r32','r48','q4'), 'L13':('r19','r19','a','a','r18'), 'b':('b',), 'q2':('q2',), 'Lp2':('r24','a','a','r21')}
def mat(word):
 P=I3
 for name in word:P=mm(P,M[name])
 return P
def quotient(P):
 Z=mm(P,B);Q=((Z[0][0]+Z[2][0],Z[0][1]+Z[2][1]),Z[0]);need(mm(B,Q)==Z,'quotient intertwining');return Q
P={name:mat(word) for name,word in words.items()};S={name:quotient(A) for name,A in P.items()}
expected={'L2':((240,339),(90,127)),'L6':((720,1101),(270,413)),'L45':((1020432,1560769),(366879,561148)),'L13':((211965,378780),(76220,136205)),'b':((5,8),(2,3)),'q2':((11,18),(4,7)),'Lp2':((4773,8628),(1789,3234))}
need(S==expected,'seven quotient matrices')
for name,A in P.items():
 z=vm(t,A);need(z[2]==0 and z[0]==-2*z[1],name+' is error-line return')
seen={I2};queue=deque([I2])
while queue:
 a=queue.popleft()
 for b in S.values():
  c=mm(a,b,30)
  if c not in seen:seen.add(c);queue.append(c)
need(all(mm(a,b,30) in seen for a in seen for b in S.values()),'closed under all generators')
zero3=[a for a in seen if reduce(a,3)==((0,0),(0,0))]
profile=Counter((span(reduce(a,2),(1,0),2),span(reduce(a,5),(0,1),5)) for a in zero3)
need((len(seen),len(zero3))==(6400,184),'finite closure counts')
need(profile=={(False,True):12,(True,False):111,(False,False):61},'joint gate profile')
print('mod30 reachable monoid:',len(seen),'zero mod3:',len(zero3))
print('gate profile: only mod5 12, only mod2 111, neither 61, both 0')
def norm_word(names):
 Q=I2
 for name in names:Q=mm(Q,S[name])
 need(all(x%3==0 for row in Q for x in row),'normalized word integral')
 return tuple(tuple(x//3 for x in row) for row in Q)
for names,expectedR,det_expected,gate2,gate5 in [
 (('Lp2','L2'),((640680,904601),(240140,339063)),-1300,False,True),
 (('Lp2','Lp2'),((12739007,23028132),(4774841,8631416)),16900,True,False)]:
 Q=norm_word(names);det=Q[0][0]*Q[1][1]-Q[0][1]*Q[1][0]
 need(Q==expectedR and det==det_expected and det%3 and span(Q,(1,0),2)==gate2 and span(Q,(0,1),5)==gate5,'separate gate examples')
 print('normalized example',names,'det',det,'passes (mod2,mod5)',gate2,gate5)
D=mat(words['L2']+words['b']+words['L2']);W=mm(M['b'],D);v=vm((31,14,25),D);w=(69,0,-21)
left=vm(w,W);right=tuple(3*x for x in v);delta=tuple(a-b for a,b in zip(left,right))
need(v==(33868343,15615527,26989950) and delta==tuple(216*x for x in t),'mixed descendant exact integer witness')
Rw=norm_word(('b','L2','b','L2'))
need(Rw==((1906550,2692023),(744750,1051577)) and Rw[0][0]*Rw[1][1]-Rw[0][1]*Rw[1][0]==100,'mixed quotient invertibility')
need(vm(w,B)==(-21,90) and vm(vm(w,B),Rw)==vm(v,B),'unique negative predecessor bulk')
print('mixed descendant witness: v=',v,'w=',w,'carry=216')
print('PASS')
