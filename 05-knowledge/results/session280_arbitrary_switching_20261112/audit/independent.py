"""Independent exact audit calculations; no imports from candidate verifiers."""
from pathlib import Path
from fractions import Fraction as Q
import itertools,json,math,hashlib,random
ROOT=Path(__file__).resolve().parents[1]
def need(x,msg):
 if not x: raise RuntimeError(msg)
def g(n):
 if n==0:return 1
 q,r=divmod(n,2)
 return g(q)+(q+q%2 if not r else 0)
def matrix(word):
 p=int(word,2);q=int(''.join('1' if c=='0' else '0' for c in word[::-1]),2)
 return ((g(p),q+1,-1),(p+1,2**len(word),0),(1,0,0))
def mv(M,x):return tuple(sum(M[i][j]*x[j] for j in range(3)) for i in range(3))
def mm(M,N):return tuple(tuple(sum(M[i][k]*N[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def det(M):return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
def minors(M):
 return [M[i][a]*M[j][b]-M[i][b]*M[j][a] for i,j in itertools.combinations(range(3),2) for a,b in itertools.combinations(range(3),2)]
def thin(x):
 a,w,d=x;return a>0 and a%2==1 and 0<=d<=a<=w<2*a and (w-d)%2==0
ws=['1'+''.join(z)+'0' for n in range(2,9) for z in itertools.product('01',repeat=n-2)]
ms=list(map(matrix,ws));need(len(ws)==127,'count')
# Exhaustive collision check across all small formal source packets AND all127 words.
seen={};count=0
for A in range(1,16,2):
 for W in range(A,2*A):
  for d in range(W%2,A+1,2):
   x=(A,W,d)
   for word,M in zip(ws,ms):
    y=mv(M,x);need(thin(y),'invariant strip')
    if y in seen:need(seen[y]==(word,x),'disjoint image collision')
    seen[y]=(word,x);count+=1
    t=2**len(word);need(3*t*y[2]<2*y[1]<6*t*y[2],'dyadic separator')
    a=(y[1]*pow(y[2],-1,t))%t;need(a==M[1][0],'unknown word residue')
# Independently erase and recover variable-length itineraries, including non-Pisot letters.
def reverse_one(y):
 N,V,D=y
 need(thin(y) and D>0 and D%2,'eligible image')
 n=2
 while V>=3*(2**n)*D:n+=1
 t=2**n
 need(3*t*D<2*V<6*t*D,'length interval')
 a=V*pow(D,-1,t)%t
 need(t//2<a<t and a%2,'active digit')
 word=bin(a-1)[2:];need(len(word)==n,'decoded length')
 M=matrix(word);W=(V-a*D)//t;x=(D,W,M[0][0]*D+M[0][1]*W-N)
 need(thin(x) and mv(M,x)==y,'decoded predecessor')
 return word,x
rng=random.Random(280901);itinerary_checks=0
for trial in range(300):
 source=(3,3,1);labels=[];x=source
 for step in range(rng.randrange(1,101)):
  n=rng.randrange(2,21);word='1'+''.join(rng.choice('01') for _ in range(n-2))+'0'
  labels.append(word);x=mv(matrix(word),x)
 rev=[]
 while x!=source:
  word,x=reverse_one(x);rev.append(word)
 need(rev[::-1]==labels,'unknown complete itinerary');itinerary_checks+=1
# All three-count determinants, independently expanded.
for M in ms:
 for N in ms:
  P=mm(N,M);O=((1,0,0),M[0],P[0]);need(det(O)==N[0][1]*M[1][1],'observation determinant')
need(det(((1,0,0),(2,3,4),(2,2,3)))==1,'chart det')
# Sharp exterior/Smith depth for100 at every k<=40.
M=matrix('100');P=((1,0,0),(0,1,0),(0,0,1))
for k in range(1,41):
 P=mm(M,P);d1=math.gcd(*(x for row in P for x in row));d2=math.gcd(*minors(P));need(d1==1 and d2==2**(k-1) and det(P)==2**(3*k),'Smith sharpness')
# Independent exact pressure row reconstruction from stored certificate vectors.
data=json.loads((ROOT/'coverage/results.json').read_text());rays=((1,1,0),(1,1,1),(0,1,0));pressure={}
for key,low,threshold in [('subcritical',True,Q(973,1000)),('supercritical',False,Q(1003,1000))]:
 D=data[key];mats=list(map(matrix,D['alphabet']));c=D['vector'];rows=[]
 for i,M in enumerate(mats):
  ratios=[Q(mv(M,r)[0],mv(M,r)[1]) for r in rays];z=min(ratios) if low else max(ratios)
  value=sum(Q(c[j],1)/(N[1][1]+N[1][0]*z) for j,N in enumerate(mats))/c[i];rows.append(value)
 need(all(v<threshold if low else v>threshold for v in rows),'pressure '+key)
 pressure[key]={'extreme':str(max(rows) if low else min(rows)),'rows':len(rows)}
D=json.loads((ROOT/'coverage/low_oscillation_pressure.json').read_text());mats=list(map(matrix,D['alphabet']));c=D['vector'];rows=[]
for i,M in enumerate(mats):
 z=max(Q(mv(M,r)[0],mv(M,r)[1]) for r in rays)
 rows.append(sum(Q(c[j],1)/(N[1][1]+N[1][0]*z) for j,N in enumerate(mats))/c[i])
need(all(v>Q(101,100) for v in rows),'low oscillation pressure');pressure['low_oscillation']={'extreme':str(min(rows)),'rows':len(rows)}
# Independently expand every switched 2x2 quadratic inequality.
C=json.loads((ROOT/'switching/metric56.json').read_text());tau=Q(*C['tau']);D=C['Q_den']
metrics=[(Q(a,D),Q(b,D),Q(c,D)) for a,b,c in C['Q_num']];mats=list(map(matrix,C['words']))
need(len(mats)==56 and tau==Q(249,250),'switch alphabet/rate')
ener=Q(0);switched_tests=0
for i,M in enumerate(mats):
 qa,qb,qc=metrics[i]
 need(qa>Q(1,2) and qc>Q(1,2) and (qa-Q(1,2))*(qc-Q(1,2))>qb*qb,'metric lower')
 rr=[mv(M,r) for r in rays]
 native=[(a,Q(w+d,2),Q(w-d,2)) for a,w,d in rr]
 for j in range(3):need(max(v[j] for v in native)<=2*min(v[j] for v in native),'positive ray comparison')
 for x,y in itertools.permutations(rr,2):
  z1=y[1]-Q(y[0],x[0])*x[1];z2=y[2]-Q(y[0],x[0])*x[2]
  ener=max(ener,qa*z1*z1+2*qb*z1*z2+qc*z2*z2)
 for j,N in enumerate(mats):
  a,b,c=metrics[j];q=N[0][1];t=N[1][1]
  for r in rays:
   A,W,d=mv(N,mv(M,r));R=Q(W,A);S=Q(d,A)
   p=t-q*R;v=R;u=-q*S;w=S
   f00=a*p*p+2*b*p*u+c*u*u;f01=a*p*v+b*(p*w+u*v)+c*u*w;f11=a*v*v+2*b*v*w+c*w*w
   l00=tau*qa-f00;l01=tau*qb-f01;l11=tau*qc-f11
   need(l00>0 and l11>0 and l00*l11>l01*l01,'switched PSD');switched_tests+=1
initial=4*96**2*ener/Q(1,2)
need(initial*tau**3291<1 and initial*tau**3290>=1,'uniform certified threshold')
# Same independent expansion for the explicit151-letter extension.
C2=json.loads((ROOT/'switching/metric151.json').read_text());tau2=Q(*C2['tau']);D2=C2['Q_den']
metrics2=[(Q(a,D2),Q(b,D2),Q(c,D2)) for a,b,c in C2['Q_num']];mats2=list(map(matrix,C2['words']))
need(len(mats2)==len(set(C2['words']))==151 and len(metrics2)==151 and tau2==Q(9991,10000),'extension alphabet/rate')
need(all(2<=len(w)<=9 and w[0]=='1' and w[-1]=='0' and sum(a!=b for a,b in zip(w,w[1:]))<=3 for w in C2['words']),'extension domain')
ener2=Q(0);extension_tests=0
for i,M in enumerate(mats2):
 qa,qb,qc=metrics2[i]
 need(qa>Q(1,4) and qc>Q(1,4) and (qa-Q(1,4))*(qc-Q(1,4))>qb*qb,'extension metric lower')
 rr=[mv(M,r) for r in rays];native=[(a,Q(w+d,2),Q(w-d,2)) for a,w,d in rr]
 for j in range(3):need(max(v[j] for v in native)<=2*min(v[j] for v in native),'extension ray comparison')
 for x,y in itertools.permutations(rr,2):
  z1=y[1]-Q(y[0],x[0])*x[1];z2=y[2]-Q(y[0],x[0])*x[2]
  ener2=max(ener2,qa*z1*z1+2*qb*z1*z2+qc*z2*z2)
 for j,N in enumerate(mats2):
  a,b,c=metrics2[j];q=N[0][1];t=N[1][1]
  for r in rays:
   A,W,d=mv(N,mv(M,r));R=Q(W,A);S=Q(d,A)
   p=t-q*R;v=R;u=-q*S;w=S
   f00=a*p*p+2*b*p*u+c*u*u;f01=a*p*v+b*(p*w+u*v)+c*u*w;f11=a*v*v+2*b*v*w+c*w*w
   l00=tau2*qa-f00;l01=tau2*qb-f01;l11=tau2*qc-f11
   need(l00>0 and l11>0 and l00*l11>l01*l01,'extension PSD');extension_tests+=1
initial2=4*96**2*ener2/Q(1,4)
need(initial2*tau2**15416<1 and initial2*tau2**15415>=1,'extension threshold')
# Verify the pressure certificate on exactly the SAME151-letter alphabet.
D151=json.loads((ROOT/'coverage/synchronized_pressure151.json').read_text())
need(D151['alphabet']==C2['words'],'joint alphabet equality')
vec151=D151['vector'];pr151=[]
for i,M in enumerate(mats2):
 z=max(Q(mv(M,r)[0],mv(M,r)[1]) for r in rays)
 pr151.append(sum(Q(vec151[j],1)/(N[1][1]+N[1][0]*z) for j,N in enumerate(mats2))/vec151[i])
need(all(r>Q(1069,1000) for r in pr151),'joint151 pressure')
need(1024*1000**128<1069**128,'joint collision exponent')
need(512*1000**4096<1003**4096 and 512*100**1024<101**1024,'earlier collision exponents')
orderpoly={}
for w in C2['words']:orderpoly[len(w)+1]=orderpoly.get(len(w)+1,0)+1
need(orderpoly=={3:1,4:2,5:4,6:8,7:15,8:26,9:37,10:58},'joint packet order GF')
pressure['synchronized151']={'extreme':str(min(pr151)),'rows':151,'delta':'1/128'}
print(json.dumps({'status':'PASS','independent_strip_image_cases':count,'distinct_images':len(seen),'all_observation_pairs':len(ms)**2,'independent_unknown_itineraries':itinerary_checks,'sharp_Smith_depths':40,'pressure':pressure,'switching_PSD_checks':switched_tests,'switching_initial_energy':str(ener),'switching_uniform_depth':3292,'extension_PSD_checks':extension_tests,'extension_initial_energy':str(ener2),'extension_uniform_depth':15417,'scope':'Algebraic and rational checks; no unrestricted count coverage'},indent=2))
