import math,json,numpy as np
checks=0;rows=[]
def req(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
for D in range(5,151):
 if math.gcd(D,6)>1:continue
 i2=pow(2,-1,D);i3=pow(3,-1,D);i4=pow(4,-1,D)
 E=lambda x:x*i2%D;O=lambda x:(3*x+1)*i2%D
 Ei=lambda x:2*x%D;Oi=lambda x:(2*x-1)*i3%D
 for x in range(D):req(E(O(Ei(Oi(x))))==(x-i4)%D,'commutator')
 U=np.zeros((D,D));V=np.zeros((D,D))
 for x in range(D):U[E(x),x]=1;V[O(x),x]=1
 req(np.array_equal(U.sum(0),np.ones(D)) and np.array_equal(U.sum(1),np.ones(D)),'permutation E')
 req(np.array_equal(V.sum(0),np.ones(D)) and np.array_equal(V.sum(1),np.ones(D)),'permutation O')
 S=(U+U.T+V+V.T)/4;L=(np.eye(D)+S)/2;P=(U+V)/2
 vals=np.linalg.eigvalsh(L);gap=1-vals[-2];bound=math.sin(math.pi/D)**2/16
 req(gap+1e-10>=bound,'symmetric numerical bound')
 eig=list(np.linalg.eigvals(P));idx=min(range(len(eig)),key=lambda j:abs(eig[j]-1));eig.pop(idx);rho=max(abs(z) for z in eig)
 req(rho<=1-math.sin(math.pi/D)**2/8+1e-9,'forward numerical bound')
 rows.append({'D':D,'lazy_gap':float(gap),'proved_lower_bound':bound,'forward_spectral_radius':float(rho)})
print(json.dumps({'status':'PASS','checks':checks,'scope':'Permutation identities exact; eigenvalue comparisons numerical with tolerances1e-10 and1e-9; the general bounds are proved in the note','rows':rows},indent=2))
