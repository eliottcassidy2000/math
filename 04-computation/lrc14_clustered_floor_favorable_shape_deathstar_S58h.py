"""The clustered floor M>=1/(2*rho) (rho=v_max/v_min) and the favorable-shape decomposition of the
AP-extraction kernel (death-star-S58h).
 - M(V) >= 1/(2*rho): at t=1/(2*v_max), every ||v t|| = v/(2 v_max) in [1/(2 rho), 1/2]. Verified sharp.
 - => M<1/13 => rho>6.5: fully-clustered families are NOT strict interior.
 - Covering-core gap (Lemma G, THM-1028): far-from-AP covering-2..12 12-cores have M(W)>=1/13+0.026 (~1/10);
   verified min 0.0265 over 10k. Clustered covering-2..13 families have M>=1/5 (0/2678 below 1/13).
 - Kernel splits: near-AP (Hamming<=6, THM-1004/5/6) + clustered (M>=1/(2rho)) DONE; spread far core = residual."""
from fractions import Fraction as F
from math import gcd
import random,time
def Mexact(V):
    mx=max(V);Qs=set()
    for i in range(len(V)):
        for j in range(i,len(V)):
            s=V[i]+V[j]
            for d in range(2,2*mx+1):
                if s%d==0:Qs.add(d)
    best=F(0)
    for q in Qs:
        for a in range(1,q):
            if gcd(a,q)!=1:continue
            m=min(min((v*a)%q,q-((v*a)%q)) for v in V)
            if F(m,q)>best:best=F(m,q)
    return best
random.seed(2);bad=0;n=0;t0=time.time();worst=(None,9.9)
while time.time()-t0<25:
    k=random.randint(3,13);V=sorted(random.sample(range(1,50),k))
    M=Mexact(V);rho=F(max(V),min(V));n+=1
    if M<1/(2*rho):bad+=1
    r=float(M*2*rho)
    if r<worst[1]:worst=(V,r)
print("M>=1/(2rho) verified on %d families: violations=%d; tightest M*2rho=%.3f at %s"%(n,bad,worst[1],worst[0]))
print("consequence: M<1/13 => rho>6.5 (fully-clustered families are not strict interior)")
