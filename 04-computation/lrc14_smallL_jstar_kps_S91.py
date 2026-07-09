# small-L (dissociated / well-distributed) => j* <= O(7)? Find the exact bound + mechanism.
# well-distributed E: no gap > Vmax/7 at j=1 (else j*=1 trivially). Among those, max j* by longest-AP.
import random
from collections import defaultdict
def maxgap_res(res,Vmax):
    pts=sorted(set(res)); n=len(pts); mg=0
    for i in range(n):
        nxt=pts[(i+1)%n]+(Vmax if i==n-1 else 0); g=nxt-pts[i]
        if g>mg: mg=g
    return mg
def jstar(E,Vmax,jmax=100):
    for j in range(1,jmax+1):
        if maxgap_res([(e*j)%Vmax for e in E],Vmax)*7>Vmax: return j
    return None
def longest_ap(E):
    E=sorted(set(E)); S=set(E); best=1
    for i in range(len(E)):
        for jj in range(i+1,len(E)):
            d=E[jj]-E[i]; L=2; nx=E[jj]+d
            while nx in S: L+=1; nx+=d
            pv=E[i]-d
            while pv in S: L+=1; pv-=d
            if L>best: best=L
    return best
random.seed(717)
print("small-L (well-distributed) => max j* by longest-AP; is it O(7) uniform in Vmax?")
for k in [11,13]:
    print(f"\n  k={k}:")
    for Vmax in [200,1001,5003]:
        byap=defaultdict(int)  # longest-AP -> max j*
        n=0
        for _ in range(4000):
            # well-distributed cluster: partition [0,Vmax) into k arcs each < Vmax/7
            E=sorted(set([0]+random.sample(range(1,Vmax),k-1)))
            if len(E)!=k: continue
            # require well-distributed (gap<=Vmax/7 at j=1, else j*=1)
            if maxgap_res(E,Vmax)*7>Vmax: continue
            js=jstar(E,Vmax); lap=longest_ap(E); n+=1
            if js: byap[lap]=max(byap[lap],js)
        mx=max(byap.values()) if byap else 0
        # split by whether longest-AP is small (<=k-3) or large
        small_max=max((v for l,v in byap.items() if l<=k-3), default=0)
        print(f"    Vmax={Vmax:5d}: {n} well-distr; MAX j*={mx}; longest-AP<=k-3 -> max j*={small_max}; by-AP: {dict(sorted(byap.items()))}")
print("\n=> if MAX j* for longest-AP<=k-3 stays <=~7 (Vmax-uniform), that's the O(7) small-L bound.")
print("   The high-longest-AP well-distributed configs (near-AP) are the O(k) embedded-AP case.")
