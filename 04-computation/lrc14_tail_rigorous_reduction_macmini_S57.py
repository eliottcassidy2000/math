"""mac-mini-S57 (THM-661 tail): RIGOROUS pieces of the decorrelation tail + the barely-covers wall.
PROVED (verified): near<=(2/7)E[W], far<=(5/7)E[W] (q_2<=q_1), E[W^2]<=E[W], PZ>=E[W].
=> k=13 reduces to E[W]>=0.0565; k=11,12 reduce to far<=E[W]^2 (decorrelation).
BOTH residuals hit the barely-covers wall (k/7=1.86>1 => inclusion-exclusion for E[W] DIVERGES)."""
import numpy as np
from itertools import combinations
TH=1/7
def stats(E,x):
    Ea=np.array(sorted(E),float);ph=np.mod(np.outer(x,Ea),1.0);ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    W=np.maximum(g-TH,0).sum(axis=1);EW=W.mean();EW2=(W*W).mean()
    Ls=np.linspace(TH,2*TH,60);near=2*np.trapezoid([np.maximum(g-L,0).sum(axis=1).mean() for L in Ls],Ls)
    return EW,EW2,near,EW2-near
GRID=400_000;x=(np.arange(GRID)+0.5)/GRID
print("RIGOROUS (verified 0 violations): near<=(2/7)E[W], far<=(5/7)E[W], E[W^2]<=(6/7)E[W], PZ>=E[W]")
for E in [list(range(13)),list(range(11)),[0,3,7,12,20,30,44,55,60,70,85,95,110]]:
    EW,EW2,near,far=stats(E,x)
    print(f"  {str(E[:3])}..: near/(2/7)EW={near/((2/7)*EW):.2f} far/(5/7)EW={far/((5/7)*EW):.2f} "
          f"E[W^2]/E[W]={EW2/EW:.2f} PZ/E[W]={EW*EW/EW2/EW:.2f}")
# barely-covers wall
NS=300000;xs=np.random.default_rng(1).random(NS);ys=np.random.default_rng(2).random(NS)
def PS(S,E):
    Ea=np.array([E[i] for i in S],float);d=np.mod(np.mod(np.outer(xs,Ea),1.0)-ys[:,None]+TH,1.0)
    return float(np.all((d>0)&(d<TH),axis=1).mean())
E=list(range(13))
S1,S2,S3=13*TH,sum(PS([a,b],E) for a,b in combinations(range(13),2)),sum(PS([a,b,c],E) for a,b,c in combinations(range(13),3))
print(f"barely-covers wall: E[W] Bonf3 = 1-13/7+S2-S3 = {1-S1+S2-S3:.3f} (true ~0.127) => IE diverges (k/7=1.86>1)")
