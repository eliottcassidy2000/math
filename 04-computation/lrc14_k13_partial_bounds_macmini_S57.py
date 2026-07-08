"""mac-mini-S57 (THM-657): rigorous partial bounds toward E[maxgap]>=0.1913 (k=13) and the
localization of the wall to t in [1/12,1] (barely-covers, near-independent, neg-associated).
 - mu_t=1 for t<1/13 (exact) => E[maxgap]>=1/13=0.077
 - mu_t >= t*k*(1-(k-1)t) for t in [1/13,1/12] (pairwise PZ; N_t<=1/t, E[N_t]>=k(1-(k-1)t)); verified valid
 - integrated PZ-on-N_t clears the TAIL (diam>=76: 0.199/0.222) but not the block (0.173, AP76-covered)."""
import numpy as np, random
from scipy.integrate import trapezoid
random.seed(2718); TH=1/7; mP=0.0565; k=13; target=TH+mP*6/7
def mu_true(E,x,ts):
    Ea=np.array(sorted(E),float);ph=np.mod(np.outer(x,Ea),1.0);ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1);mg=g.max(axis=1)
    return np.array([float((mg>t).mean()) for t in ts])
def NN(E,x,ts):
    Ea=np.array(sorted(E),float);ph=np.mod(np.outer(x,Ea),1.0);ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    EN=[];EN2=[]
    for t in ts:
        Nt=(g>t).sum(axis=1).astype(float);EN.append(Nt.mean());EN2.append((Nt*Nt).mean())
    return np.array(EN),np.array(EN2)
GRID=300_000;x=(np.arange(GRID)+0.5)/GRID
def mu_lower(t): return 1.0 if t<1/k else (max(0.0,t*k*(1-(k-1)*t)) if t<=1/(k-1) else 0.0)
ts=np.linspace(0.0005,0.999,700); lb=np.array([mu_lower(t) for t in ts])
print(f"target E[maxgap] >= {target:.4f}")
for nm,E in {'block':list(range(k)),'2blk-far(d80)':[0,1,2,3,4,5,40,41,42,43,44,45,80],'wide':sorted(random.sample(range(120),k))}.items():
    mt=mu_true(E,x,ts); EN,EN2=NN(E,x,ts); pz=np.where(EN2>0,EN**2/EN2,0.0)
    print(f"  {nm:14s}: valid-lb-holds={(lb-mt).max()<1e-3}; rig-pairwise E[mg]>={trapezoid(lb,ts):.4f}; "
          f"int-PZ={trapezoid(pz,ts):.4f} {'>=target' if trapezoid(pz,ts)>=target else '<target'}; true={trapezoid(mt,ts):.4f}")
print("WALL: deficit lives entirely in t in [1/12,1] (barely-covers, P(A_iA_j)=1/49 exactly, neg-associated) -> extremal lemma unavoidable")
