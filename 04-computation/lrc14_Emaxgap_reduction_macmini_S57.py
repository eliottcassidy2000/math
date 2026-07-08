"""mac-mini-S57 (THM-657): the E[maxgap] first-moment reduction for k=13 + stochastic-dominance negative result.
mu >= (7/6)(E[maxgap]-1/7); block minimizes E[maxgap]=0.2114 => mu >= 1.41x m_P.
NEGATIVE: block does NOT minimize mu_t at all t (perf-AP beats at t=0.30)."""
import numpy as np, random
random.seed(112358); TH=1/7; mP=0.0565
def stats(E,x,ts=None):
    Ea=np.array(sorted(E),float);ph=np.mod(np.outer(x,Ea),1.0);ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1);mg=g.max(axis=1)
    r=(float(mg.mean()), float((g*g).sum(axis=1).mean()))
    if ts is not None: r=r+(np.array([float((mg>t).mean()) for t in ts]),)
    return r
GRID=500_000;x=(np.arange(GRID)+0.5)/GRID;k=13;block=list(range(k))
Eb,Sb=stats(block,x)
print(f"E[maxgap](block)={Eb:.4f} => mu>=(7/6)(E[maxgap]-1/7)={7/6*(Eb-TH):.4f}={7/6*(Eb-TH)/mP:.2f}x m_P; target E[maxgap]>={TH+mP*6/7:.4f}")
gmin=Eb;below=0
for _ in range(200):
    E=sorted(random.sample(range(60),k))
    if len(set(E))<k:continue
    e=stats(E,x)[0]; gmin=min(gmin,e); below+= e<Eb-1e-3
print(f"200 random: min E[maxgap]={gmin:.4f}, #below block={below}; E[sum g^2](block)={Sb:.4f} (<target, too weak)")
ts=[0.10,1/7,0.20,0.30]; perf=[t for t in range(15) if t not in(4,9)][:k]
_,_,mb=stats(block,x,ts);_,_,mp2=stats(perf,x,ts)
print(f"stochastic dominance CHECK t={ts}: block {mb}, perf-AP {mp2} (perf < block at t=0.30 => NO dominance)")
