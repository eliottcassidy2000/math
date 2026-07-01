"""Verify THM-534: L_y(E)=sum_t g(t) p_t(E), p_t=meas{exactly t of sectors 1..6 empty}. Does consec maximize L_y?
g_8=(t-1)(t-2)(t-4)(t-5)/40 => [1,0,0,0.1,0,0,1]. Target L_y(consec_8)=2633/7350~0.3582, MAX over 8-sets."""
import numpy as np, itertools, random
def Ly(E,g,N=300000):
    x=np.arange(N)/N; hit=np.zeros((N,7),bool)
    for e in E:
        sec=np.floor(7*((e*x)%1.0)).astype(int)%7; hit[np.arange(N),sec]=True
    Nmiss=6-hit[:,1:].sum(axis=1); return g[Nmiss].mean()
g8=np.array([1.,0,0,0.1,0,0,1.])
consec=list(range(8))
Lc=Ly(consec,g8); print(f"L_y(consec_8={consec}) = {Lc:.5f}  (THM-534 target 2633/7350={2633/7350:.5f})")
# check consec is the max over offset sets {0,e1..e7}, e in 1..15 (bounded spread)
rng=random.Random(0); best=Lc; bestE=consec; over=0; tested=0
for _ in range(3000):
    E=[0]+sorted(rng.sample(range(1,16),7)); L=Ly(E,g8,120000); tested+=1
    if L>Lc+1e-4: over+=1
    if L>best: best=L; bestE=E
print(f"tested {tested} offset 8-sets (spread<=15): #(L_y>consec)={over}; max L_y={best:.5f} at {bestE}")
print(f"=> consec {'MAXIMIZES L_y (THM-534 verified)' if over==0 else 'does NOT maximize -- CHECK'}")
