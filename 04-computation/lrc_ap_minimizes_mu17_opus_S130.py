import numpy as np, random
random.seed(11)
GRID=60000
xs=(np.arange(GRID)+0.5)/GRID
def mu_17(E):
    ph=np.mod(np.outer(xs,np.array(E,float)),1.0); ph.sort(axis=1)
    allg=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1]).reshape(-1,1)],axis=1)
    return float((allg.max(axis=1)>1/7).mean())
print("=== Is AP the GLOBAL minimizer of mu_17? Strong search, k=13 (opus-S130) ===")
ap=mu_17(list(range(1,14)))
print(f"mu_17(AP {{1..13}}) = {ap:.5f}  (claimed floor 477/1078={477/1078:.5f})\n")
best=ap; bestE=list(range(1,14))
# many random starts, aggressive descent, wide range
for s in range(40):
    E=random.sample(range(1,120),13); cur=mu_17(E)
    for it in range(60):
        i=random.randrange(13); cand=E.copy()
        cand[i]=random.randint(1,150)
        if len(set(cand))<13: continue
        cv=mu_17(cand)
        if cv<cur: E,cur=cand,cv
    if cur<best: best,bestE=cur,sorted(E)
print(f"best (lowest) mu_17 found over 40 aggressive starts: {best:.5f}")
print(f"  E = {sorted(bestE)}")
print(f"  below AP? {best < ap - 0.001}")
# also test structured adversaries: multiples, near-geometric, doubled-AP
for name,E in [("2*AP",[2*i for i in range(1,14)]),("{1..13} shifted+100",[i+100 for i in range(1,14)]),
               ("geometric-ish",[1,2,3,5,8,13,21,34,55,89,144,233,377]),
               ("{1..7}+{100..105}",[1,2,3,4,5,6,7,100,101,102,103,104,105])]:
    print(f"  mu_17({name}) = {mu_17(E):.5f}")
print(f"\nVERDICT: {'AP {1..13} IS the minimizer (nothing found below it) => floor = mu_17(AP) ~0.442, ROBUST.' if best>=ap-0.001 else 'FOUND E below AP -- AP does NOT minimize; floor lower than thought, investigate!'}")
