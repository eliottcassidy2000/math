import numpy as np, random
random.seed(3)
GRID = 120000
xs = (np.arange(GRID) + 0.5) / GRID
def mu_17(E):
    ph = np.mod(np.outer(xs, np.array(E, dtype=float)), 1.0); ph.sort(axis=1)
    gaps = np.diff(ph, axis=1)
    allg = np.concatenate([gaps, (ph[:,0]+1.0-ph[:,-1]).reshape(-1,1)], axis=1)
    return float((allg.max(axis=1) > 1.0/7.0).mean())
print("=== witnessG2 1/7-floor probe: mu_17(E)=meas{maxgap>1/7}, adversarially minimized (opus-S130) ===")
print(f"{'k':>3} {'consec':>8} {'AP2':>7} {'AP3':>7} {'rand-min':>9} {'ADV-min':>8}   adversarial E")
omin=1.0; worst=None
for k in range(8,14):
    cons=mu_17(list(range(1,k+1))); ap2=mu_17([2*i for i in range(1,k+1)]); ap3=mu_17([3*i for i in range(1,k+1)])
    rmin=1.0
    for _ in range(60):
        v=mu_17(random.sample(range(1,61),k)); rmin=min(rmin,v)
    amin=1.0; abest=None
    for _ in range(5):
        E=random.sample(range(1,50),k); cur=mu_17(E)
        for _ in range(25):
            i=random.randrange(k); cand=E.copy(); cand[i]=random.randint(1,90)
            if len(set(cand))<k: continue
            cv=mu_17(cand)
            if cv<cur: E,cur=cand,cv
        if cur<amin: amin,abest=cur,E
    km=min(cons,ap2,ap3,rmin,amin)
    if km<omin: omin=km; worst=(k,abest,amin)
    print(f"{k:>3} {cons:>8.4f} {ap2:>7.4f} {ap3:>7.4f} {rmin:>9.4f} {amin:>8.4f}   {sorted(abest)}", flush=True)
print(f"\nOVERALL MIN mu_17 (k=8..13, all tested/adversarial E): {omin:.4f}  at k={worst[0]}")
print(f"claimed rhoGlobFloorRat 0.336..0.443 ; m_P=0.0565")
print("=> mu_17 >> m_P and bounded well above 0 => 1/7 floor ROBUST (unlike refuted 2/7)." if omin>0.3 else "=> LOW mu_17 found -- floor AT RISK, investigate.")
