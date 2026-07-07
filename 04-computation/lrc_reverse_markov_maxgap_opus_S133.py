import numpy as np
GRID=500000; xs=(np.arange(GRID)+0.5)/GRID
def gaps_mat(E):
    ph=np.mod(np.outer(xs,np.array(E,float)),1.0); ph.sort(axis=1)
    return np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1]).reshape(-1,1)],axis=1)
def stats(E):
    g=gaps_mat(E); mg=g.max(axis=1)
    return float(mg.mean()), float(mg.max()), float((mg>1/7).mean())  # E[maxgap], B=max maxgap, mu_17

thr=1/7
print("=== REVERSE MARKOV reduction: mu_17 >= (E[maxgap]-1/7)/(B-1/7) (opus-S133, owner hint) ===")
print(f"    threshold 1/7={thr:.5f}; if E[maxgap] > 1/7 with margin => positive density floor\n")
print(f"{'family (k)':<26} {'E[maxgap]':>10} {'B=max':>7} {'mu_17':>7} {'revMarkov LB':>12} {'>1/7?':>6}")
def show(name,E):
    Em,B,mu=stats(E)
    lb=(Em-thr)/(B-thr) if B>thr else 0
    print(f"  {name:<24} {Em:>10.4f} {B:>7.4f} {mu:>7.4f} {lb:>12.4f}   {'YES' if Em>thr else 'no'}")
# AP (the binding minimizer) for k=8..13, and a few others
for k in range(8,14): show(f"AP {{1..{k}}}", list(range(1,k+1)))
print("  --- non-AP (should have E[maxgap] even larger) ---")
show("{2..14} consec", list(range(2,15)))
show("{1..4,10..18} sat", [1,2,3,4,10,11,12,13,14,15,16,17,18])
import random; random.seed(9)
mn=1.0
for _ in range(200):
    E=random.sample(range(1,60),13); Em,B,mu=stats(E); mn=min(mn,Em)
print(f"\n  min E[maxgap] over 200 random k=13 families: {mn:.4f}  (all > 1/7={thr:.4f}? {mn>thr})")
print(f"  AP is the binding (smallest E[maxgap]) case, like mu_17.")
print(f"\n  KEY: if inf_E E[maxgap] > 1/7 provable (mean max-gap of orbit), reverse-Markov => density floor.")
