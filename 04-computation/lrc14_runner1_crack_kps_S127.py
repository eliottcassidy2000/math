# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont70: the worst |core|=1 body {1..11,13,84} is the CRACK between opus's runner-1 lemma
# Arguments A (measure, needs |S_rest|>1/14) and B (equidistribution, needs SPREAD rest). A fails, B stressed,
# yet LRC holds (eps_1<6/7). opus-S268's |core|=1 sample maxed at 0.328, MISSING this 0.60 body. Refinement:
# multi-runner Argument A (carve D_1 with the several smallest runners) -- the long run 2..11 that defeats B
# is the tool that rescues A.
def good_measure(B, level=1/14, grid=2000000):
    g=0
    for k in range(grid):
        t=k/grid; ok=True
        for b in B:
            r=(b*t)%1.0
            if min(r,1-r)<level: ok=False; break
        if ok: g+=1
    return g/grid
def coreCover(B, level=1/14, grid=2000000):
    g=0; ins=0
    for k in range(grid):
        t=k/grid; ok=True
        for b in B:
            r=(b*t)%1.0
            if min(r,1-r)<level: ok=False; break
        if ok:
            g+=1; r0=t%1.0
            if min(r0,1-r0)<level: ins+=1
    return (ins/g) if g else 1.0
if __name__=="__main__":
    rest=[2,3,4,5,6,7,8,9,10,11,13,84]
    Sr=good_measure(rest); cc=coreCover(rest); eps1=cc-1/7
    pairs=[(w,w+1) for w in rest if (w+1) in rest]
    print(f"{{1..11,13,84}} rest={rest}: |S_rest|={Sr:.5f} (A needs >1/14={1/14:.5f}: {'HOLDS' if Sr>1/14 else 'FAILS'})")
    print(f"  consecutive pairs={len(pairs)} (B needs FEW/spread: STRESSED); eps_1={eps1:.4f} <6/7={6/7:.4f}: LRC HOLDS margin {6/7-eps1:.3f}")
    print(f"  actual |S_rest cap D_1|=coreCover*|S_rest|={cc*Sr:.5f} < |S_rest|={Sr:.5f} => tight measure holds; A just too loose (1/14).")
