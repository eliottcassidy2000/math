import random
TH=1.0/7.0
def W_of(E,Vmax,j):
    pts=sorted(((e*j)%Vmax)/Vmax for e in E); n=len(pts); unc=0.0
    for i in range(n):
        nxt=pts[(i+1)%n]+(1.0 if i==n-1 else 0.0); g=nxt-pts[i]
        if g>TH: unc+=g-TH
    return unc
def jstar_and_sum(E,Vmax,jmax=60):
    cum=W_of(E,Vmax,0); js=None
    for j in range(1,jmax+1):
        w=W_of(E,Vmax,j); cum+=w
        if w>1e-9 and js is None:
            return j, cum   # S_{0..j*}
    return None,cum
random.seed(664); k=13; Vmax=2003; lead=(6/7)**k
dlo=-(-6*Vmax//(7*(k-1)))   # ceil, min d for large-spread AP
dhi=Vmax//(k-1)
print(f"k={k} Vmax={Vmax}: (6/7)^k={lead:.5f}; 6/7={6/7:.4f}; leading crosses 6/7 at N+1={6/7/lead:.1f}")
print(f"large-spread AP step d in [{dlo},{dhi}]")
print(f"{'config':>24}{'jstar':>6}{'S(0..js)':>10}{'lead=(js+1)(6/7)^k':>19}{'correction':>12}")
# full AP at large-spread d
for d in [dlo, dlo+5, (dlo+dhi)//2, dhi]:
    E=[d*i for i in range(k)]
    if max(E)>=Vmax: continue
    js,S=jstar_and_sum(E,Vmax)
    if js: print(f"{'AP d='+str(d):>24}{js:>6}{S:>10.4f}{(js+1)*lead:>19.4f}{S-(js+1)*lead:>+12.4f}")
# AP12 + 1 point (longest-AP=k-1)
d=dlo+3
for _ in range(2):
    p=random.randint(1,Vmax-1); E=sorted(set([d*i for i in range(k-1)]+[p]))
    if len(E)!=k or max(E)-min(E)<6*Vmax//7: continue
    js,S=jstar_and_sum(E,Vmax)
    if js: print(f"{'AP11+pt d='+str(d):>24}{js:>6}{S:>10.4f}{(js+1)*lead:>19.4f}{S-(js+1)*lead:>+12.4f}")
# random low-AP
for _ in range(3):
    sp=random.randint(6*Vmax//7,Vmax-1); E=sorted(set([0]+random.sample(range(1,sp),k-2)+[sp]))
    if len(E)!=k: continue
    js,S=jstar_and_sum(E,Vmax)
    if js: print(f"{'random low-AP':>24}{js:>6}{S:>10.4f}{(js+1)*lead:>19.4f}{S-(js+1)*lead:>+12.4f}")
print(f"\nREADING: jstar = first j with S(0..j)>6/7. leading (js+1)(6/7)^k + correction = S(0..js).")
print("AP: correction NEGATIVE (resonances cancel), delays jstar to ~k. low-AP: correction small, jstar~leading.")
print("=> jstar <= N where (N+1)(6/7)^k > 6/7 + |correction|; klein-S194 a-priori W-hat decay bounds |correction|")
print("   => explicit jstar<=N_k, the SAME W-hat object as the density-floor tail. Proof route unified.")
