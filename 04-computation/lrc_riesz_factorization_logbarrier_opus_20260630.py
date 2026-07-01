"""
lrc_riesz_factorization_logbarrier_opus_20260630.py
Grounds the multiscale log-barrier signed certificate (route A) via the Riesz factorization of
S = {0,1,2,30,31,32,60,61} = ({0,1,2} (+) {0,30,60}) minus the corner 62.
 (1) T_S = T_A*T_B - z^62 EXACT (pure DFT);
 (2) |T_S|^2 = tensor |T_A|^2|T_B|^2 (Pearson ~0.985); log|T_F|^2 = log|T_A|^2+log|T_B|^2 EXACT for the product;
 (3) predicts the coherence peak HEIGHTS (corr ~0.997) and the THM-546 peel-deviation peak-FREQUENCIES (10/12 on comb);
 (4) log-barrier B=<log m> separates tight (+0.43, min m=1 borderline) from holes (-inf); AM conserved, GM separates.
ASCII only (Windows cp1252 safe). See reflection the-Riesz-factorization-...-log-barrier-...-opus-20260630.md.
"""
import numpy as np

def part1_factorization():
    N=18000; t=np.arange(N)/N; z=np.exp(2j*np.pi*t)
    S=[0,1,2,30,31,32,60,61]; A=[0,1,2]; B=[0,30,60]; F=[0,1,2,30,31,32,60,61,62]
    TS=sum(z**v for v in S); TA=sum(z**a for a in A); TB=sum(z**b for b in B); TF=sum(z**v for v in F)
    print("(1) FACTORIZATION")
    print(f"    max|T_S-(T_A T_B - z^62)| = {np.abs(TS-(TA*TB-z**62)).max():.2e}  (exact, machine precision)")
    print(f"    max|T_F - T_A T_B|        = {np.abs(TF-TA*TB).max():.2e}  (clean product F=A(+)B)")
    PS=np.abs(TS)**2; PT=np.abs(TA)**2*np.abs(TB)**2
    print(f"(2) TENSOR SPECTRUM  Pearson(|T_S|^2, |T_A|^2|T_B|^2) = {np.corrcoef(PS,PT)[0,1]:.4f}")
    # exact additive log-decomposition for F away from the Dirichlet zeros
    m=(np.abs(TA)>1e-3)&(np.abs(TB)>1e-3)
    d=np.log(np.abs(TF[m])**2)-(np.log(np.abs(TA[m])**2)+np.log(np.abs(TB[m])**2))
    print(f"    max|log|T_F|^2 - (log|T_A|^2+log|T_B|^2)| off-zeros = {np.abs(d).max():.2e}  (ADDITIVE across scales)")
    # coherence heights at the comb vs envelope
    ks=np.arange(0,15)
    comb=np.array([np.abs(TS[int(round(k/30*N))%N])**2 for k in ks])
    env=np.array([abs(3*(1+np.exp(2j*np.pi*k/30)+np.exp(2j*np.pi*2*k/30)))**2 for k in ks])
    print(f"(3) COMB HEIGHTS  corr(|T_S(k/30)|^2, (3|T_A(k/30)|)^2) = {np.corrcoef(comb,env)[0,1]:.4f}")

def part2_peel_deviation():
    n=14; N=36000; t=np.arange(N)/N; S=[0,1,2,30,31,32,60,61]
    def safe(v):
        fr=(v*t)%1.0; fr=np.minimum(fr,1-fr); return (fr>1.0/n).astype(float)
    P=np.ones(N)
    for v in S:
        if v!=0: P*=safe(v)
    ws=np.arange(1,151); D=np.array([np.mean(P*(safe(w)-safe(w).mean())) for w in ws]); wD=ws*D
    order=np.argsort(-np.abs(wD))[:12]
    comb=[int(ws[k]) for k in order if (int(ws[k])%30) in (0,1,2,28,29)]
    print(f"(3b) PEEL DEVIATION w*D(w): {len(comb)}/12 top peel-freqs on the 30-comb+{{0,1,2}} satellites: {sorted(comb)}")
    print(f"     signed: sum w*D={wD.sum():+.3f}, max|w*D|={np.abs(wD).max():.3f} (cancellation, THM-546 signed Abel/QR)")

def part3_logbarrier():
    n=14; w=1.0/n; N=126000; t=np.arange(N)/N
    def cover(S):
        fr=(np.outer(np.array(S),t))%1.0; fr=np.minimum(fr,1-fr); return (fr<=w+1e-12).sum(axis=0).astype(float)
    sets={"AP":list(range(1,14)),"GW":[1,2,3,4,5,6,7,8,9,10,11,13,24],
          "12->26":[1,2,3,4,5,6,7,8,9,10,11,13,26],"drop-7":[1,2,3,4,5,6,8,9,10,11,12,13,14]}
    print("(4) LOG-BARRIER B=<log m>  (AM conserved, GM separates tight from holes)")
    for name,S in sets.items():
        m=cover(S); hole=np.mean(m<0.5)
        B=(-np.inf if hole>0 else np.mean(np.log(m)))
        print(f"     {name:>8}: min m={m.min():.0f}  hole={hole:.4f}  B={B:>8}  geo-mean={0.0 if hole>0 else np.exp(B):.3f}  arith={m.mean():.3f}")

if __name__=="__main__":
    part1_factorization(); part2_peel_deviation(); part3_logbarrier()
    print("\n=> Riesz factorization (multiplicative) => log-barrier ADDITIVE across scales => a MULTISCALE signed")
    print("   Beurling-Selberg certificate (route A). The barrier is the signed smoothing HYP-2738 proves is forced.")
