"""Growth analysis of the E-game exceptional sets (forward Z_2 and backward Z_3).

Inputs: the logs and bad-class dumps written by collatz_procgen_20260922_dim_forward.c (P0=30 run) and
collatz_procgen_20260922_dim_backward.c (P0=17 run).
Statistics:
  raw count N_m = |Bad_m|; local minima (levels right after a Beatty drop);
  projected count P_m(M) = |Bad_M mod 2^m| (classes mod 2^m that still have an exceptional descendant at
     level M) -- an UPPER bound for the covering number of Bad_inf at scale 2^-m, exact in the limit M->inf;
  Horton-Strahler number of the splitting tree of the projected set;
  fits on P_m, m in a window: exponential c*2^(d m), power law c*m^e, linear a+b m (RMS of log residuals),
  and local exponents d_loc = log2(P_{m+w}/P_m)/w, e_loc = log(P_{m+w}/P_m)/log((m+w)/m).
usage: python3 ..._dim_analysis.py FWD_LOG FWD_DUMPDIR BWD_LOG BWD_DUMPDIR [closure_file]
"""
import sys, math
import numpy as np
from collections import defaultdict
from fractions import Fraction as F

def read_log(fn, key):
    out={}
    for l in open(fn):
        if l.startswith(key+'='):
            f=dict(kv.split('=',1) for kv in l.split() if '=' in kv)
            out[int(f[key])]=f
    return out

def strahler(leaves, L, base):
    S={x:1 for x in leaves}
    for k in range(L-1,-1,-1):
        par=defaultdict(list)
        for x,s in S.items(): par[x % base**k].append(s)
        S={}
        for p,ch in par.items():
            ch.sort(reverse=True)
            if len(ch)==1: S[p]=ch[0]
            else:
                o=ch[0] if ch[0]>ch[1] else ch[0]+1
                S[p]=o
    return max(S.values())

def fits(ms, Ps, label):
    ms=np.array(ms,float); y=np.log(np.array(Ps,float))
    A=np.vstack([np.ones_like(ms),ms]).T; ce,res_e,_,_=np.linalg.lstsq(A,y,rcond=None)
    B=np.vstack([np.ones_like(ms),np.log(ms)]).T; cp,res_p,_,_=np.linalg.lstsq(B,y,rcond=None)
    C=np.vstack([np.ones_like(ms),ms]).T; cl,_,_,_=np.linalg.lstsq(C,np.exp(y),rcond=None)
    rms=lambda pred: float(np.sqrt(np.mean((y-pred)**2)))
    re=rms(A@ce); rp=rms(B@cp); lin=C@cl; rl=rms(np.log(np.maximum(lin,1e-9)))
    print(f"  {label}: window m=[{int(ms[0])},{int(ms[-1])}] ({len(ms)} pts)")
    print(f"    exponential: P ~ {math.exp(ce[0]):.1f} * 2^({ce[1]/math.log(2):.4f} m)   rms(log)={re:.4f}")
    print(f"    power law  : P ~ {math.exp(cp[0]):.2f} * m^{cp[1]:.3f}            rms(log)={rp:.4f}")
    print(f"    linear     : P ~ {cl[0]:.1f} + {cl[1]:.2f} m                rms(log)={rl:.4f}")
    return ce[1]/math.log(2), cp[1]

def local(ms, Ps, w, label):
    d=[]; e=[]
    for i in range(len(ms)):
        for k in range(i+1,len(ms)):
            if ms[k]-ms[i]==w:
                d.append((ms[i], math.log2(Ps[k]/Ps[i])/w, math.log(Ps[k]/Ps[i])/math.log(ms[k]/ms[i])))
    print(f"  {label} local exponents over windows of {w} levels (m: d_loc bits/level, e_loc):")
    print("    "+"  ".join(f"{m}:{a:.3f},{b:.2f}" for m,a,b in d))

def main():
    flog,fdir,blog,bdir=sys.argv[1:5]; clos=sys.argv[5] if len(sys.argv)>5 else None
    F_=read_log(flog,'m'); B_=read_log(blog,'r')
    Mf=max(int(m) for m in F_); Mb=max(int(r) for r in B_)
    print("=== FORWARD (Z_2): raw counts N_m ===")
    print("  "+" ".join(f"{m}:{F_[m]['bad']}" for m in sorted(F_)))
    mins=[m for m in sorted(F_) if 2<m<Mf and int(F_[m]['bad'])<int(F_[m-1]['bad']) and int(F_[m]['bad'])<=int(F_[m+1]['bad'])]
    print("  local minima:", ", ".join(f"m={m}:{F_[m]['bad']}" for m in mins))
    print("  smallest nonnegative representative of a non-trivial exceptional class (minrep) and min |signed rep| (minabs):")
    print("  "+"  ".join(f"{m}:{F_[m]['minrep']}/{F_[m]['minabs']}" for m in sorted(F_) if m>=30))
    bad={}
    # compare the top level with the last "clean" level (a local minimum at least 4 levels lower)
    Mclean=max([m for m in mins if m<=Mf-4] or [Mf-5])
    for M in (Mf, Mclean):
        try: bad[M]=[int(l.split()[0]) for l in open(f"{fdir}/fwd_bad_m{M}.txt")]
        except FileNotFoundError: pass
    Ms=sorted(bad)
    P={M:{m:len({x%(1<<m) for x in bad[M]}) for m in range(1,M+1)} for M in Ms}
    print(f"=== FORWARD projected counts P_m(M), M in {Ms} ===")
    print("  "+" ".join(f"{m}:{P[Ms[-1]][m]}" + (f"({P[Ms[0]][m]})" if m<=Ms[0] and P[Ms[0]][m]!=P[Ms[-1]][m] else "") for m in range(2,Ms[-1]+1,2)))
    conv=max(m for m in range(1,Ms[0]+1) if P[Ms[0]][m]==P[Ms[-1]][m])
    print(f"  P_m(M) identical for M={Ms[0]} and M={Ms[-1]} for all m <= {conv}")
    Mtop=Ms[-1]
    print("  Strahler number of the splitting tree:", ", ".join(f"L={L}:{strahler({x%(1<<L) for x in bad[Mtop]},L,2)}" for L in range(10,Mtop-4,6)))
    hi=conv if conv>=26 else Mtop-10
    if conv<26: print(f"  WARNING: short convergence window; fitting m<= {hi} anyway")
    ms=list(range(16,hi+1,2)); Ps=[P[Mtop][m] for m in ms]
    if len(ms)>=4:
        de,ep=fits(ms,Ps,"FORWARD P_m")
        local(ms,Ps,10,"FORWARD")
    if clos:
        pts=[F(l.split()[0]) for l in open(clos)]
        cl={m:len({(x.numerator*pow(x.denominator,-1,1<<m))%(1<<m) for x in pts}) for m in range(8,57,4)}
        print("  PROVED lower bound (classes meeting the Theorem-P closure):", " ".join(f"{m}:{c}" for m,c in cl.items()))
        fits(list(cl)[2:],list(cl.values())[2:],"FORWARD proved closure")
    print("=== BACKWARD (Z_3 units): raw counts N_r (classes mod 3^(r+1)) ===")
    print("  "+" ".join(f"{r}:{B_[r]['bad']}" for r in sorted(B_)))
    print("  minrep (smallest positive rep of a non-trivial exceptional class):")
    print("  "+"  ".join(f"{r}:{B_[r]['minrep']}" for r in sorted(B_) if r>=24))
    bb=[int(l.split()[0]) for l in open(f"{bdir}/bwd_bad_r{Mb}.txt")]
    bb2=[int(l.split()[0]) for l in open(f"{bdir}/bwd_bad_r{Mb-4}.txt")]
    Q={r:len({x%3**(r+1) for x in bb}) for r in range(1,Mb+1)}
    Q2={r:len({x%3**(r+1) for x in bb2}) for r in range(1,Mb-3)}
    print(f"=== BACKWARD projected counts Q_r(R), R={Mb} (R={Mb-4} in parentheses where different) ===")
    print("  "+" ".join(f"{r}:{Q[r]}"+(f"({Q2[r]})" if r in Q2 and Q2[r]!=Q[r] else "") for r in range(1,Mb+1)))
    convb=max([r for r in Q2 if Q2[r]==Q[r]] or [Mb-8])
    print(f"  Q_r identical for R={Mb-4},{Mb} for all r <= {convb}")
    print("  Strahler number:", ", ".join(f"L={L}:{strahler({x%3**(L+1) for x in bb},L+1,3)}" for L in range(6,Mb-2,4)))
    rs=list(range(8,max(convb,Mb-8)+1)); Qs=[Q[r] for r in rs]
    if len(rs)>=4:
        fits(rs,Qs,"BACKWARD Q_r (exponential rate in bits/level; divide by log2(3) for trits)")
        local(rs,Qs,8,"BACKWARD")
    # backward dyadic census: positive p/2^e (e<=40, value in [0.3,2]) in classes of Bad_Mb
    MODB=3**(Mb+1); pts=set()
    for e in range(0,41):
        q=2**e
        for c in bb:
            p=(c*q)%MODB
            if q*3//10<=p<=2*q and (p%2==1 or e==0): pts.add(F(p,q))
    print(f"=== BACKWARD dyadic census from Bad_{Mb}: {len(pts)} points p/2^e (e<=40) with value in [0.3,2];"
          f" range [{float(min(pts)):.4f},{float(max(pts)):.4f}]")
    from collections import Counter as _C
    print("  count by e:", dict(sorted(_C(x.denominator.bit_length()-1 for x in pts).items())))
    for L in (10,15,20,25):
        if L>Mb: continue
        proj={c%3**(L+1) for c in bb}
        cov={(x.numerator*pow(x.denominator,-1,3**(L+1)))%3**(L+1) for x in pts}
        print(f"  r={L}: projected classes {len(proj)}, containing a dyadic point: {len(proj&cov)}")
    # forward: census coverage lower bound and survival of |x|>3/2 candidates
    import glob, os
    cfiles=sorted(glob.glob(f"{fdir}/hostile_z13_M*_J*.txt"))
    if cfiles:
        cen=[F(l.strip()) for l in open(cfiles[-1]) if l.strip()]
        Ls=list(range(8,Mtop-3,2))
        print(f"=== FORWARD certified lower bound from {os.path.basename(cfiles[-1])}: classes mod 2^L containing a certified hostile point")
        print("  "+" ".join(f"{L}:{len({(x.numerator*pow(x.denominator,-1,1<<L))%(1<<L) for x in cen})}" for L in Ls))
        M0=int(os.path.basename(cfiles[-1]).split('_M')[1].split('_')[0]); J=int(os.path.basename(cfiles[-1]).split('_J')[1].split('.')[0])
        b0=[int(l.split()[0]) for l in open(f"{fdir}/fwd_bad_m{M0}.txt")]; MOD0=1<<M0; big=set()
        for j in range(0,J+1):
            q=3**j
            for x in b0:
                p=(-x*q)%MOD0
                if q*3//2<p<=q*8//5 and p%2==1 and (j==0 or p%3): big.add(F(-p,q))
        print(f"=== FORWARD candidates -p/3^j (j<={J}) with 3/2<|x|<=1.6 in Bad_{M0}: {len(big)}; still exceptional at higher levels:")
        for m in range(M0+1, Mf+1):
            fn=f"{fdir}/fwd_bad_m{m}.txt"
            if not os.path.exists(fn): continue
            S=set(int(l.split()[0]) for l in open(fn))
            print(f"  m={m}: {sum(1 for x in big if (x.numerator*pow(x.denominator,-1,1<<m))%(1<<m) in S)}")
main()
