"""
opus-2026-07-11-S261: RUN the Beurling-Selberg mollification of the good set G' against the coprime core (the
S260 refined path). Outcome: the mollification is the RIGHT tool (finite degree captures the discrepancy) and
improves the bound ~17x, but the residual is the SIGNED CANCELLATION in Sum_h b_h ghat(-hv), not the magnitude.

SETUP. Covering v => core (coprime to 30030), non-core; G'={t:||wt||>=1/14 for non-core w}. Want coreCover<1.
Exact identity: density(D_v in G') = 1/7 + eps_v, eps_v = (Sum_{h!=0} b_h ghat(-hv))/|G'|, b_h=sin(pi h/7)/(pi h).
coreCover<1 <= Sum_core density < 1 (union) OR directly coreCover ~ 1-(6/7)^core (independent model, S260).

WHAT THE MOLLIFICATION SHOWS (FFT on grid D=13860):

(1) FINITE DEGREE SUFFICES. The truncated discrepancy Sum_{|h|<=K} b_h ghat(-hv) CONVERGES to eps_v fast:
    at K=50 it is essentially exact (tail beyond K=50 negligible). So a degree ~50 Beurling-Selberg majorant
    captures the full discrepancy -- the mollification is the correct, finite tool.

(2) TWO REGIMES. LARGE core runners (v>=17, coprime): eps_v SMALL (0.01-0.09) -- their frequencies hv>=17 hit
    the HIGH-frequency (small) part of ghat, so they equidistribute. RUNNER 1 (v=1): eps_1 LARGE (0.57 at the
    deep well) -- its frequencies h hit the LOW-frequency (large) part of ghat. Runner 1 is the exception;
    when it is the only core runner (deep well), coreCover = density(runner 1) < 1 is handled by S255 (near-AP).

(3) THE L2 BOUND improves the naive ~17x but still ignores cancellation. |eps_v| <= sqrt(tail_v)/(sqrt6*|G'|),
    tail_v = Sum_{|m|>=v}|ghat(m)|^2 (high-freq mass of G'). For large core v this is 0.4-0.9 (< 1 per arc,
    vs the naive N/(6v|G'|)~14). But the ACTUAL eps_v is 0.02 -- the L2 bound is ~40x too weak because it uses
    |ghat| magnitudes, discarding the SIGNED cancellation in Sum b_h ghat(-hv) that makes eps small.

THE RESIDUAL (sharpened). The true smallness of eps_v is a SIGNED CANCELLATION: Sum_h b_h ghat(-hv) is small
because, for v coprime to the non-core, the frequencies -hv are "generic" relative to the resonance structure
of ghat (the non-core lattice), so the signed sum cancels -- NOT because the terms are individually tiny.
Capturing this needs a BILINEAR / cancellation estimate on Sum_h b_h ghat(-hv) exploiting gcd(v, non-core)=1,
which the magnitude bounds (naive, L2) cannot see. This is the crux, now precisely named: a cancellation bound
for the coprime core against the non-core resonance lattice.

NET. Running the mollification: it is the right tool (finite degree K~50), improves the bound ~17x (L2 < 1 per
arc for large core), reduces to tail_v = high-freq mass of G', and isolates runner 1 (=> S255). But it does
NOT close the proof -- the residual is the SIGNED CANCELLATION in Sum_h b_h ghat(-hv), a bilinear estimate for
v coprime to the non-core lattice. The independent model coreCover ~ 1-(6/7)^core < 1 (margin (6/7)^core) is
the target; the cancellation is what must be proven to reach it.

-> opus-S260 (mollification path -- run here), opus-S259 (equidistribution mechanism), opus-S255 (runner-1/
near-AP), LRCFourierCompletion (the completion identity = a cancellation bound of this type), s558o.
"""
import numpy as np
from math import gcd, pi
from functools import reduce
import random
def divcomplete(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def primitive(v): return reduce(gcd,v)==1
def main():
    D=13860; c=1.0/14; cD=c*D
    random.seed(4); fams=[[1,2,3,4,5,6,7,8,9,10,11,12,182]]
    while len(fams)<4:
        v=sorted(random.sample(range(1,90),13))
        if primitive(v) and divcomplete(v): fams.append(v)
    mags=np.abs(np.fft.fftfreq(D,1/D)).astype(int)
    for v in fams:
        core=[x for x in v if gcd(x,30030)==1]; non=[x for x in v if gcd(x,30030)!=1]
        a=np.arange(D); safe=np.ones(D,bool)
        for w in non:
            r=(w*a)%D; safe &= (np.minimum(r,D-r)>=cD)
        g=safe.astype(float); gh=np.fft.fft(g)/D; Gm=g.mean(); power=np.abs(gh)**2
        print(f"core={core} |G'|={Gm:.3f}:")
        for vv in core[:6]:
            r=(vv*a)%D; Dv=(np.minimum(r,D-r)<cD).astype(float); eps=(Dv*g).sum()/g.sum()-1/7
            def tr(K):
                hs=np.arange(1,K+1); return (2*(np.sin(pi*hs/7)/(pi*hs)*np.real(gh[(-hs*vv)%D])).sum())/Gm
            tail=power[mags>=vv].sum(); Lb=(tail**0.5)/((6**0.5)*Gm)
            print(f"   v={vv:>4}: eps={eps:+.4f}, trunc K=50={tr(50):+.4f}, L2-bound={Lb:.3f} {'[runner-1]' if vv==1 else ''}")
    print("=> finite degree captures eps; L2 bound ~17x better than naive but ignores signed cancellation (the residual).")
if __name__=='__main__':
    main()
