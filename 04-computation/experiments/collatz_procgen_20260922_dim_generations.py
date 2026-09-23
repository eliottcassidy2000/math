"""Closure of {-1} under the perturbation lemma (Theorem P of collatz_procgen_20260922_exceptional_dimension.md).

Theorem P: let h = -p/3^{j_h} < 0 be hostile with all E-orbit values negative (H1), delta_h = |h|-1 < 1/2.
If i >= 1, j >= j_h, c >= 1 satisfy
   (B) j <= alpha_h(i)   [every E-path from h has made >= j x3-moves at its i-th halving]
   (A) (1 - 3^-j)/2 > delta_h + 2^i c / 3^j,
then x = h - 2^i c/3^j is hostile and satisfies (H1).  All points produced here are therefore PROVED
hostile (alpha_h(i) is computed exactly by collatz_procgen_20260922_dim_alpha.c).
usage: python3 ..._dim_generations.py ALPHA_BIN P0 IMAX   (prints generation statistics)
"""
import sys, subprocess
from fractions import Fraction as F
from collections import defaultdict, Counter

def alphas(binpath, P0, IMAX, hs):
    args=[binpath,str(P0),str(IMAX)]
    for h in hs: args += [str(h.numerator), str(h.denominator)]
    out=subprocess.run(args,capture_output=True,text=True,check=True).stdout.strip().splitlines()
    res={}
    for line,h in zip(out,hs):
        vals=[int(t) for t in line.split('alpha:')[1].split()]
        res[h]=vals          # vals[i-1] = alpha_h(i)
    return res

def v3(q):
    e=0
    while q%3==0: q//=3; e+=1
    return e

def children(h, scale_h, al, IMAX):
    """perturbations h - 2^i c/3^j allowed by Theorem P with i > scale_h (finer than h's own scale)."""
    dh=-h-1; jh=v3(h.denominator); out=[]
    for i in range(1, IMAX+1):
        a=al[i-1]
        for j in range(max(jh,1), a+1):
            # (A): (1-3^-j)/2 - dh > 2^i c/3^j  -> c < ((1-3^-j)/2 - dh) 3^j / 2^i
            room=(F(1,2)-F(1,2*3**j)-dh)
            if room<=0: continue
            cmax=room*3**j/2**i
            c=1
            while c < cmax:
                if c%2==1: out.append((h - F(2**i*c, 3**j), i, c, j))
                c+=1
    return out

if __name__=="__main__":
    binpath=sys.argv[1]; P0=int(sys.argv[2]); IMAX=int(sys.argv[3])
    gen={F(-1):(0,0)}          # point -> (generation, scale)
    frontier=[F(-1)]; g=0
    while frontier:
        al=alphas(binpath,P0,IMAX,frontier)
        new=[]
        for h in frontier:
            for x,i,c,j in children(h, gen[h][1], al[h], IMAX):
                if x not in gen:
                    gen[x]=(gen[h][0]+1, i); new.append(x)
        g+=1
        byscale=Counter(gen[x][1] for x in new)
        print(f"generation {g}: {len(new)} new PROVED hostile points; |x| range "
              + (f"[{float(-max(new)):.4f},{float(-min(new)):.4f}]" if new else "-")
              + f"; by scale i: {dict(sorted(byscale.items()))}", flush=True)
        frontier=new
    pts=sorted(gen, key=lambda x: -x)
    with open(sys.argv[4] if len(sys.argv)>4 else '/dev/null','w') as f:
        for x in pts: f.write(f"{x} {gen[x][0]} {gen[x][1]}\n")
    # covering counts: number of classes mod 2^L containing a generated point
    for L in range(8, IMAX+1, 4):
        cls={ (x.numerator*pow(x.denominator,-1,1<<L))%(1<<L) for x in gen }
        print(f"L={L}: classes mod 2^L meeting the generated set: {len(cls)}")
