from fractions import Fraction as F
import itertools
def covered(W,Vstar):
    R=F(1,2*Vstar); h=F(1,14); forb=[]
    for w in set(W):
        j=0
        while True:
            c=F(j,w); lo=c-h/F(w)
            if lo>R: break
            forb.append((max(F(0),lo),min(R,c+h/F(w)))); j+=1
    forb.sort(); merged=[]
    for a,b in forb:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    prev=F(0)
    for a,b in merged:
        if a>prev: return False  # found a gap
        prev=max(prev,b)
    return prev>=R  # fully covered up to R
# For V*=14..30, look for ANY W (size 2..4, all w>V*/7) fully covering window.
# allow w up to 4*Vstar (teeth must reach into (0,R]; large w give fine teeth).
for Vstar in range(14,31):
    wmin=Vstar//7+1
    found=None
    pool=list(range(wmin,4*Vstar+1))
    # size 2 first
    for r in (2,3):
        for W in itertools.combinations(pool,r):
            if all(7*w>Vstar for w in W) and covered(W,Vstar):
                found=(r,W); break
        if found: break
    print(f"V*={Vstar}: wmin={wmin} minimal covering W = {found}")
