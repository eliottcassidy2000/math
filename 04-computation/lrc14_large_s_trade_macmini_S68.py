#!/usr/bin/env python3
"""mac-mini-S68: close opus-S253's 'large-s trade' escape (single-killer covering-min).
M = M_core*v_f/(v_f+s); does a FAST binding runner s push M(S) < 14/183? Directly test M(S)
for single-killer covering 13-sets with fast binding runners. EARLY-EXIT: a config clears as
soon as ANY base q gives best_q >= 14/183 (then M(S) >= 14/183)."""
from fractions import Fraction as F
from math import gcd

target=F(14,183)
def clears(S, Qmax):
    """True if some q<=Qmax gives best_q(S) >= 14/183 (=> M(S)>=14/183). Also return best seen."""
    best=F(0)
    for q in range(2,Qmax+1):
        # only a coprime to q; but for min-distance, scan a in 1..q-1 (skip gcd>1)
        thr = target*q  # need min residue-dist >= thr
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            ok=True; mind=q
            for v in S:
                r=(a*v)%q; d=r if r<=q-r else q-r
                if d<mind: mind=d
                if F(d,q)<target: ok=False; break
            if ok:
                return True, F(mind,q)
            if F(mind,q)>best: best=F(mind,q)
    return False, best

def is_covering(S,n=14): return all(any(v%q==0 for v in S) for q in range(2,n+1))

print(f"target 14/183={float(target):.6f}. Testing single-killer covering configs w/ fast binding runners.")
print("A config CLEARS if some base q<=Qmax gives all residues >= 14q/183 from 0 (=> M>=14/183).\n")
below=[]; tested=0; mins=[]
def run(cores, killers, Qmax=260, label=""):
    global tested
    loc=[]
    for core in cores:
        for killer in killers:
            S=sorted(set(list(core)+[killer]))
            if len(S)!=13 or not is_covering(S): continue
            tested+=1
            cl,best=clears(S,Qmax)
            if not cl:
                # deeper check before declaring below-floor
                cl2,best2=clears(S,3*max(S)+5)
                if not cl2: below.append((S,best2))
                else: loc.append((best2,S))
            else: loc.append((best,S))
    return loc

# Family 1: {1..11, fs} + killer, fs = fast binding runner 13..70
c1=[list(range(1,12))+[fs] for fs in range(13,71)]
r1=run(c1,[182,364,182+13,182+14,546], label="F1")
# Family 2: drop a small runner, add a fast one (raise s), keep covering
c2=[]
for drop in range(1,12):
    for fs in range(13,90):
        c2.append([x for x in range(1,13) if x!=drop]+[fs])
r2=run(c2,[182,364,546],label="F2")
# Family 3: two fast runners in the core (push s high), minimal small part
c3=[]
for f1 in range(14,40):
    for f2 in range(f1+1,60):
        c3.append([1,2,3,4,5,6,7]+[f1,f2]+[11,13])  # keep some smalls for covering
r3=run(c3,[182,364],Qmax=220,label="F3")

print(f"tested {tested} single-killer covering configs.")
if below:
    print(f"*** {len(below)} configs with M < 14/183 (would LOWER the floor!):")
    for S,b in sorted(below,key=lambda x:x[1])[:12]:
        print(f"    best_q(M)~{float(b):.6f}  S={S}")
else:
    print("NONE below 14/183: every fast-binding-runner single-killer covering config CLEARS.")
    print("=> the large-s trade does NOT beat the deep well (opus-S253's open escape: evidence closed).")

allmins=sorted(r1+r2+r3)[:8]
print("\nClosest approaches to the floor (smallest cleared M) -- how close does large-s get?:")
for b,S in allmins:
    tag=" <-- deep-well-like" if b==target else ""
    print(f"    M~{float(b):.6f}={b}  S={S}{tag}")
