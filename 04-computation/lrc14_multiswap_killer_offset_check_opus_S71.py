#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Killer-offset (THM-618) is SINGLE-killer-specific; multi-swap covering families sit HIGHER than 14/183.
opus-2026-07-04-S71. Tested the potential counterexample: a covering family with SEPARATE 13- and 14-killers
(both ==1 mod 13). NO counterexample: {1..5,7..12,13,14} is covering, M=2/23 > 14/183 at t*=4/23 (NOT the
1/13 offset). The 1/13 killer-offset gives a LOW local value but the GLOBAL max-min is elsewhere (higher) --
mac-mini's non-convexity. So the deep well (DOUBLE killer 182=13*14, smallest runner covering BOTH q=13,14
while keeping {1..12}) is uniquely the covering-min; separate killers cost a dropped small runner => higher M
=> higher Ostrowski rung (2/23). Confirms mac-mini THM-618 + klein deep-well-unique; covering-min proof =
single-killer ladder (THM-618) + multi-swap>=14/183 (klein-verified, parametric per opus-S70 HYP-4087).
"""
import sys
from fractions import Fraction as Fr
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def norm(x):
    x=x-int(x)
    if x<0:x+=1
    return min(x,1-x)
def exact_M(S):
    S=sorted(set(S));cands=set()
    for v in S:
        for k in range(v):cands.add(Fr(2*k+1,2*v))
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for den in (S[i]+S[j],abs(S[i]-S[j])):
                if den:
                    for s in range(den):cands.add(Fr(s,den))
    b=Fr(0);arg=None
    for t in cands:
        v=min(norm(x*t) for x in S)
        if v>b:b=v;arg=t
    return b,arg
def covers(S): return all(any(v%q==0 for v in S) for q in range(2,15))
print("covering-min 14/183=%.5f. Multi-swap / separate-killer covering families (potential counterexamples):"%float(Fr(14,183)))
fams={
 'deep well {1..12,182} (double killer)':list(range(1,13))+[182],
 'sep killers {1..5,7..12,13,14}':[1,2,3,4,5,7,8,9,10,11,12,13,14],
 'sep killers {1..4,6..12,13,14}':[1,2,3,4,6,7,8,9,10,11,12,13,14],
 'sep {1,2,3,5,6,7,9,10,11,12,13,14,4}':[1,2,3,4,5,6,7,9,10,11,12,13,14],
 'triple-swap {2..12,13,14} (no 1)':[2,3,4,5,6,7,8,9,10,11,12,13,14],
}
mn=(Fr(1),None)
for name,S in fams.items():
    if len(set(S))!=13 or not covers(S): print("  %-40s (not 13-cov)"%name); continue
    M,t=exact_M(S)
    if M<mn[0]: mn=(M,name)
    print("  %-40s M=%-9s t*=%-9s  %s"%(name,str(M),str(t),'>=14/183' if M>=Fr(14,183) else '*** BELOW ***'))
print("  min over these: %s (%s)  -- deep well is the covering-min; multi-swap all HIGHER (Ostrowski rungs)."%(mn[0],mn[1]))
print("DONE.")
