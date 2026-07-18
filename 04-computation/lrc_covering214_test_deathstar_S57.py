# death-star-S57: DECISIVE test. Is there a covering-2..14 13-family with M<1/13 and NON-AP core?
# If NO -> boxeph THM-1017 (covering=2..14) stands; my cont22 (covering=2..13, mults-of-13) was the error.
# If YES -> THM-1017 refuted (court case).
# Search AP-core-adjacent: W = perturbation of {1..12} (covers 2..12, misses 13,14), v_max = mult of 182
# (forced to cover 13,14). Also W covering 13 or 14 handled separately below.
from fractions import Fraction as F
from math import gcd
from itertools import combinations
TH=F(1,13)
def M_lt_113(fam):   # returns M(V) if <1/13 else None (early bail once a witness gives >=1/13)
    Q=2*max(fam)+2; best=F(0)
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            r=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if 13*r>=q: return None       # >=1/13 witness -> not a counterexample
            if F(r,q)>best: best=F(r,q)
    return best
def covers(fam,lo,hi): return all(any(v%q==0 for v in fam) for q in range(lo,hi+1))
def is_dilated_AP(W):
    W=sorted(W); d=W[0]; return all(W[i]==d*(i+1) for i in range(12))
base=list(range(1,13))
found=[]; tested=0
pool=[x for x in range(1,35) if x%13 and x%14]   # replacement pool (miss 13,14)
# perturbations: replace r=1,2,3 elements of {1..12}
for r in range(0,4):
    for pos in combinations(range(12),r):
        repls=combinations([x for x in pool if x not in base],r) if r>0 else [()]
        for repl in repls:
            W=base[:]
            for idx,p in enumerate(pos): W[p]=repl[idx]
            if len(set(W))<12: continue
            W=sorted(W)
            if not covers(W,2,12): continue
            if any(v%13==0 or v%14==0 for v in W): continue   # W must miss 13,14 (AP-core case)
            # v_max = multiple of 182 (covers 13,14), > max(W)
            for k in range(1,12):
                vmax=182*k
                if vmax<=max(W): continue
                V=W+[vmax]
                if not covers(V,2,14): continue
                tested+=1
                M=M_lt_113(V)
                if M is not None:            # M(V) < 1/13
                    if not is_dilated_AP(W):
                        found.append((tuple(W),vmax,M))
print("covering-2..14 families (AP-core-adjacent, v_max mult of 182) with M<1/13 tested: %d"%tested)
print("NON-AP cores among them: %d"%len(found))
if found:
    for w in found[:20]: print("  *** NON-AP core with M<1/13, covering 2..14:",list(w[0]),"vmax=",w[1],"M=",w[2],"=%.5f"%float(w[2]))
    print("=> THM-1017 (covering 2..14) is REFUTED. Court case needed.")
else:
    print("=> ALL covering-2..14 M<1/13 families (in this search) have AP cores. THM-1017 stands for covering=2..14.")
    print("=> My cont22 'candidate=mults of 13' claim was the ERROR: it dropped the covering-14 requirement.")
# sanity: confirm the deep wells are found and are AP
dw=[base+[182]]
for V in dw:
    print("  deep well",V[:3],"...182: M<1/13?",M_lt_113(V) is not None,"core AP?",is_dilated_AP(V[:12]))
