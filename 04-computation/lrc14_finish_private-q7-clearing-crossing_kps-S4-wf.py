# LRC(14) — private-q7-clearing-crossing angle (kps-S4-wf)
# ANGLE: make codex's CLEARING-CROSSING route (HYP-2578/2579) rigorous.
# By THM-524, M(S)=g(tau*) attained at a binding-pair crossing tau*=j/D, D=v_a +/- v_b.
# Non-tautological target: the parked multiple-of-14 runner w=14m, which is uniquely
# responsible for q=14 (private q-obligation), forces the binding crossing index j>=D/14,
# hence M=j/D>=1/14.
#
# We EXACTLY measure, over many covering S3 sets:
#  (Q1) Is tau* always of the form j/D for a binding pair (v_a,v_b), D=v_a+v_b or |v_a-v_b|?
#  (Q2) Is the parked runner (the multiple of 14) ITSELF in the binding pair? Or merely
#       a "clearing" constraint that sets the level?
#  (Q3) Does j>=D/14 always hold, and is it tight?
#  (Q4) The PRIVATE-q claim: when 14m is the unique multiple of 14 (so it privately covers q=14),
#       does removing it ever DROP M below 1/14 while a non-parked binding pair still governs?
#  (Q5) Mechanism: at tau*, which runners are "binding" (||v tau*||=M)? Is the parked runner
#       always binding, or is the level set by small runners?
#
# All decisions exact via fractions.Fraction.

from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r

def safe_components(A,h=F(1,14)):
    iv=[]
    for u in A:
        for j in range(0,u):
            c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
            if a<b: iv.append((a,b))
            else: iv.append((a,F(1))); iv.append((F(0),b))
    iv.sort(); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe

def Wwidth(A):
    sc=safe_components(A)
    if not sc: return F(0)
    ws=[b-a for a,b in sc]
    if sc[0][0]==0 and sc[-1][1]==1 and len(sc)>1: ws.append(sc[0][1]+(1-sc[-1][0]))
    return max(ws)

def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C

def Mval(S):
    b=F(0); bt=None
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v; bt=t
    return b,bt

def is_cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))

def is_prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1

def binding_set(S, t, M):
    """Runners v with ||v t|| == M (the binding runners at tau*=t)."""
    return [v for v in S if nrm(v*t)==M]

def classify_tau(S, t, M):
    """Return list of binding pairs (va,vb,kind,D,j) consistent with t=j/D being a crossing.
       A crossing is tau*=k/(va+vb) or k/(va-vb) with ||va t||=||vb t||=M.
       We find all unordered pairs of binding runners and report D, and j such that t=j/D.
    """
    B=binding_set(S,t,M)
    out=[]
    for va,vb in combinations(sorted(set(B)),2):
        for D in (va+vb, abs(va-vb)):
            if D==0: continue
            # t must equal j/D exactly for some integer j
            jt = t*D
            if jt.denominator==1:
                j=int(jt)
                out.append((va,vb,D,j))
    return B, out

def private_q_owner(S):
    """For each q in 2..14, the runners that are multiples of q. Return dict q->list.
       'private' for q means exactly one runner is a multiple of q.
    """
    d={}
    for q in range(2,15):
        d[q]=[v for v in S if v%q==0]
    return d

# --------- S3 covering-set generator ---------
# S3 = k>=2 (at least 2 speeds > 13) AND Vmax >= 13*Vmin.
# Covering: must contain a multiple of every q in 2..14, in particular a multiple of 14.
# We build sets: a small part P subset of {1..13} plus a large cluster L of speeds>13,
# ensuring covering. We force a unique multiple of 14 (the parked runner) when possible.

def gen_covering_S3(num=4000, seed=1, maxbig=400):
    random.seed(seed)
    out=[]
    tries=0
    while len(out)<num and tries<num*200:
        tries+=1
        # choose small part size, then large part to reach 13 total
        ksmall=random.randint(7,12)
        P=sorted(random.sample(range(1,14), ksmall))
        klarge=13-ksmall
        if klarge<2:   # need k>=2 large for S3
            continue
        L=sorted(random.sample(range(14,maxbig+1), klarge))
        S=sorted(set(P+L))
        if len(S)!=13: continue
        if not is_prim(S): continue
        if not is_cov(S): continue
        Vmin=min(S); Vmax=max(S)
        kbig=sum(1 for v in S if v>13)
        if kbig<2: continue
        if Vmax < 13*Vmin: continue   # S3 requires Vmax>=13 Vmin
        out.append(S)
    return out

# Also a targeted family: parked multiple-of-14 + small unit-residue core.
def gen_parked14_family(num=2500, seed=7, maxm=30):
    """S = P ∪ {large cluster} where exactly one element is a multiple of 14.
       Bias toward sets where 14m privately covers q=14 (and possibly q=7)."""
    random.seed(seed)
    out=[]
    tries=0
    while len(out)<num and tries<num*300:
        tries+=1
        m=random.randint(1,maxm)
        w=14*m
        if w<=13: continue
        # small part avoiding extra multiples of 14 (trivial since <14)
        ksmall=random.randint(8,11)
        P=sorted(random.sample(range(1,14), ksmall))
        # other large members, none a multiple of 14, to make k>=2
        klarge=13-ksmall-1
        if klarge<1: continue
        pool=[x for x in range(14,400) if x%14!=0 and x!=w]
        L=sorted(random.sample(pool, klarge))
        S=sorted(set(P+L+[w]))
        if len(S)!=13: continue
        if not is_prim(S): continue
        if not is_cov(S): continue
        # ensure w is the UNIQUE multiple of 14
        if sum(1 for v in S if v%14==0)!=1: continue
        Vmin=min(S); Vmax=max(S); kbig=sum(1 for v in S if v>13)
        if kbig<2: continue
        if Vmax<13*Vmin: continue
        out.append((S,w))
    return out

def analyze(S, parked=None):
    M,t=Mval(S)
    B,pairs=classify_tau(S,t,M)
    # j>=D/14 for crossing reps
    info=[]
    for va,vb,D,j in pairs:
        info.append((va,vb,D,j, 14*j>=D))   # 14j>=D <=> j>=D/14
    # is M>=1/14?
    ok = M>=F(1,14)
    # private q owners
    owners=private_q_owner(S)
    privq=[q for q in range(2,15) if len(owners[q])==1]
    # which runner privately owns q=14?
    own14 = owners[14][0] if len(owners[14])==1 else None
    own7  = owners[7][0]  if len(owners[7])==1 else None
    parked_in_binding = (parked in B) if parked is not None else None
    own14_in_binding = (own14 in B) if own14 is not None else None
    return dict(M=M,t=t,B=B,pairs=info,ok=ok,privq=privq,
                own14=own14,own7=own7,
                parked_in_binding=parked_in_binding,
                own14_in_binding=own14_in_binding)

if __name__=="__main__":
    print("="*70)
    print("PART 0: sanity on named cases")
    for S in [[1,2,3,5,7,8,9,10,11,12,13,38,42],
              [1,2,3,4,5,6,7,8,9,10,11,12,182],
              [1,2,3,4,5,6,7,8,9,10,11,13,84]]:
        a=analyze(S)
        print(f"S={S}\n  M={a['M']} t={a['t']} binding={a['B']} privq={a['privq']} own14={a['own14']}")
        print(f"  crossing pairs (va,vb,D,j, j>=D/14): {a['pairs']}")

    print("="*70)
    print("PART 1: random covering S3 sets")
    S3=gen_covering_S3(num=4000, seed=11, maxbig=400)
    print(f"generated {len(S3)} covering primitive S3 sets")
    n_break=0
    n_tau_is_crossing=0
    n_all_j_ok=0
    n_no_crossing=0
    minM=F(1)
    minMset=None
    # over all sets, does EVERY binding-pair crossing rep have j>=D/14?
    # (this is essentially M>=1/14 restated; we record honestly)
    n_parkedbinding=0; n_parkedtotal=0
    for S in S3:
        a=analyze(S)
        if a['M']<F(1,14): n_break+=1
        if a['M']<minM: minM=a['M']; minMset=S
        if len(a['pairs'])>0:
            n_tau_is_crossing+=1
            if all(flag for *_,flag in a['pairs']):
                n_all_j_ok+=1
        else:
            n_no_crossing+=1
        # parked = the (a) multiple of 14 if unique
        own14=a['own14']
        if own14 is not None:
            n_parkedtotal+=1
            if own14 in a['B']: n_parkedbinding+=1
    print(f"LRC breaks (M<1/14): {n_break}")
    print(f"min M = {minM} on {minMset}  (1/14={F(1,14)})  M*14={minM*14}")
    print(f"tau* expressible as binding-pair crossing j/D: {n_tau_is_crossing}/{len(S3)} (no-crossing {n_no_crossing})")
    print(f"  of those, all crossing reps have j>=D/14: {n_all_j_ok}/{n_tau_is_crossing}")
    print(f"unique-mult-of-14 sets: {n_parkedtotal}; parked runner is BINDING at tau*: {n_parkedbinding} ({n_parkedbinding}/{n_parkedtotal})")

    print("="*70)
    print("PART 2: parked-14 family (unique multiple of 14 = parked runner)")
    fam=gen_parked14_family(num=2500, seed=23, maxm=30)
    print(f"generated {len(fam)} parked-14 covering S3 sets")
    n_break=0; minM=F(1); minMset=None
    n_parkedbinding=0
    n_own14priv=0          # 14m privately owns q=14
    n_drop_below=0         # removing parked drops M below 1/14 (so parked is essential)
    n_drop_below_butpair=0 # ...yet a non-parked pair still binds in full set
    j_over_D14_fail=0
    examples_parked_not_binding=[]
    for S,w in fam:
        a=analyze(S, parked=w)
        if a['M']<F(1,14): n_break+=1
        if a['M']<minM: minM=a['M']; minMset=S
        if w in a['B']: n_parkedbinding+=1
        else:
            if len(examples_parked_not_binding)<6:
                examples_parked_not_binding.append((S,w,a['M'],a['t'],a['B']))
        if 14 in a['privq']: n_own14priv+=1
        for *_,flag in a['pairs']:
            if not flag: j_over_D14_fail+=1
        # essential test: M of S\{w}
        Sm=[v for v in S if v!=w]
        Mm,_=Mval(Sm)
        if Mm<F(1,14):
            n_drop_below+=1
    print(f"LRC breaks: {n_break}")
    print(f"min M = {minM} on {minMset}  M*14={minM*14}")
    print(f"parked runner 14m is BINDING at tau*: {n_parkedbinding}/{len(fam)}")
    print(f"14m privately owns q=14: {n_own14priv}/{len(fam)}")
    print(f"removing parked 14m drops M<1/14 (parked essential to LRC): {n_drop_below}/{len(fam)}")
    print(f"crossing reps with j<D/14 (would be LRC fail): {j_over_D14_fail}")
    print("examples where parked NOT binding:")
    for ex in examples_parked_not_binding:
        print("   ", ex)
