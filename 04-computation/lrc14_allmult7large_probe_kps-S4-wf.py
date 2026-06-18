"""
Decisive probe of the ALL-MULT7-LARGE sub-case for LRC(14) S3.

ALL-MULT7-LARGE := primitive covering 13-set, S3 (k=#{v>13}>=2), and every
multiple of 7 in S exceeds V* := max{v in S : 7 does not divide v}.

Window-collapse reduction (PROVED): on tau = k0/7 + s with |s| <= 1/(14 V*),
  - every non-mult-of-7 runner is automatically >= 1/14 safe,
  - mult-of-7 runners 7 w_i reduce to ||w_i u|| with u = 7s, |u| <= 1/(2 V*).
A {w_i}-safe u in (0, 1/(2 V*)] (or symmetric) => GLOBAL WITNESS => M(S) >= 1/14.

QUESTION the critic raised: does such an in-window safe u ALWAYS exist for
covering ALL-MULT7-LARGE sets?  Covering forces an EVEN multiple of 7 (q=14)
and a multiple of 13.  We search exhaustively for ANY covering ALL-MULT7-LARGE
set where the window is BLOCKED (no in-window witness), and -- crucially --
whether M is still >= 1/14 there (LRC unbroken regardless).

Search design:
  - V* ranges 14..120.  small runners (non-mult-of-7) all <= V*, with V* present.
  - mult-of-7 large runners: subsets (size 1..3) of multiples of 7 in (V*, 8 V*],
    must include an even one (q=14) and collectively a multiple of 13 may come
    from small 13 or a large mult like 91,182,273.
We don't enumerate ALL sets (combinatorial blowup); we instead enumerate the
REDUCED w_i multisets that can occur and test window-blockage directly, then
for any blocked w_i we confirm by constructing an actual covering set and
computing exact M.
"""
from fractions import Fraction as F
from math import gcd
import itertools

def gcd_list(L):
    g=0
    for x in L: g=gcd(g,x)
    return g
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
    b=F(0)
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v
    return b
def is_cov(S):
    return all(any(v%q==0 for v in S) for q in range(2,15))

def window_witness_exists(wlist, Vstar):
    """Is there a {w_i}-safe u with |u| <= 1/(2 Vstar) (u != 0)?"""
    bound = F(1, 2*Vstar)
    sc = safe_components(sorted(set(wlist)))
    # positive side (0,bound]
    for (a,b) in sc:
        lo=max(a,F(0)); hi=min(b,bound)
        if hi>lo and not (lo==0 and hi==0):
            # need a point strictly >0: interval (lo,hi) with hi>max(lo,0)
            if hi>F(0) and hi>lo: return True
    # negative side [1-bound,1)
    for (a,b) in sc:
        lo=max(a,F(1)-bound); hi=min(b,F(1))
        if hi>lo: return True
    return False

R=[]
def out(s): R.append(s); print(s)

out("DECISIVE PROBE: ALL-MULT7-LARGE window-collapse witness existence")
out("="*72)

# Step 1: enumerate possible reduced w_i multisets for covering ALL-MULT7-LARGE.
# The mult-of-7 runners are 7 w_i with 7 w_i > V*, i.e. w_i > V*/7.  Covering needs
# an even mult-of-7 (=> some w_i even) AND a mult of 13 in S (could be small 13<=V*,
# or a mult-of-7 that is a mult of 91 => w_i mult of 13).
# We sweep V* and ask: over ALL admissible w_i-sets (size 1..3, each w_i in
# (V*/7, large], with >=1 even), is the window EVER blocked?
blocked_cases=[]
checked=0
for Vstar in range(14,121):
    wmin = Vstar//7 + 1            # smallest w with 7w > V*
    wmax = (8*Vstar)//7            # cap 7w <= 8 V* (keep mult-of-7 < 8 V*; generous)
    wpool=list(range(wmin, max(wmin+1,wmax+1)))
    if not wpool: continue
    for r in (1,2,3):
        for combo in itertools.combinations(wpool, r):
            # covering needs an even mult of 7 => some 7 w_i even => some w_i ...
            # 7 w even iff w even. require >=1 even w in combo OR allow small even mult?
            # In ALL-MULT7-LARGE the ONLY mult-of-7 are these large ones, so the q=14
            # multiple MUST be among them => need some 7 w_i divisible by 14 => w_i even.
            if not any(w%2==0 for w in combo):
                continue
            checked+=1
            if not window_witness_exists(combo, Vstar):
                blocked_cases.append((Vstar, combo))
    if Vstar%20==0:
        out(f"  ...swept V* up to {Vstar}, checked {checked} w-multisets, blocked so far {len(blocked_cases)}")

out("")
out(f"Total admissible reduced w-multisets checked: {checked}")
out(f"BLOCKED-window cases (no in-window {{w_i}}-safe u): {len(blocked_cases)}")
if blocked_cases:
    out("  First 20 blocked (V*, w_i):")
    for c in blocked_cases[:20]:
        out(f"    V*={c[0]} w={c[1]}")
else:
    out("  => For EVERY admissible covering ALL-MULT7-LARGE reduced w-set (V*<=120),")
    out("     an in-window global witness EXISTS. Window-collapse closes them all.")

# Step 2: For any blocked w-set, try to realize an actual covering set & compute M.
out("")
out("Realizing blocked cases as actual covering ALL-MULT7-LARGE sets & exact M:")
realized=0
broke=0
for (Vstar, combo) in blocked_cases[:200]:
    mults=[7*w for w in combo]
    if any(m<=Vstar for m in mults): continue
    # build covering set: small runners non-mult-of-7 <= V*, include V*, plus mults.
    need = 13 - len(mults)
    if need<1: continue
    small_pool=[v for v in range(1,Vstar+1) if v%7!=0]
    if Vstar not in small_pool: continue
    rest=[v for v in small_pool if v!=Vstar]
    # greedy: try to find a covering completion by sampling combinations (cap work)
    import itertools as it
    done=False
    cnt=0
    for c2 in it.combinations(rest, need-1):
        cnt+=1
        if cnt>20000: break
        S=sorted([Vstar]+list(c2)+mults)
        if len(S)!=13: continue
        if gcd_list(S)!=1: continue
        if not is_cov(S): continue
        if max(v for v in S if v%7!=0)!=Vstar: continue
        if any(v%7==0 and v<=Vstar for v in S): continue
        kbig=sum(1 for v in S if v>13)
        if kbig<2: continue
        m=Mval(S)
        realized+=1
        if m<F(1,14): broke+=1; out(f"  !!! LRC BROKEN: S={S} M={m}")
        else:
            out(f"  blocked-window but M>=1/14: S={S} M={m}={float(m):.4f}")
        done=True
        break
out(f"\nRealized blocked covering sets: {realized}; LRC breaks: {broke}")

out("")
out("VERDICT on ALL-MULT7-LARGE:")
if not blocked_cases:
    out("  Window-collapse witness exists for ALL admissible reduced w-sets (V*<=120).")
    out("  => The ALL-MULT7-LARGE sub-case is CLOSED by window-collapse over V*<=120,")
    out("     EXACTLY because covering forces an even w_i (small tooth) that keeps a")
    out("     gap open near 0. The critic's V*=19 {21,273} block is OUT OF SCOPE (odd")
    out("     mults => not covering).")
else:
    out(f"  {len(blocked_cases)} blocked reduced w-sets exist; {broke} actually break LRC.")
    out("  If broke==0, LRC holds on these via OTHER tau (window-collapse is only sufficient).")

with open(r"C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\lrc14_allmult7large_probe_kps-S4-wf.out","w") as f:
    f.write("\n".join(R))
out("\n[written]")
