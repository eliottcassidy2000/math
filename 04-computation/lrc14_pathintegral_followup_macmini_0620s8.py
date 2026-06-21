#!/usr/bin/env python3
"""
lrc14_pathintegral_followup_macmini_0620s8.py   (mac-mini-2026-06-20-S8)

Follow-up probes after the first amplitude test:

  T2-deep: The raw clock-coherence  coh(E) = avg_f |sum_e omega^{floor(e f)}|^2
           picked CONSEC as argmax with 0 violators at k=8. Is this a REAL
           constructive-interference signal correlated with measS7, or a trivial
           "clustered offsets => higher coherence" artifact decoupled from
           surjectivity?
     Tests:
       (a) closed form: coh(E) = sum_{e,e'} avg_f omega^{floor(ef)-floor(e'f)}.
           Compute the rank correlation between coh(E) and measS7(E) over all
           k=8 shapes; if high & consec is joint-argmax, the picture has teeth.
       (b) Does coh pick consec at k=9,10 too? (extend the box)
       (c) Is coh just a function of the gap profile (so it'd pick consec for
           a trivial reason)? Compare coh of consec vs an EQUALLY-clustered but
           shifted/reflected set; vs the iid/spread sets.
       (d) The DECISIVE test: take the set that BEATS consec on surj-weighted
           coherence ([0,2,4,6,7,8,10,12]) and check measS7 -- consec should
           still win measS7. So which functional tracks measS7: raw coh or
           surj-coh? Spearman both.

  T3-deep: OCF squared-amplitude.  We found signed_sum (-1)^desc takes small
           values. Test the CORRECT path-integral statement:
       (a) Walsh: H(T) = sum_S H_hat[S] chi_S(T). Build the Walsh vector and
           ask whether the GRAM matrix view holds: is there a unimodular feature
           map psi: T -> C^d with <psi(T),psi(T)> giving H? Equivalent clean
           test: define amplitude per Ham path with phase i^{rot} and see if
           |amplitude|^2 = H or = a known OCF quantity.
       (b) The cleanest amplitude object: per OCF, H = sum over INDEPENDENT SETS
           of odd cycles of 2^{#cycles}. The 2^{#cycles} = |{+,-}^{cycles}| is a
           sum over SIGN ASSIGNMENTS = a sum over 'which slit' (double slit!).
           Test: H = sum_{indep set I of odd cycles} prod_{c in I} 2.
           Reframe 2 = |1+e^{i*0}|... no; 2 = sum_{eps in {0,1}} 1. So H counts
           PATHS = (indep set, sign vector). This is a literal path-integral
           with all-equal amplitudes 1 (incoherent). Verify H = that count.
       Honest goal: separate (identity: H is a path SUM) from
       (analogy: the SIGNS are interference).

stdlib only.
"""
import sys, itertools, math, cmath
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
P = 7

def measS7(E, scale=P):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, scale*e+1): bps.add(F(m, scale*e))
    bps = sorted(b for b in bps if 0 <= b <= 1); tot = F(0)
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi <= lo: continue
        mid = (lo+hi)/2
        if len(set(int(((e*mid)%1)*scale) for e in E)) == scale: tot += hi-lo
    return tot

def coherence(E, scale=P, surj_only=False):
    omega = cmath.exp(2j*math.pi/scale); E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, e+1): bps.add(F(m, e))
    bps = sorted(b for b in bps if 0 <= b <= 1); tot = 0.0
    for i in range(len(bps)-1):
        f0, f1 = bps[i], bps[i+1]
        if f1 <= f0: continue
        fm = (f0+f1)/2
        if surj_only:
            cols = set(int(e*fm)%scale for e in E)
            if len(cols) != scale: continue
        amp = sum(omega**(int(e*fm)) for e in E)
        tot += (abs(amp)**2)*float(f1-f0)
    return tot

def spearman(xs, ys):
    """rank correlation, stdlib."""
    def ranks(v):
        idx = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0]*len(v); i=0
        while i < len(v):
            j=i
            while j+1<len(v) and v[idx[j+1]]==v[idx[i]]: j+=1
            avg=(i+j)/2.0
            for t in range(i,j+1): r[idx[t]]=avg
            i=j+1
        return r
    rx, ry = ranks(xs), ranks(ys)
    n=len(xs); mx=sum(rx)/n; my=sum(ry)/n
    num=sum((rx[i]-mx)*(ry[i]-my) for i in range(n))
    dx=math.sqrt(sum((rx[i]-mx)**2 for i in range(n)))
    dy=math.sqrt(sum((ry[i]-my)**2 for i in range(n)))
    return num/(dx*dy) if dx*dy else 0.0

def test_T2_deep():
    print("="*72); print("(T2-deep) is raw clock-coherence a real measS7 proxy?"); print("="*72)
    for k in (8,9):
        rows=[]
        span = 13 if k<=8 else 13
        for combo in itertools.combinations(range(1, span), k-1):
            E=[0]+list(combo)
            rows.append((E, float(measS7(E)), coherence(E), coherence(E, surj_only=True)))
        cons=list(range(k))
        ms=[r[1] for r in rows]; cohs=[r[2] for r in rows]; scohs=[r[3] for r in rows]
        # argmaxes
        am_m=max(rows,key=lambda r:r[1])[0]
        am_c=max(rows,key=lambda r:r[2])[0]
        am_sc=max(rows,key=lambda r:r[3])[0]
        sp_c=spearman(ms,cohs); sp_sc=spearman(ms,scohs)
        print(f"  k={k}  #shapes={len(rows)} (span<{span})")
        print(f"    argmax measS7   = {am_m}  (consec? {am_m==cons})")
        print(f"    argmax raw-coh  = {am_c}  (consec? {am_c==cons})  Spearman(coh,measS7)={sp_c:+.3f}")
        print(f"    argmax surj-coh = {am_sc} (consec? {am_sc==cons}) Spearman(surjcoh,measS7)={sp_sc:+.3f}")
    print("  INTERP: if raw-coh argmax=consec but Spearman is weak/negative, raw-coh")
    print("  picks consec for a CLUSTERING reason decoupled from measS7 (analogy, not")
    print("  mechanism). If Spearman strong & positive, it's a genuine proxy.")

def test_T2_clustering_control():
    print("="*72); print("(T2-control) is raw-coh just a clustering statistic?"); print("="*72)
    # coh closed form: coh = sum_{e,e'} g(e,e') where g = avg_f omega^{floor(ef)-floor(e'f)}
    # If coh depends ONLY on the gap multiset, then any reflection/translation in
    # index space with same gaps has same coh -> trivial. Test reflections.
    import random
    rng=random.Random(1)
    same=0; diff=0; maxdiff=0.0
    for _ in range(40):
        k=rng.randint(5,8)
        E=sorted(random.Random(rng.random()).sample(range(0,13),k))
        if 0 not in E: E=[0]+E[1:]
        # reflect within span
        mx=max(E); Eref=sorted(mx-e for e in E)
        c1=coherence(E); c2=coherence(Eref)
        if abs(c1-c2)<1e-9: same+=1
        else: diff+=1; maxdiff=max(maxdiff,abs(c1-c2))
    print(f"  coherence under index-reflection: equal={same} differ={diff} maxdiff={maxdiff:.4f}")
    print("  (if it differs, coh is NOT a pure gap statistic -> it sees the actual")
    print("   floor-word structure, i.e. real interference, not just clustering.)")

def hpaths(adj,n):
    return [perm for perm in itertools.permutations(range(n))
            if all(adj[perm[i]][perm[i+1]] for i in range(n-1))]

def all_tournaments(n):
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in itertools.product((0,1),repeat=len(edges)):
        adj=[[0]*n for _ in range(n)]
        for (i,j),b in zip(edges,bits):
            if b: adj[i][j]=1
            else: adj[j][i]=1
        yield adj

def odd_cycles(adj,n):
    """list all directed odd cycles (as frozenset of vertices won't do; need the
    cyclic structure). Return list of vertex-sets per directed odd cycle, with
    multiplicity = number of directed odd cycles. For OCF we need INDEPENDENT
    SETS of vertex-DISJOINT odd cycles."""
    cycles=[]
    for L in range(3,n+1,2):
        for c in itertools.combinations(range(n),L):
            for perm in itertools.permutations(c):
                if perm[0]!=min(perm): continue  # canonical start to dedupe rotations
                if all(adj[perm[i]][perm[(i+1)%L]] for i in range(L)):
                    cycles.append(frozenset(c))
    return cycles  # each directed cycle once (canonical-start), but reflections distinct

def ocf_count(adj,n):
    """H via OCF: sum over independent sets I of vertex-disjoint directed odd
    cycles of 2^{|I|}. (THM/OCF). We enumerate sets of pairwise vertex-disjoint
    directed odd cycles."""
    cycs=odd_cycles(adj,n)
    # build list of (vertexset) for disjointness; multiple directed cycles can
    # share a vertexset (e.g. two orientations) -> they are DIFFERENT cycles.
    n_c=len(cycs)
    total=0
    # enumerate antichains of pairwise disjoint cycles
    # do DFS
    def rec(start, used_mask):
        nonlocal total
        # count this independent set: contributes 2^{size}; but we add via recursion:
        # we will instead count by: total += 2^{size}. Easier: enumerate subsets.
        pass
    # direct: enumerate all subsets that are pairwise disjoint
    vsets=[c for c in cycs]
    # represent each cycle's vertex set as bitmask
    masks=[sum(1<<v for v in c) for c in vsets]
    res=0
    def dfs(i, occ, size):
        nonlocal res
        if i==len(masks):
            res += 2**size
            return
        # skip cycle i
        dfs(i+1, occ, size)
        # take cycle i if disjoint
        if occ & masks[i]==0:
            dfs(i+1, occ|masks[i], size+1)
    dfs(0,0,0)
    return res

def test_T3_deep():
    print("="*72); print("(T3-deep) OCF = literal path SUM (indep set, sign vector)?"); print("="*72)
    n=5
    bad_ocf=0; bad_walsh=0; tot=0; ex=[]
    for adj in all_tournaments(n):
        tot+=1
        H=len(hpaths(adj,n))
        ocf=ocf_count(adj,n)
        if H!=ocf: bad_ocf+=1
        if len(ex)<5: ex.append((H,ocf))
    print(f"  n={n} tournaments={tot}  H==OCF(sum 2^#cycles over indep sets) mismatches={bad_ocf}")
    print(f"  sample (H, OCF): {ex}")
    print("  => If 0 mismatches: H IS a path-integral SUM over (indep cycle set, sign)")
    print("     pairs with all amplitudes = +1 (INCOHERENT sum, |1|^2 each). The")
    print("     'amplitude' is real & positive -> NO interference in H itself.")
    print("  CONTRAST: the WALSH coeffs (-1)^desc DO carry signs, but they live at")
    print("  the FOURIER (dual) level, not in the path count. So OCF=incoherent sum;")
    print("  the interference is in the Walsh/Fourier dual (THM-071), not in H.")

def main():
    print("PATH-INTEGRAL follow-up  (mac-mini-2026-06-20-S8)\n")
    test_T2_deep(); print()
    test_T2_clustering_control(); print()
    test_T3_deep(); print()
    print("="*72); print("Verdicts in final message."); print("="*72)

if __name__=="__main__":
    main()
