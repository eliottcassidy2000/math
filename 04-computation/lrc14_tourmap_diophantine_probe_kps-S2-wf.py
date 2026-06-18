"""
lrc14_tourmap_diophantine_probe_kps-S2-wf.py

FOLLOW-UP PROBE.  The control proved M2/M5 forbidding is a RULE artifact.
Two finer questions remain, both honest:

 (Q1) ARGMAX-vs-ALLCAND sub-image:  Even if a map realizes the same set under
      "all binding crossings", does restricting to the SINGLE lonely-optimum tau
      (argmax) realize a STRICTLY SMALLER set?  If yes, the lonely optimum
      itself selects a sub-class-set -- a weak but genuine loneliness signal.
      We report argmax-image vs allcand-image for the live maps (M2,M3,M3b,M5).

 (Q2) A new candidate map M7 built directly from THM-524 BINDING structure on
      PAIRS, designed to (a) be non-transitive and (b) couple to loneliness via
      *which pair binds* at the optimum.  Vertices = runners; arc i->j iff at the
      lonely optimum tau the runner i is strictly closer to 0 than j AND they lie
      on opposite sides ("i passes j toward the integer"), else by a CF tiebreak.
      We measure its image + run the same loneliness-vs-rule control.

We also report, for the live maps at n=4,5, the per-iso-class HIT COUNT under the
argmax (true lonely) family, to expose whether the lonely optimum is biased
toward the transitive class (H=1) -- the project's central object.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def signed_frac(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else r - 1
def nrm_frac01(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r
def cont_frac(p, q):
    a = []
    while q:
        a.append(p // q); p, q = q, p % q
    return a
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def g(S, t): return min(nrm(v*t) for v in S)
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at
def gcd_list(xs):
    g0 = 0
    for x in xs: g0 = gcd(g0, x)
    return g0
def speed_sets(n, maxspeed):
    return [S for S in combinations(range(1, maxspeed+1), n) if gcd_list(S) == 1]

PAIRS = {n: list(combinations(range(n), 2)) for n in (3,4,5)}
def arcs_to_mask(arcs, n):
    idx = {pr:k for k,pr in enumerate(PAIRS[n])}; mask=0
    for (i,j) in arcs:
        a,b=(i,j) if i<j else (j,i); k=idx[(a,b)]
        if (i,j)==(a,b): mask|=(1<<k)
    return mask
def mask_to_adj(mask,n):
    A=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(PAIRS[n]):
        if (mask>>k)&1: A[i][j]=1
        else: A[j][i]=1
    return A
def adj_num_3cyc(A,n):
    c=0
    for a,b,cc in combinations(range(n),3):
        outs={a:0,b:0,cc:0}
        for (x,y) in combinations([a,b,cc],2):
            if A[x][y]==1: outs[x]+=1
            else: outs[y]+=1
        if sorted(outs.values())==[1,1,1]: c+=1
    return c
def adj_score(A,n): return tuple(sorted(sum(A[i]) for i in range(n)))
def relabel_mask(mask,perm,n):
    A=mask_to_adj(mask,n); B=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n): B[perm[i]][perm[j]]=A[i][j]
    m2=0
    for k,(i,j) in enumerate(PAIRS[n]):
        if B[i][j]==1: m2|=(1<<k)
    return m2
def build_canon_table(n):
    perms=list(permutations(range(n))); npairs=len(PAIRS[n])
    canon_of={}; class_info={}
    for mask in range(1<<npairs):
        if mask in canon_of: continue
        orbit=set(relabel_mask(mask,p,n) for p in perms)
        rep=min(orbit)
        if rep not in class_info:
            A=mask_to_adj(rep,n); class_info[rep]=(adj_num_3cyc(A,n),adj_score(A,n))
        for mm in orbit: canon_of[mm]=rep
    reps=sorted(set(canon_of.values())); repid={r:i for i,r in enumerate(reps)}
    return {m:repid[canon_of[m]] for m in canon_of}, {repid[r]:class_info[r] for r in reps}, len(reps)

# ---- maps ----
def method2(S, tau):
    n=len(S); dist=[nrm(v*tau) for v in S]
    side=[1 if signed_frac(v*tau)>0 else (-1 if signed_frac(v*tau)<0 else 0) for v in S]
    q=[v*nrm(v*tau) for v in S]; arcs=[]
    for i in range(n):
        for j in range(i+1,n):
            if side[i]==side[j]:
                if dist[i]==dist[j]: return None
                arcs.append((i,j) if dist[i]<dist[j] else (j,i))
            else:
                if q[i]==q[j]: return None
                arcs.append((i,j) if q[i]<q[j] else (j,i))
    return arcs
def method3(S, tau):
    n=len(S); arcs=[]
    for i in range(n):
        for j in range(i+1,n):
            a,b=max(S[i],S[j]),min(S[i],S[j]); P=sum(cont_frac(a,b))
            big=i if S[i]>S[j] else j; sml=j if big==i else i
            arcs.append((sml,big) if P%2==0 else (big,sml))
    return arcs
def method3b(S, tau):
    n=len(S); depth=[]
    for v in S:
        d=nrm(v*tau)
        depth.append(None if d==0 else int(F(1,1)//d))
    arcs=[]
    for i in range(n):
        for j in range(i+1,n):
            di,dj=depth[i],depth[j]
            if di is None and dj is None: return None
            if di is None: arcs.append((i,j)); continue
            if dj is None: arcs.append((j,i)); continue
            if di>dj: arcs.append((i,j))
            elif dj>di: arcs.append((j,i))
            else:
                a,b=max(S[i],S[j]),min(S[i],S[j]); P=sum(cont_frac(a,b))
                big=i if S[i]>S[j] else j; sml=j if big==i else i
                arcs.append((sml,big) if P%2==0 else (big,sml))
    return arcs
def method5(S, tau):
    n=len(S); third=[]
    for v in S:
        f=nrm_frac01(v*tau); third.append(0 if f<F(1,3) else (1 if f<F(2,3) else 2))
    dist=[nrm(v*tau) for v in S]; arcs=[]
    for i in range(n):
        for j in range(i+1,n):
            ti,tj=third[i],third[j]
            if ti==tj:
                if dist[i]==dist[j]: return None
                arcs.append((i,j) if dist[i]<dist[j] else (j,i))
            else:
                arcs.append((i,j) if (ti-tj)%3==2 else (j,i))
    return arcs

# ---- M7: binding side-crossing on runners (THM-524 flavored) ----
# i->j iff: signed_frac(v_i tau) and signed_frac(v_j tau) have opposite sign and
#   |v_i tau| (dist) < |v_j tau| (i is the one "binding tighter on the integer
#   from its side"); if same sign, the FARTHER one beats the nearer (the one
#   still chasing); ties -> CF parity of speeds. Mixing "opposite-side: near wins"
#   with "same-side: far wins" is the deliberate cycle generator.
def method7(S, tau):
    n=len(S); dist=[nrm(v*tau) for v in S]
    side=[1 if signed_frac(v*tau)>0 else (-1 if signed_frac(v*tau)<0 else 0) for v in S]
    arcs=[]
    for i in range(n):
        for j in range(i+1,n):
            if side[i]!=side[j]:
                if dist[i]==dist[j]:  # tie -> CF
                    a,b=max(S[i],S[j]),min(S[i],S[j]); P=sum(cont_frac(a,b))
                    big=i if S[i]>S[j] else j; sml=j if big==i else i
                    arcs.append((sml,big) if P%2==0 else (big,sml)); continue
                arcs.append((i,j) if dist[i]<dist[j] else (j,i))
            else:
                if dist[i]==dist[j]:
                    a,b=max(S[i],S[j]),min(S[i],S[j]); P=sum(cont_frac(a,b))
                    big=i if S[i]>S[j] else j; sml=j if big==i else i
                    arcs.append((sml,big) if P%2==0 else (big,sml)); continue
                arcs.append((i,j) if dist[i]>dist[j] else (j,i))  # FAR wins
    return arcs

LIVE = {"M2":method2,"M3":method3,"M3b":method3b,"M5":method5,"M7":method7}

def image_counts(method, n, maxspeed, mode, Qmax=30):
    """mode: 'argmax','allcand','unc'. Returns Counter id->hits."""
    from collections import Counter
    mt, info, nc = TAB[n]
    cnt=Counter()
    for S in speed_sets(n, maxspeed):
        if mode=='argmax':
            taus=[M(S)[1]]
        elif mode=='allcand':
            taus=sorted(cand(S))
        else:
            taus=[F(a,Q) for Q in range(2,Qmax+1) for a in range(1,Q)]
        for tau in taus:
            arcs=method(S,tau)
            if arcs is None: continue
            cnt[mt[arcs_to_mask(arcs,n)]]+=1
    return cnt, info, nc

TAB={}
def main():
    print("PROBE: argmax-vs-allcand sub-image + new M7 + hit-count bias\n")
    for n in (4,5):
        mt,info,nc=build_canon_table(n); TAB[n]=(mt,info,nc)

    for name,method in LIVE.items():
        print("="*64); print(name); print("="*64)
        for n in (4,5):
            ms=12 if n==5 else 14
            mt,info,nc=TAB[n]
            ca,_,_=image_counts(method,n,ms,'argmax')
            cc,_,_=image_counts(method,n,ms,'allcand')
            cu,_,_=image_counts(method,n,ms,'unc')
            ima=set(ca); imc=set(cc); imu=set(cu)
            allids=set(range(nc))
            # loneliness-forbidden = realized unconstrained but NOT at argmax
            lonely_forb = (imu - ima)
            rule_forb = allids - imu
            print(f" n={n}: free={nc} | argmax-image={len(ima)} "
                  f"allcand-image={len(imc)} unc-image={len(imu)}")
            print(f"      forbidden@argmax not by rule (unc-realized but argmax-misses): "
                  f"{len(lonely_forb)}")
            for cid in sorted(lonely_forb):
                c3,ss=info[cid]
                # is it realized under allcand (binding crossings)?
                in_allcand = cid in imc
                print(f"        argmax-misses class #3cyc={c3} scores={ss} "
                      f"(allcand-realized={in_allcand})")
            # transitive-class hit fraction at argmax (the H=1 trivial class)
            # transitive class = the one with c3=0
            trans_id=[cid for cid in info if info[cid][0]==0][0]
            tot=sum(ca.values()); th=ca.get(trans_id,0)
            print(f"      argmax hit-count: total={tot}, transitive(H=1)={th} "
                  f"({F(th,tot) if tot else 0})")
        print()

if __name__=="__main__":
    main()
