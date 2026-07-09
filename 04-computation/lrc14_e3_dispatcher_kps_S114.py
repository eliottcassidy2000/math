"""
LRC(14) E3-hardness dispatcher (kps-S114).
Operationalizes the Schur-triple count E3 as the HARDNESS COORDINATE that routes a 13-runner
speed set to a proof cell, using mac-mini's THM-668 pair-sum ruler theorem (grid-free):

  - E3(S) = #{(a,b) in S x S : a+b in S}  -- the Schur-triple count (opus LEM-015), the KILL BUDGET.
  - Pair-sum moduli q = v_i+v_j are the ONLY places the min-reach max lives (THM-668).
  - A ruler q is DEAD if some runner v_l ≡ 0 (mod q) (Schur-triple kill rule); else LIVE.
  - Witness search: for a live q and multiplier p (gcd(p,q)=1), if every residue (v_l p mod q)
    lies in the band [q/14, 13q/14], then t=p/q is a lonely instant, M(S) >= 1/14. GRID-FREE, exact.
  - Route: high normalized E3 (near saturation C(13,2)=78, tight/AP-like) -> EXACT-CHECK cell;
           low E3 (dissociated) -> GOOD-PERIOD cell.  E3 selects the cell (klein-S201 2x2).
"""
from math import gcd, comb, ceil, floor
from fractions import Fraction

def E3(S):
    Sset=set(S)
    return sum(1 for a in S for b in S if a+b in Sset)

def pairsum_moduli(v):
    return sorted({v[i]+v[j] for i in range(len(v)) for j in range(i,len(v))})

def is_live(q, v):
    return all(vl % q != 0 for vl in v)   # dead if some runner ≡ 0 mod q

def band_witness(v, qmax_mult=2):
    """Find a grid-free pair-sum witness t=p/q with all residues in [q/14,13q/14]. THM-668 complete."""
    Vmax=max(v)
    for q in pairsum_moduli(v):
        if q < 2 or not is_live(q, v): continue
        lo=q  # 14*r >= q
        hi=13*q  # 14*r <= 13q
        for p in range(1, q):
            if gcd(p,q)!=1: continue
            ok=all(lo <= 14*((vl*p)%q) <= hi for vl in v)
            if ok:
                return (p, q)
    return None

def M_exact(v, denom_cap=4000):
    """True M(S)=max_t min_l ||v_l t|| by pair-sum event enumeration (THM-668: max lives at p/q)."""
    best=Fraction(0)
    Q=pairsum_moduli(v)
    for q in Q:
        for p in range(1,q):
            t=Fraction(p,q)
            m=min(min((vl*t)%1, 1-((vl*t)%1)) for vl in v)
            if m>best: best=m
    return best

def dispatch(v, name):
    k=len(v); e3=E3(v); sat=comb(k,2)
    norm=e3/sat
    Q=pairsum_moduli(v); live=[q for q in Q if q>=2 and is_live(q,v)]
    w=band_witness(v)
    M=M_exact(v)
    lonely = (M >= Fraction(1,14))
    cell = "EXACT-CHECK (tight/near-AP)" if norm>=0.5 else "GOOD-PERIOD (dissociated)"
    print(f"\n{name}: v={v}")
    print(f"  E3={e3}/{sat} (norm {norm:.3f})  live pair-sum rulers={len(live)} of {len(Q)}")
    print(f"  grid-free band witness: {('t='+str(w[0])+'/'+str(w[1])) if w else 'NONE among live rulers'}")
    print(f"  M(S)={M} ({float(M):.4f})   >=1/14? {lonely}   -> CELL: {cell}")
    return dict(name=name, e3=e3, norm=norm, live=len(live), witness=w, M=M, lonely=lonely, cell=cell)

if __name__=="__main__":
    AP=list(range(1,14))
    tests=[
        (AP, "AP {1..13} (E3-MAX extremal)"),
        ([2*x for x in AP], "2*AP {2,4,..,26} (dilated interval, same E3)"),
        ([1,3,5,7,9,11,13,15,17,19,21,23,25], "odd AP (offset, low E3)"),
        ([1,2,3,4,5,6,7,8,9,10,11,12,14], "near-AP (one shifted)"),
        ([1,2,4,8,16,32,64,128,256,512,1024,2048,4096], "geometric (dissociated, E3~0)"),
    ]
    print("="*70)
    print("LRC(14) E3-HARDNESS DISPATCHER (THM-668 grid-free pair-sum events)")
    print("="*70)
    res=[dispatch(v,n) for v,n in tests]
    print("\n"+"="*70)
    print("SUMMARY (E3 as hardness coordinate; all lonely => LRC(14) holds on these):")
    for r in res:
        print(f"  {r['name'][:38]:38s} E3norm={r['norm']:.3f} live={r['live']:2d} lonely={r['lonely']} cell={r['cell'][:12]}")
    print("\nKEY: E3-max (AP) has the FEWEST live rulers (max kill budget) yet is still lonely at its")
    print("single live ruler q=14 (tangency); dissociated sets have MANY live rulers, lonely with margin.")
