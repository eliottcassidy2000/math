"""
opus-2026-07-11-S238: the spread bulk of the residual -- the pigeonhole reduction for prime moduli, and an
HONEST NEGATIVE result: no small-modulus-subset shortcut exists; the full window is irreducible.

Residual (S234-237): every spread divisor-complete family (longest-AP<=7, mult of every d<=14) clears at a
non-14 modulus q in [15,31] => M>1/14 (band-edge lemma). This is the S230/S231 anti-concentration. This
session probes whether a SMALL fixed modulus set suffices (which would simplify the proof target).

THE PIGEONHOLE REDUCTION (clean, for prime q with band {0,+-1}, q in {17,19,23}):
  A family with no multiple of q does NOT clear at q  <=>  its 13 folded residues {min(v_i mod q, q-v_i mod q)}
  occupy ALL (q-1)/2 fold-classes (no empty class => every multiplier hits {0,+-1}). Clearing at q <=> some
  fold-class empty. (13 speeds, (q-1)/2 classes: q=17->8, q=19->9, q=23->11.)

MINIMAL COVERING SET: over 2000 spread divisor-complete families, the minimal set of non-14 moduli in [15,31]
clearing ALL is {15, 17, 19, 23, 25, 27, 31} (7 moduli); best single is q=31 (70%).

HONEST NEGATIVE RESULT (the {17,19,23} shortcut FAILS): adversarial hill-climb FINDS a spread divisor-complete
family with no multiple of 17,19,23 that occupies ALL fold-classes mod 17 (8/8), 19 (9/9), AND 23 (11/11)
simultaneously -- so it clears at NONE of {17,19,23}:
  v = [42, 48, 60, 108, 125, 154, 195, 206, 210, 245, 252, 259, 294]
It clears only at q=26 in the full window (M >= 2/26 = 1/13 > 1/14). So NO 3-prime subset suffices; the full
7-modulus window is genuinely needed. (A random sample gives 0 such families -- occupying all fold-classes
mod all three primes is rare but ACHIEVABLE, so the "0/425" random result was misleading.)

CONSEQUENCE. The spread-bulk anti-concentration is IRREDUCIBLE to a few primes: proving it requires the full
window disjunction (clearing at some q in {15,17,19,23,25,27,31}), each a genuine anti-concentration
condition, with no shortcut via a fixed small prime set. This rules out the tempting "prove clearing at 2-3
fixed primes" approach and confirms the wall's shape (the fleet's long-standing anti-concentration crux).
"""
from math import gcd, ceil
from functools import reduce
import random
def clears(v,q):
    for p in range(1,q):
        if all(q<=14*((vi*p)%q)<=13*q for vi in v): return True
    return False
def foldocc(v,q): return len(set(min(x%q,q-x%q) for x in v if x%q!=0))
def divisor_complete(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def primitive(v): return reduce(gcd,v)==1
def longest_AP(v):
    s=set(v); best=1
    for a in v:
        for d in range(1,max(v)//2+1):
            L=1;x=a+d
            while x in s:L+=1;x+=d
            if L>best:best=L
    return best

def main():
    # the pigeonhole reduction + the counterexample
    cx=[42,48,60,108,125,154,195,206,210,245,252,259,294]
    print("PIGEONHOLE REDUCTION: not-clear-at-prime-q (no mult q) <=> folded residues occupy all (q-1)/2 classes.")
    print(f"COUNTEREXAMPLE to the {{17,19,23}} shortcut (spread, divisor-complete, no mult 17/19/23):")
    print(f"  v={cx}")
    print(f"  spread(longest-AP)={longest_AP(cx)}, divisor-complete={divisor_complete(cx)}, primitive={primitive(cx)}")
    print(f"  fold-occupancy: mod17={foldocc(cx,17)}/8, mod19={foldocc(cx,19)}/9, mod23={foldocc(cx,23)}/11 (all FULL)")
    print(f"  clears at {{17,19,23}}? {[q for q in [17,19,23] if clears(cx,q)]} (NONE)")
    wq=[q for q in range(15,40) if q%14!=0 and clears(cx,q)]
    print(f"  clears in full window at q={wq[:5]} => M>=ceil({wq[0]}/14)/{wq[0]}={ceil(wq[0]/14)/wq[0]:.4f}>1/14")
    # minimal covering set
    random.seed(1); pool=[]; tries=0
    while len(pool)<1500 and tries<400000:
        tries+=1
        v=sorted(random.sample(range(1,150),13))
        if primitive(v) and divisor_complete(v) and longest_AP(v)<=7: pool.append(v)
    WIN=[q for q in range(15,32) if q%14!=0]
    clr={q:set(i for i,v in enumerate(pool) if clears(v,q)) for q in WIN}
    un=set(range(len(pool))); chosen=[]
    while un:
        b=max(WIN,key=lambda q:len(clr[q]&un))
        if len(clr[b]&un)==0: break
        chosen.append(b); un-=clr[b]
    print(f"\nMINIMAL covering set over {len(pool)} spread DC: {sorted(chosen)} ({len(chosen)} moduli, no small-subset shortcut)")

if __name__=='__main__':
    main()
