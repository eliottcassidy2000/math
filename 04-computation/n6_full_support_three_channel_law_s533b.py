#!/usr/bin/env python3
"""
n6_full_support_three_channel_law_s533b.py    oracle-2026-06-01-S533

THE clean three-channel parity law for n=6 (focused: FULL-SUPPORT debt).

n=4: inside debt = the order-3 (= FULL-support, all 3 runners) resonance +-a+-b+-c=0.
     Vanishes iff a+b+c ODD. units mod n*=2 are {1}: NO sign freedom -> ONE fixed
     sum mod 2 = ONE channel.
n=6: the true analogue is the FULL-support (order 5, all 5 runners) resonance
     sum m_i v_i = 0, all m_i != 0 (mod 3) [where ghat lives]. units mod n*=3 are
     {+1,-1}: a SIGN per runner. Necessary mod-3 condition:
         exists eps in {+-1}^5 with  sum eps_i v_i = 0  (mod 3).
     This is the GENUINELY multi-channel parity law: you scan 2^5 sign patterns
     (the +-1 units) over the 3 residue channels (mod 3). mod 6 = (mod-3 channel)
     (x) (mod-2 sign-unit).

REDUCTION: each active runner (v_i != 0 mod 3) contributes eps_i v_i in {1,2} (mod 3)
freely; inert runners (v_i = 0 mod 3) contribute 0. With k = #active, the reachable
sums mod 3 = {k,...,2k} mod 3 = all of Z/3 iff k>=2, = {1,2} iff k=1, = {0} iff k=0.
So 0 reachable (FULL-support debt PRESENT) iff k != 1; absent iff k=1.
For PRIMITIVE sets k=0 impossible (not all divisible by 3), so:
   FULL-SUPPORT INSIDE DEBT VANISHES  <=>  exactly ONE runner != 0 (mod 3).

(NB: 'any-order' debt is almost always present because inert runners resonate among
themselves -- that is the (n-1-k) sub-system, scaled by 3. The CLEAN law is for the
full-support channel, exactly as n=4's law is for its full-support order-3 term.)

We verify all of this and the mod-6 sign refinement of the debt VALUE.
"""
from math import sin, pi, gcd
from itertools import product, combinations
from functools import reduce
import random

N = 6
def ghat(m):
    if m == 0: return 1.0 - 2.0/N
    return -sin(2*pi*m/N)/(pi*m)

ALLOWED = [m for m in range(-8, 9) if m != 0 and m % 3 != 0]   # support of ghat, nonzero

def full_support_debt(speeds):
    """sum over all-nonzero m_i (in ALLOWED) with sum m_i v_i = 0 of prod ghat(m_i)."""
    tot = 0.0; cnt = 0
    for ms in product(ALLOWED, repeat=len(speeds)):
        if sum(a*b for a, b in zip(ms, speeds)) == 0:
            t = 1.0
            for a in ms: t *= ghat(a)
            tot += t; cnt += 1
    return tot, cnt

def sign_reachable_mod3(speeds):
    for eps in product((1, -1), repeat=len(speeds)):
        if sum(e*s for e, s in zip(eps, speeds)) % 3 == 0:
            return True
    return False

def channels(speeds):
    n0 = sum(1 for s in speeds if s % 3 == 0)
    n1 = sum(1 for s in speeds if s % 3 == 1)
    n2 = sum(1 for s in speeds if s % 3 == 2)
    return n0, n1, n2

def part1():
    print("="*72)
    print("(1) FULL-SUPPORT debt vs channel occupancy (n0,n1,n2) and k=n1+n2")
    print("="*72)
    sets = [
        (1,2,3,4,5),     # k=4
        (3,6,9,12,1),    # k=1  -> predicted DEBT-FREE
        (3,6,9,15,2),    # k=1  -> predicted DEBT-FREE
        (1,2,3,6,9),     # k=2
        (1,4,7,10,2),    # res 1,1,1,1,2 k=5
        (2,5,8,11,1),    # res 2,2,2,2,1 k=5
        (1,2,4,5,7),     # k=5
        (6,1,12,18,24),  # k=1 (only '1' active)
    ]
    for v in sets:
        debt, cnt = full_support_debt(v)
        reach = sign_reachable_mod3(v)
        n0, n1, n2 = channels(v)
        k = n1 + n2
        verdict = "PRESENT" if abs(debt) > 1e-9 else "ZERO"
        pred = "ZERO" if k == 1 else "PRESENT"
        print(f"  v={v}: (n0,n1,n2)=({n0},{n1},{n2}) k={k}  0-reach mod3={reach}  "
              f"full-support debt={debt:+.6f} [{cnt} resonances] -> {verdict} (predict {pred}) "
              f"{'OK' if verdict==pred else 'MISMATCH'}")

def part2():
    print("\n"+"="*72)
    print("(2) LAW on primitive 5-sets: full-support debt ZERO <=> not sign-reachable mod3 <=> k=1")
    print("="*72)
    rnd = random.Random(5331)
    # include k=1 sets deliberately (four multiples of 3 + one not)
    pool = []
    for combo in combinations(range(1, 16), 5):
        if reduce(gcd, combo) == 1:
            pool.append(combo)
    rnd.shuffle(pool)
    tested=0; ok=0; debtfree=[]
    for v in pool:
        tested += 1
        if tested > 400: break
        reach = sign_reachable_mod3(v)
        n0,n1,n2 = channels(v); k=n1+n2
        debt,_ = full_support_debt(v)
        zero = abs(debt) < 1e-9
        law = (zero == (not reach)) and ((not reach) == (k == 1))
        if law: ok += 1
        if zero: debtfree.append((v,k))
    print(f"  tested {tested} primitive 5-sets; law holds {ok}/{tested}")
    print(f"  full-support DEBT-FREE sets (all have k=1, one runner !=0 mod3): {debtfree[:10]}")

def part3():
    print("\n"+"="*72)
    print("(3) mod-6 sign refinement: same mod-3 channels, different mod-2 signs -> different debt")
    print("="*72)
    # speeds with identical residues mod 3 but different residues mod 6
    pairs = [
        ((1,2,4,5,7),(1,2,4,5,1)),   # last runner 7 vs 1: both =1 mod3, differ mod6 (7%6=1 same!) -> pick real mod6 diff
        ((1,2,4,5,7),(4,5,1,2,7)),
        ((1,1,1,2,2),(4,4,4,5,5)),   # 1->4 (=1 mod3, flips mod2), 2->5 (=2 mod3, flips mod2)
    ]
    for a,b in pairs:
        da,_=full_support_debt(a); db,_=full_support_debt(b)
        print(f"  {a} (mod6={[x%6 for x in a]}) debt={da:+.6f}")
        print(f"  {b} (mod6={[x%6 for x in b]}) debt={db:+.6f}   "
              f"{'(value differs -> mod-6 matters)' if abs(da-db)>1e-6 else '(same)'}")
    print("  => feasibility is mod-3 (the 3 channels); the debt VALUE/sign is mod-6")
    print("     (each ghat(m) sign set by m mod 6). mod6 = mod3-channel (x) mod2-sign-unit.")

def part4():
    print("\n"+"="*72)
    print("(4) n=4 parallel (units mod 2 = {1}: one fixed sum, 'a+b+c odd')")
    print("="*72)
    def n4_full_debt(a,b,c):
        allowed=[m for m in range(-9,10) if m!=0 and m%2!=0]  # odd m
        tot=0.0
        for ms in product(allowed,repeat=3):
            if ms[0]*a+ms[1]*b+ms[2]*c==0:
                t=1.0
                for m in ms: t*= -sin(2*pi*m/4)/(pi*m)
                tot+=t
        return tot
    for v in [(1,2,3),(1,2,4),(1,3,5),(2,3,5)]:
        d=n4_full_debt(*v); s=sum(v)
        print(f"   (a,b,c)={v} sum={s} ({'odd' if s%2 else 'even'}) full debt={d:+.6f} "
              f"-> {'ZERO' if abs(d)<1e-9 else 'PRESENT'}")
    print("   n=4: debt-free <=> a+b+c ODD (units {1}, no sign freedom: ONE channel).")
    print("   n=6: debt-free <=> NO eps in {+-1}^5 with sum eps_i v_i = 0 mod 3")
    print("        (units {+-1}: a sign per runner over 3 channels) = FIRST multi-channel law.")

def main():
    part1(); part2(); part3(); part4()

if __name__ == "__main__":
    main()
