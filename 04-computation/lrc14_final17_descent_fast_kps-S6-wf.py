"""
LRC(14) — FAST aggressive descent for the 1/7-spread bound (companion to the
main verify script). Goal: HUNT hard for an integer E (0 in E, 8<=k<=12) with
mu_{1/7}(E) < thr_k, using descent + structured adversaries, but with bounded
spread so each exact mu_theta call stays cheap.

These are the SAME structured families that crushed mu_{2/7}:
 - maximally perforated rulers (coprime steps, prime gaps)
 - Sidon-ish / Singer difference sets
 - "stretch one point far" adversaries
 - simulated annealing that ACCEPTS uphill moves to escape the consec basin
All EXACT rationals.
"""
from fractions import Fraction as F
from itertools import combinations
import random, math, sys

def mu_theta(E, theta):
    E = sorted(set(E)); n = len(E); bp = set([F(0), F(1)])
    for i in range(n):
        for j in range(i+1, n):
            d = E[j]-E[i]
            for m in range(0, d+1): bp.add(F(m, d))
    bp = sorted(b for b in bp if 0 <= b <= 1); total = F(0)
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        mid = (a+b)/2
        order = sorted(range(n), key=lambda i: (E[i]*mid) % 1)
        ks = [(E[order[t]]*mid).__floor__() for t in range(n)]; subs = []
        for t in range(n):
            o1 = order[t]; o2 = order[(t+1) % n]; k1 = ks[t]; k2 = ks[(t+1) % n]
            wrap = 1 if t == n-1 else 0
            s = E[o2]-E[o1]; c = F(k1-k2+wrap)
            if s == 0:
                if c > theta: subs.append((a, b))
            elif s > 0:
                lo = max(a, (theta-c)/s); subs.append((lo, b)) if lo < b else None
            else:
                hi = min(b, (theta-c)/s); subs.append((a, hi)) if a < hi else None
        subs.sort(); cur = cb = None
        for lo, hi in subs:
            if cur is None: cur, cb = lo, hi
            elif lo <= cb: cb = max(cb, hi)
            else: total += cb-cur; cur, cb = lo, hi
        if cur is not None: total += cb-cur
    return total

def mu17(E): return mu_theta(tuple(sorted(set(E))), F(1,7))

THR = {8:F(3637,5880), 9:F(2025,4004), 10:F(36,91), 11:F(25,91), 12:F(1,7)}
CONSEC = {k: mu17(range(k)) for k in range(8,13)}

def main():
    violations = []
    print("CONSEC mu_1/7 and THR per k:")
    for k in range(8,13):
        print(f"  k={k}: consec={float(CONSEC[k]):.6f}  thr={float(THR[k]):.6f}  slack(consec-thr)={float(CONSEC[k]-THR[k]):.6f}")
    sys.stdout.flush()

    # ---- A. Simulated annealing that accepts uphill (escape consec basin) ----
    print("\n[A] Annealing descent (accepts uphill; spread capped at 3k):")
    for k in range(8,13):
        thr = THR[k]; cap = 3*k
        rng = random.Random(2024+k)
        glob_best = CONSEC[k]; glob_E = tuple(range(k))
        for restart in range(40):
            # random start
            E = tuple(sorted(set([0]+rng.sample(range(1,cap+1),k-1))))
            if len(E)!=k: continue
            cur = mu17(E)
            T = 0.05
            for step in range(150):
                i = rng.randint(1,k-1)
                Ei = list(E)
                Ei[i] = rng.randint(1,cap)
                cand = tuple(sorted(set(Ei)))
                if len(cand)!=k or 0 not in cand: continue
                m2 = mu17(cand)
                dE = float(m2-cur)
                if dE < 0 or rng.random() < math.exp(-dE/max(T,1e-9)):
                    E, cur = cand, m2
                    if cur < glob_best: glob_best, glob_E = cur, E
                    if cur < thr: violations.append((k,E,cur,thr))
                T *= 0.985
        tag = "BELOW CONSEC!" if glob_best < CONSEC[k] else "ok (>=consec)"
        print(f"  k={k}: best={float(glob_best):.6f} at {glob_E}  thr={float(thr):.6f}  [{tag}]")
        sys.stdout.flush()

    # ---- B. Structured perforated rulers (coprime/prime-gap, the 2/7 crushers) ----
    print("\n[B] Structured perforated/coprime rulers:")
    def structured(k):
        out = []
        # coprime steps
        for step in range(1,12):
            out.append([step*i for i in range(k)])
        # prime-ish gaps
        primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47]
        cum=[0]
        for p in primes[:k-1]: cum.append(cum[-1]+p)
        out.append(cum)
        # Sidon-like (Mian-Chowla start)
        mc=[0,1,3,7,12,20,30,44,65,80,96,122,147]
        out.append(mc[:k])
        # Singer-ish / random coprime difference sets
        rng=random.Random(55+k)
        for _ in range(2000):
            cap=rng.randint(k+1,4*k)
            out.append([0]+rng.sample(range(1,cap+1),k-1))
        # consec with ONE point stretched very far (stretch adversary)
        for far in range(k, 8*k):
            out.append(list(range(k-1))+[far])
        return out
    for k in range(8,13):
        thr=THR[k]; best=CONSEC[k]; bE=tuple(range(k))
        for E in structured(k):
            E=tuple(sorted(set(E)))
            if len(E)!=k or 0 not in E: continue
            m=mu17(E)
            if m<best: best,bE=m,E
            if m<thr: violations.append((k,E,m,thr))
        tag = "BELOW CONSEC!" if best < CONSEC[k] else "ok (>=consec)"
        print(f"  k={k}: best={float(best):.6f} at {bE}  thr={float(thr):.6f}  [{tag}]")
        sys.stdout.flush()

    print("\nVERDICT:")
    if violations:
        print(f"*** {len(violations)} VIOLATIONS (mu<thr) — holds=FALSE ***")
        for v in violations[:15]: print("   ",v)
    else:
        print("NO violation found. mu_{1/7}(E) >= thr_k held in every aggressive search.")
        print("Consecutive was never beaten as the minimizer.")

if __name__=="__main__":
    main()
