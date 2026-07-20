#!/usr/bin/env python3
"""kind-pasteur-S128c91 -- HYP-8000 triage: randomized-greedy ALIVE detection.
A found <=13 cover is a certificate that p is ALIVE (c(p) <= 13).  Greedy with
randomized top-t choice + restarts finds covers fast when they are plentiful;
failures are INCONCLUSIVE (the exhaustive C search owns the extinct side).
Output: the alive-certified set and the frontier candidates (greedy-fails)."""
import sys, random, time
random.seed(91)

def primes_in(lo, hi):
    sieve = bytearray([1])*(hi+1); sieve[0:2] = b'\x00\x00'
    for i in range(2, int(hi**0.5)+1):
        if sieve[i]: sieve[i*i::i] = b'\x00'*len(sieve[i*i::i])
    return [i for i in range(lo, hi+1) if sieve[i]]

def triage(p, n=14, capsize=13, restarts=400, topt=3):
    h = (p-1)//2; dk = p//n
    def fold(x):
        x %= p
        return x if x <= h else p - x
    maskA = [0]*(h+1)
    for w in range(1, h+1):
        m = 0
        for j in range(1, dk+1): m |= 1 << (fold(j*w) - 1)
        maskA[w] = m
    FULL = (1 << h) - 1
    ws = list(range(1, h+1))
    for _ in range(restarts):
        U = FULL; chosen = 0
        for step in range(capsize):
            # rank remaining w by coverage of U; random among top-t
            best = []
            for w in ws:
                c = (maskA[w] & U).bit_count()
                if c:
                    best.append((c, w))
            if not best: break
            best.sort(reverse=True)
            c, w = random.choice(best[:topt])
            U &= ~maskA[w]; chosen += 1
            if U == 0: return True, chosen
        if U == 0: return True, chosen
    return False, None

def main():
    lo, hi = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) >= 3 else (321, 1200)
    n = int(sys.argv[3]) if len(sys.argv) >= 4 else 14
    cap = 13 if n == 14 else 12
    print(f"== greedy triage, n={n} cap={cap}, p in [{lo},{hi}] ==", flush=True)
    frontier = []
    for p in primes_in(lo, hi):
        if p % (7 if n == 14 else 13) == 0: continue
        t0 = time.time()
        ok, sz = triage(p, n=n, capsize=cap)
        if ok:
            print(f"  p={p} ALIVE (greedy cover, size {sz}) [{time.time()-t0:.1f}s]", flush=True)
        else:
            frontier.append(p)
            print(f"  p={p} greedy-FAIL (frontier candidate) [{time.time()-t0:.1f}s]", flush=True)
    print(f"== frontier candidates (need exhaustive search): {frontier}", flush=True)

if __name__ == "__main__":
    main()
