"""Affine bank on root-response triples (h, e, B) = (h(T), e_t(T), B_t(T)).
Maps realised by explicit constructions (root t kept unless stated):
  S: add a sink z (everyone -> z):            (h, e, B) -> (h, 0, 2B)
  X: add x with sole exit x -> t, re-root at x: (h, e, B) -> (h - e + B, h - e, B + h)
  R: re-root at a new sink:                    (h, e, B) -> (h, h, h)          [then X gives (h, 0, 2h)]
Start: the single vertex (1, 1, 1).  Values: h itself and the two-step h + 2B (= S then X).
BFS to a given depth; report odd numbers below LIMIT not reached as h."""
import sys
depth=int(sys.argv[1]); LIMIT=int(sys.argv[2]); maps=sys.argv[3] if len(sys.argv)>3 else "SX"
start=(1,1,1)
seen={start}; frontier=[start]
for d in range(depth):
    nxt=[]
    for (h,e,B) in frontier:
        cands=[]
        if 'S' in maps: cands.append((h,0,2*B))
        if 'X' in maps: cands.append((h-e+B,h-e,B+h))
        if 'R' in maps: cands.append((h,h,h))
        if 'P' in maps: cands.append((3*h,3*e,3*B))      # prefix join by C3: C3 -> T (h, e_t, B_t all x3)
        if 'Q' in maps: cands.append((3*h,0,12*B))       # suffix join by C3: T -> C3 (h x3, e_t = 0, B_t x Q(C3) = 12)
        if 'T' in maps: cands.append((h,0,2*B))          # transitive suffix of order 1 = S (kept for clarity)
        for c in cands:
            if c[0] > 4*LIMIT: continue   # prune: h only grows under X (h' = h + (B - e) >= h since B >= e? B >= h(T-t) >= e); S, R keep h
            if c not in seen:
                seen.add(c); nxt.append(c)
    frontier=nxt
hs={h for (h,e,B) in seen}
odd_missing=[v for v in range(1,LIMIT,2) if v not in hs]
print(f"maps {maps}, depth {depth}: {len(seen)} states; h-values below {LIMIT}: {sum(1 for v in range(1,LIMIT,2) if v in hs)} of {LIMIT//2}; first missing {odd_missing[:25]}; count missing {len(odd_missing)}")
# largest X with all odd in [a, X] covered
def interval(S,a):
    X=a
    while X+2 in S: X+=2
    return X
for a in (1,69,307):
    print(f"   from {a}: every odd in [{a}, {interval(hs,a)}] is an h-value")
