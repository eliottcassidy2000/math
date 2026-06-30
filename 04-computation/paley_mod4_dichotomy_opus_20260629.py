"""
How the D_{2p} HP-decomposition DIFFERS by p mod 4 -- and whether it's predictable.

Prediction:
  p=1 mod4: negation is an AUTOMORPHISM => reflection acts as negate-only (no reversal).
            v->-v+a fixes a unique vertex a/2, never a whole path  =>  f=0  =>  a=b=H/(2p).
  p=3 mod4: negation is ANTI-aut => reflection = negate+REVERSE; anti-palindromic HPs exist
            =>  f>0  =>  a-b=f>0.
  Chirality ratio b/(a+b)=(1-pf/H)/2:  =1/2 EXACTLY for p=1mod4;  <1/2 (->1/2) for p=3mod4.
"""
def QR(p): return set((x*x)%p for x in range(1,p))

def H_count(p):
    qr=QR(p); adj=[[ (j-i)%p in qr for j in range(p)] for i in range(p)]
    dp=[[0]*p for _ in range(1<<p)]
    for v in range(p): dp[1<<v][v]=1
    for mask in range(1<<p):
        row=dp[mask]
        for last in range(p):
            c=row[last]
            if not c or not (mask>>last)&1: continue
            al=adj[last]
            for nxt in range(p):
                if (mask>>nxt)&1 or not al[nxt]: continue
                dp[mask|(1<<nxt)][nxt]+=c
    return sum(dp[(1<<p)-1][v] for v in range(p))

def is_HP(P,qr,p): return all((P[i+1]-P[i])%p in qr for i in range(p-1))

def f_count(p):
    """count HPs fixed by THE valid reflection for this p (negate+reverse if 3mod4, negate if 1mod4)."""
    qr=QR(p); m=(p-1)//2
    cnt=0
    if p%4==3:
        # anti-palindromic: v_i=-v_{p-1-i}, middle 0; DFS first half
        pairs=[(a,(p-a)%p) for a in range(1,m+1)]; seq=[]; used=[False]*m
        def arc(i,j): return (j-i)%p in qr
        def dfs(pos):
            nonlocal cnt
            if pos==m:
                full=seq+[0]+[(-x)%p for x in reversed(seq)]
                if all(arc(full[i],full[i+1]) for i in range(p-1)): cnt+=1
                return
            for pi in range(m):
                if used[pi]: continue
                for val in pairs[pi]:
                    if pos>0 and not arc(seq[-1],val): continue
                    used[pi]=True; seq.append(val); dfs(pos+1); seq.pop(); used[pi]=False
        dfs(0)
    else:
        # negate-fixed: P=(-v_0,...,-v_{p-1}) => v_i=-v_i => v_i=0 for all i => impossible
        # (and any reflection v->-v+a: fixes a unique vertex, never a path). f=0 structurally.
        cnt=0
    return cnt

print(f"{'p':>3} {'p%4':>4} {'neg=':>8} {'H':>10} {'f':>7} {'a=triv':>9} {'b=sign':>9} {'c':>9} {'b/(a+b)':>9}")
for p in [3,5,7,11,13,19]:
    H=H_count(p); f=f_count(p)
    a=(H+p*f)//(2*p) if (H+p*f)%(2*p)==0 else (H+p*f)/(2*p)
    b=(H-p*f)//(2*p) if (H-p*f)%(2*p)==0 else (H-p*f)/(2*p)
    c=H//p if H%p==0 else H/p
    negkind = "anti-aut" if p%4==3 else "AUT"
    intok = isinstance(a,int) and isinstance(b,int) and isinstance(c,int)
    print(f"{p:>3} {p%4:>4} {negkind:>8} {H:>10} {f:>7} {str(a):>9} {str(b):>9} {str(c):>9} {(1-p*f/H)/2:>9.4f}  int={intok}")

print("\nDICHOTOMY: p=1mod4 -> f=0, a=b, b/(a+b)=1/2 exactly (reflection has no fixed path).")
print("           p=3mod4 -> f>0, a>b by exactly f, b/(a+b)<1/2 ->1/2 (anti-palindromic HPs).")
print("The discriminant is the REVERSAL forced by the anti-automorphism -- the same p mod 4")
print("that gives Borsuk-Ulam(odd) vs Brouwer(even), imaginary vs real, Q(sqrt-p) vs Q(sqrt p).")
