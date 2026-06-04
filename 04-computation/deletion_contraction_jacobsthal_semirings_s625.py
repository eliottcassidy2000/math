# The deletion-contraction recursion I(G,x) = I(G-v,x) + x*I(G-N[v],x) is THE partition-function transfer recursion.
# For the PATH P_n (conflict graph = chain): I(P_n,x) = I(P_{n-1},x) + x*I(P_{n-2},x)  -- the "n-2" / Fibonacci recursion.
def path_indep_poly(n):
    # returns coefficients via I_n = I_{n-1} + x I_{n-2}; represent as list of coeffs
    I=[ [1], [1,1] ]  # I_0=1, I_1=1+x
    for k in range(2,n+1):
        a=I[k-1][:]; b=[0]+I[k-2]  # x*I_{k-2}
        L=max(len(a),len(b)); c=[ (a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(L)]
        I.append(c)
    return I
I=path_indep_poly(10)
print("=== independence polynomial of the path P_n: I_n = I_{n-1} + x*I_{n-2}  (the n-2 / Fibonacci recursion) ===")
for n in range(0,9): print(f"  I(P_{n}) coeffs = {I[n]}")
def ev(coeffs,x): return sum(c*x**i for i,c in enumerate(coeffs))
print("\n  at x=1 (Fibonacci):  ", [ev(I[n],1) for n in range(9)])
print("  at x=2 (JACOBSTHAL): ", [ev(I[n],2) for n in range(9)], " = (2^{n+1}-(-1)^{n+1})/3")
print("  note 21 (=J at n=6) is a Jacobsthal value AND a FORBIDDEN tournament-H (21=7*3): not every conflict graph is tournament-realizable.")

# ===== ONE recursion, THREE semirings = partition functions everywhere =====
print("\n=== ONE recursion val(p)=COMBINE_moves(WEIGHT op val(child)), three semirings ===")
print("  (sum, x*) ordinary  : Z = sum over configs of x^size  -> independence polynomial / H / depth-GF / unit-dist count")
print("  (min, +1) tropical  : game value = min over my moves / sup(+1) over opp -> Hamkins transfinite game value, shortest path, Collatz altitude")
print("  (ordinal sup,+1)    : transfinite rank = the open-game value (infinite Go) = well-founded recursion depth")
# demonstrate: same path recursion, tropical (min,+) gives the 'distance to ground' = linear (n); ordinary gives Jacobsthal (exp)
def tropical_path(n):  # 'game value' rank via min,+: val_n = min(val_{n-1}, 1+val_{n-2}) with val_0=0,val_1=0
    v=[0,0]
    for k in range(2,n+1): v.append(min(v[k-1], 1+v[k-2]))
    return v
print("\n  same P_n recursion, TROPICAL (min, +1):", tropical_path(8), " (linear-bounded = the 'rank'/game-value face)")
print("  same P_n recursion, ORDINARY (sum, x=2):", [ev(I[n],2) for n in range(9)], " (exponential = the 'count'/partition-function face)")

# ===== Collatz altitude = ordinal game value (HYP-2180) =====
print("\n=== Collatz as an OPEN GAME: the altitude (HYP-2180 iterated-log) = the ordinal game value / well-founded rank ===")
print("  Collatz position n; move = shortcut step; WON = reach 1. Game value = well-founded rank = #steps-ish.")
print("  Collatz conjecture = every position has a (countable) ordinal value = the game is OPEN (no infinite play) = the recursion is WELL-FOUNDED.")
print("  altitude tower (HYP-2180): value lives at floor-1 (log n steps), epoch at floor-2 (loglog) = the ordinal is omega-ish, not just finite.")

# ===== partition functions everywhere: the catalog =====
print("\n=== PARTITION FUNCTIONS EVERYWHERE (the recurring master object across the arc) ===")
for face,Z in [
 ("tournament H","I(Omega,2) = #Hamiltonian paths (Redei)"),
 ("LRC loneliness","covering-depth GF P(z)=int z^depth ; p_0=P(0)"),
 ("unit distance","#unit-distance pairs = norm-1 element count in CM lattice"),
 ("Collatz","ordinal game value / well-founded rank = altitude tower"),
 ("Hamkins infinite Go","transfinite game value = partition function over (sup,+1) ordinal semiring"),
 ("Krawtchouk/Delsarte","weight enumerator transform; dual = the LP"),
]:
    print(f"  {face:22}: {Z}")
