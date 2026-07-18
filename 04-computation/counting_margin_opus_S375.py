from math import gcd
def W(q): return set([0]+[r%q for j in range(1,(q-1)//14+1) for r in (j,-j)])
def phi(n): return sum(1 for k in range(1,n+1) if gcd(k,n)==1)
def k_q(q): return sum(1 for w in W(q) if w!=0 and gcd(w,q)==1)
print("(7) THE COUNTING THEOREM DOES FIRE -- just not at 13 speeds.")
print("    THEOREM: if s speeds, none divisible by q, and s*k_q < phi(q),")
print("    then some p/q is lonely (bad numerators < phi(q), so one survives).")
print("    best modulus = the one maximising phi(q)/k_q:")
best=[]
for q in range(15,200):
    k=k_q(q)
    if k: best.append((phi(q)/k, q, phi(q), k))
best.sort(reverse=True)
for r,q,ph,k in best[:8]:
    print(f"      q={q:4d}  phi={ph:4d}  k_q={k:3d}  phi/k = {r:6.2f}  -> fires for s <= {-(-ph//k)-1 if ph%k==0 else ph//k}")
top=best[0]
smax=(top[2]-1)//top[3]
print()
print(f"    MAXIMUM s for which the counting proof works: s = {smax}  (at q = {top[1]})")
print(f"    LRC(14) needs s = 13.  The proof covers s <= {smax}, missing by {13-smax}.")
print(f"    i.e. this argument closes the analogue for up to {smax+1} runners and")
print(f"    falls exactly {13-smax} speeds short of the open case.")
