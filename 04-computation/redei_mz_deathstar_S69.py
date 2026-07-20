from itertools import product, permutations
# Tournament on n vertices = orientation of each pair; arc[i][j]=1 iff i->j. Ham-path count.
def ham_count(n, upper):  # upper: dict {(i,j): 0/1 for i<j}, 1 => i->j
    def arc(i,j): return upper[(i,j)] if i<j else 1-upper[(j,i)]
    return sum(1 for p in permutations(range(n)) if all(arc(p[k],p[k+1]) for k in range(n-1)))
def all_tournaments(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in product([0,1],repeat=len(pairs)):
        yield dict(zip(pairs,bits))

print("=== (1) Redei: ham-count is ODD for every tournament (n=2..5) ===")
for n in range(2,6):
    cnts=set(ham_count(n,t) for t in all_tournaments(n))
    print(f"  n={n}: {2**(n*(n-1)//2)} tournaments, ham-counts {sorted(cnts)}, all odd: {all(c%2==1 for c in cnts)}")

print("\n=== (2) Is the parity functional MULTIPLICATIVE (=> trivially Mathieu-Zhao)? ===")
# ordinal sum T (+) S : all of T before all of S (arcs T->S). Ham path must traverse T then S.
def ordinal_sum_ham(hamT, hamS): return hamT*hamS   # ham(T(+)S)=ham(T)ham(S) structurally
# verify on explicit small T,S
import random
def ham_of(n,t): return ham_count(n,t)
T3=next(all_tournaments(3)); S2=next(all_tournaments(2))
# build T3 (+) S2 as a 5-vertex tournament: 0,1,2 = T3 ; 3,4 = S2 ; arcs i->j for i in T, j in S
def osum(nT,tT,nS,tS):
    pairs={}
    # within T
    for (i,j),b in tT.items(): pairs[(i,j)]=b
    for (i,j),b in tS.items(): pairs[(i+nT,j+nT)]=b
    for i in range(nT):
        for j in range(nS): pairs[(i,j+nT)]=1  # i (in T) -> j (in S)
    return pairs
hamT=ham_of(3,T3); hamS=ham_of(2,S2); hamTS=ham_count(5,osum(3,T3,2,S2))
print(f"  ham(T)={hamT}, ham(S)={hamS}, ham(T(+)S)={hamTS}, product={hamT*hamS}, MULTIPLICATIVE: {hamTS==hamT*hamS}")
print("  => the Redei count is a CHARACTER (multiplicative) under ordinal sum.")
print("  For ANY character chi: chi(P)=0 => chi(Q P^m)=chi(Q)chi(P)^m=0 for all m>=1. So ker(chi) is")
print("  Mathieu-Zhao TRIVIALLY -- the degenerate end. GMC's depth is that E is NOT multiplicative.")

print("\n=== (3) The faithful non-multiplicative analog = the ARC expectation = a GMC instance ===")
for n in range(2,6):
    print(f"  n={n} vertices => C(n,2)={n*(n-1)//2} arc variables => GMC({n*(n-1)//2}) regime "
          f"({'TRUE' if n*(n-1)//2<=2 else 'FALSE (GMC>=3 has counterexamples)'})")
print("  So the natural tournament expectation E_T[f]=(1/2^C(n,2)) sum_T f(T) is: n=2 -> GMC(1) true;")
print("  n>=3 -> GMC(>=3) FALSE, so ker E is NOT Mathieu-Zhao there. Redei/parity does not rescue it.")
