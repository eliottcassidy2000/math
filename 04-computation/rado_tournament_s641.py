#!/usr/bin/env python3
"""
S641 — the Rado graph as a tournament = the universal homogeneous (random) tournament; Paley
tournaments are its explicit finite approximants. "A loop of the input swaps the two outputs":
multiplication by a quadratic NON-RESIDUE g swaps T <-> T^op (self-complementation), and g generates
the Z/2 of the QR double cover (the monodromy / deck transformation of sqrt). The convergence rate to
the Rado limit is the Weil/Gauss-sum bound sqrt(p) = the character-ratio spectrum (S638).
We verify: (1) the non-residue swap T<->T^op (the 'loop swaps two outputs'); (2) k-existential-closure
(the extension property) of Paley tournaments and its sqrt(p) threshold (-> Rado in the limit).
"""
import itertools, math

def legendre(a,p):
    a%=p
    return 0 if a==0 else (1 if pow(a,(p-1)//2,p)==1 else -1)
def isprime(n):
    if n<2: return False
    i=2
    while i*i<=n:
        if n%i==0: return False
        i+=1
    return True

def paley_adj(p):   # i->j iff (j-i) is a QR mod p   (p≡3 mod4 => tournament)
    return [[legendre((j-i)%p,p)==1 for j in range(p)] for i in range(p)]

print("(1) THE LOOP SWAPS THE TWO OUTPUTS: x|->g*x (g a non-residue) maps T_p -> T_p^op")
for p in [7,11,19,23]:
    if p%4!=3: continue
    g=next(a for a in range(2,p) if legendre(a,p)==-1)   # a non-residue = the 'loop'
    adj=paley_adj(p)
    # check: i->j in T  iff  g*i <- g*j  (arc reversed) for all i,j
    ok=all(adj[i][j]==adj[(g*j)%p][(g*i)%p] for i in range(p) for j in range(p) if i!=j)
    print(f"   p={p}: non-residue g={g}; x->g*x sends T->T^op (arc-reversal/self-complement)? {ok}; g^2 a QR? {legendre(g*g,p)==1} (g generates the Z/2 QR cover)")

print("\n(2) EXTENSION PROPERTY (k-existential closure): the finite Rado/random structure.")
print("    k-EC = for every disjoint U,V with |U|+|V|<=k, EXISTS x with x->U and V->x.")
def min_extension_count(adj,p,s_u,s_v):
    """min over disjoint U(|U|=s_u),V(|V|=s_v) of #{x notin U,V: x beats all U, loses to all V}."""
    best=None
    verts=range(p)
    for U in itertools.combinations(verts,s_u):
        Us=set(U)
        for V in itertools.combinations([v for v in verts if v not in Us],s_v):
            Vs=set(V)
            cnt=sum(1 for x in verts if x not in Us and x not in Vs
                    and all(adj[x][u] for u in U) and all(adj[v][x] for v in V))
            if best is None or cnt<best: best=cnt
    return best
for p in [7,11,19,23,31,43]:
    if p%4!=3: continue
    adj=paley_adj(p)
    row=[]
    for s in (1,2,3):      # patterns with |U|+|V|=s (worst split)
        worst=min(min_extension_count(adj,p,su,s-su) for su in range(s+1))
        row.append((s,worst))
    kec=max(s for s,w in row if w>0) if any(w>0 for s,w in row) else 0
    pred=[(s, round(p/2**s,1), round((s)*math.sqrt(p),1)) for s,_ in row]   # Weil: ~p/2^s ± s√p
    print(f"   p={p}: min #extenders by |U|+|V|: {row}  => {kec}-EC ;  Weil pred p/2^s (±s√p): {pred}")

print("\n  => Paley tournaments are k-EC for p >> 4^k (Weil/√p); as p->inf they converge to the unique")
print("     countable homogeneous RANDOM TOURNAMENT (the Rado graph's tournament analogue), which is")
print("     k-EC for ALL k. The √p error = the Gauss sum = the character-ratio spectrum (S638) = the")
print("     monodromy: the non-residue 'loop' that swaps T<->T^op is the Z/2 whose Weil bound is the")
print("     rate of approach to the universal limit.")
