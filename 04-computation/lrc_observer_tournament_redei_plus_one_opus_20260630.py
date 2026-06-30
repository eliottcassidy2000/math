"""
The OBSERVER on the tournament side. The marked vertex v = the observer; it INSERTS itself into each
Ham path P' of T-v. Its view of P' = the signature s (v->p or p->v). inshat(v,P')=1+2*#TypeII = ODD
(Redei), TypeII = the observer's directed 3-cycles. H(T)-H(T-v)=2*sum mu(C) over odd cycles THROUGH v
(even => parity preserved => Redei). The transitive observer (no odd cycles) = the baseline H=1.
"""
from itertools import permutations
def ham_paths(verts, beats):
    cnt=0
    for p in permutations(verts):
        if all(beats(p[i],p[i+1]) for i in range(len(p)-1)): cnt+=1
    return cnt
def signature_and_typeII(v, Ppath, beats):
    s=[1 if beats(v,x) else 0 for x in Ppath]
    typeII=sum(1 for j in range(len(s)-1) if s[j]==1 and s[j+1]==0)
    typeI =sum(1 for j in range(len(s)-1) if s[j]==0 and s[j+1]==1)
    b = s[0] + (1 - s[-1])
    inshat = b + typeI + typeII
    return s, typeII, inshat
# example tournaments on n=5 (adjacency as beats-sets)
def make(adjpairs,n):
    beatset={(i,j) for (i,j) in adjpairs}
    return lambda a,b:(a,b) in beatset
# transitive T5: i->j iff i<j
nverts=5; V=list(range(nverts))
trans=lambda a,b: a<b
# a 3-cycle-rich tournament: rotational (i -> i+1,i+2 mod 5) = the quadratic-residue/Paley-like
rot=lambda a,b: ((b-a)%5) in (1,2)
for name,beats in [("transitive T5",trans),("rotational T5 (cyclic)",rot)]:
    H=ham_paths(V,beats)
    print(f"\n{name}: H(T)={H}  (Redei: odd = {H%2==1})")
    v=0  # the observer
    rest=[x for x in V if x!=v]
    Hmv=ham_paths(rest,beats)
    # sum over P' of insact via signed inshat (observer self-insertion)
    tot_inshat=0; details=[]
    for p in permutations(rest):
        if all(beats(p[i],p[i+1]) for i in range(len(p)-1)):
            s,t2,ins=signature_and_typeII(v,p,beats)
            tot_inshat+=ins; details.append((p,s,t2,ins))
    print(f"  observer v={v}: H(T-v)={Hmv}; paths P' of T-v: {len(details)}")
    print(f"  each inshat=1+2*#TypeII (odd): {[ (d[2],d[3]) for d in details][:6]}  [(#TypeII, inshat)]")
    print(f"  all inshat odd: {all(d[3]%2==1 for d in details)}; sum inshat over P' = {tot_inshat} = H(T)? {tot_inshat==H}")
    print(f"  H(T)-H(T-v) = {H-Hmv} (EVEN = {((H-Hmv)%2==0)}): the observer's odd-cycle correction 2*sum mu")
print("\nTHE OBSERVER READING (Redei/OCF):")
print(" * the observer v inserts into each path P'; its count inshat = 1 (BASELINE) + 2*#TypeII (its 3-cycles).")
print(" * the '+1' is the observer's irreducible self-insertion; the 2*#TypeII is its odd-cycle interaction.")
print(" * deleting v changes H by an EVEN 2*sum mu(C) over odd cycles THROUGH v => parity preserved => Redei.")
print(" * TRANSITIVE observer = NO odd cycles => pure baseline H=1 (the metagraph's identity/origin class).")
print()
print("UNIFICATION with the LRC observer -- the same '+1' BASELINE:")
print("  LRC observer:    escape = 1/(n(n-1)+1) form; the +1 in Phi_6=n(n-1)+1 = the Farey hair (clears 1/n).")
print("  tournament obs:  inshat = 1 + 2*#TypeII; the +1 = the irreducible self-insertion (clears parity, Redei).")
print("  => both observers carry an IRREDUCIBLE +1 (baseline) + an even correction (2 x odd-cycle/blocking).")
print("     LRC: +1 makes the escape POSITIVE (>1/n). Tournament: +1 makes the count ODD (Redei). Same atom.")
