#!/usr/bin/env python3
"""
lrc_halftiling_recursion_invariance_macmini_0620s4.py  (mac-mini-2026-06-20-S4)

(1) even/odd half-size recursions; (2) geometric 7-region inclusion-exclusion;
(3) COMPLEMENT-INVARIANCE of tournament invariants (c3, HP, score-multiset-up-to-reversal)
    => any such invariant is computable on the floor((n-1)^2/4) half-region (2x saving),
    and complement pairs iso classes (SC = self-paired = the diagonal spine).
"""
import itertools, sys
from math import comb
sys.stdout.reconfigure(line_buffering=True)
def h(n): return (n-1)**2//4   # half-tiling size

print("(1) RECURSIONS:")
print(f"   {'n':>3}{'h(n)':>6}{'even 2h(n-1)-h(n-2)':>22}{'odd 2h(n-1)-2h(n-3)+h(n-4)':>28}")
for n in range(2,15):
    ev = 2*h(n-1)-h(n-2) if n>=4 else None
    od = 2*h(n-1)-2*h(n-3)+h(n-4) if n>=5 else None
    ev_ok = (ev==h(n)) if (n>=4 and n%2==0) else '-'
    od_ok = (od==h(n)) if (n>=5 and n%2==1) else '-'
    print(f"   {n:>3}{h(n):>6}{str(ev)+(' OK' if ev_ok==True else (' FAIL' if ev_ok==False else '')):>22}{str(od)+(' OK' if od_ok==True else (' FAIL' if od_ok==False else '')):>28}")

print("\n(2) GEOMETRIC 7-region incl-excl for ODD n: h(n) =? A+B-C+D-E-F+G")
print("    A=B=h(n-1), C=D=h(n-2), E=F=h(n-3), G=h(n-4)")
for n in range(5,15,2):
    A=B=h(n-1); C=D=h(n-2); E=F=h(n-3); G=h(n-4)
    val=A+B-C+D-E-F+G
    print(f"    n={n}: A+B-C+D-E-F+G = {val}  h(n)={h(n)}  {'OK' if val==h(n) else 'FAIL'}")

print("\n(3) COMPLEMENT-INVARIANCE of tournament invariants (reverse all arcs):")
def all_tournaments(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in itertools.product([0,1],repeat=len(pairs)):
        # bit=1 -> i->j, bit=0 -> j->i
        adj=[[0]*n for _ in range(n)]
        for (b,(i,j)) in zip(bits,pairs):
            if b: adj[i][j]=1
            else: adj[j][i]=1
        yield adj
def comp(adj,n):  # reverse all arcs
    return [[adj[j][i] for j in range(n)] for i in range(n)]
def c3(adj,n):
    c=0
    for i,j,k in itertools.combinations(range(n),3):
        # count cyclic triangles among i,j,k
        e=[(i,j),(j,k),(i,k)]
        # cyclic iff out-degrees within triple are all 1
        for a,b,cc in [(i,j,k)]:
            pass
        # directed: check 3-cycle
        def arc(u,v): return adj[u][v]==1
        outs=[sum(adj[x][y] for y in (i,j,k) if y!=x) for x in (i,j,k)]
        if sorted(outs)==[1,1,1]: c+=1
    return c
def scores(adj,n): return tuple(sorted(sum(adj[i]) for i in range(n)))
def ham_paths(adj,n):
    cnt=0
    for p in itertools.permutations(range(n)):
        if all(adj[p[t]][p[t+1]]==1 for t in range(n-1)): cnt+=1
    return cnt
for n in [4,5,6]:
    tot=0; c3_inv=0; hp_inv=0; sc_revinv=0
    for adj in all_tournaments(n):
        tot+=1; co=comp(adj,n)
        if c3(adj,n)==c3(co,n): c3_inv+=1
        if ham_paths(adj,n)==ham_paths(co,n): hp_inv+=1
        # scores under complement: s_i -> (n-1)-s_i, so multiset reverses
        if scores(co,n)==tuple(sorted((n-1)-s for s in scores(adj,n))): sc_revinv+=1
    print(f"   n={n} ({tot} tournaments): c3(T)=c3(T^op) in {c3_inv}/{tot}; HP(T)=HP(T^op) in {hp_inv}/{tot}; "
          f"score-complement law in {sc_revinv}/{tot}")
print("\n   => c3, HP complement-INVARIANT (so summable over the half-region); scores complement to (n-1)-s.")
print("   => SC tournaments (T iso T^op) = self-paired = the diagonal {x+y=n+1} spine.")
