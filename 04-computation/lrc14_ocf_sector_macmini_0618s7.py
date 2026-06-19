#!/usr/bin/env python3
"""
lrc14_ocf_sector_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

THE OCF / ODD-CYCLE ANGLE on the LRC residual, plus the AP=transitive extremality test.

Established (this session, exact, 0 failures):
  R2: maxgap(x)>1/2  <=>  T(x) has a Condorcet winner (score k-1) and loser (score 0)
      <=>  T(x) is NOT strongly connected.
  R3: maxgap(x)>1/7  <=>  T(x) has a SCALE-1/7 local sink (empty 1/7-arc after some pt).
  R4: c3(T(x)) = C(k,3) - sum_i C(s_i,2).
So mu_{1/7}(E) = P_x[ scale-1/7 local sink ].   And measN = P[1/7-net] = P[NO 1/7-sink].

NOW develop:
 (A) STRONG-CONNECTIVITY MEASURE.  measSC = P[T(x) strongly connected] = P[maxgap<1/2].
     A round tournament is strongly connected  <=>  it has a Hamiltonian CYCLE  <=>  every
     open semicircle is nonempty.  For LRC: maxgap<1/2 always (k>=3 distinct phases can't
     all fit a semicircle for generic x) — but the FRACTION matters.  Compute measSC.
 (B) OCF on the round tournament.  H(T) = #directed Ham PATHS (Redei: odd).  For a round
     tournament, count odd cycles and test the OCF gradient.  We compute E_x[H], E_x[c3],
     E_x[#5-cycles] and ask: is meas(S7) a LINEAR functional of these cycle-moments?
     meas(S7) is INCLUSION-EXCLUSION over sector-occupancy; cycle counts are I-E over
     winding.  Test the linear fit meas(S7) ~ a + b*E[c3] + ... across shapes.
 (C) AP = TRANSITIVE EXTREMALITY.  Claim (the angle's hypothesis): among all k-shapes the
     AP {0..k-1} produces the tournament ensemble {T(x)} that is "MOST TRANSITIVE on
     average" (max E[transitive-triples] = min E[c3]) AND simultaneously fills sectors
     most (max meas(S7)).  Test:  does AP minimise E[c3]  (= maximise E[#transitive
     triples] = sum_i C(s_i,2))?  And does AP maximise measN / meas(S7)?  This would say:
     "the EXTREMAL LRC cluster is the one whose winding tournament is closest to transitive"
     — a clean parity-side statement of AP-extremality.
 (D) THE WALSH / SECTOR-FOURIER LINK.  meas(S7) main term M7(k) = E_x[indicator all 7
     sectors hit] for INDEPENDENT phases.  The 7-vanishing (THM-503: hat a_T(7m)=0) is the
     statement that the sector-occupancy is BLIND to the 7-divisible winding modes — exactly
     the modes that would create a 7-fold symmetric (odd) winding.  Verify the sector
     character orthogonality numerically.

Exact Fraction.  H/c3 via permutation enumeration (k<=8).
"""
import sys, itertools
from fractions import Fraction as F
from math import comb
sys.stdout.reconfigure(line_buffering=True)

def round_tournament(E,x):
    n=len(E); A=[[0]*n for _ in range(n)]; tf=True
    for i in range(n):
        for j in range(n):
            if i==j: continue
            rel=((E[i]-E[j])*x)%1
            if 0<rel<F(1,2): A[i][j]=1
            elif rel>F(1,2): A[i][j]=0
            else: A[i][j]=1 if E[i]<E[j] else 0; tf=False
    return A,tf
def scores(A): return [sum(r) for r in A]
def c3_count(A):
    n=len(A); c=0
    for a,b,cc in itertools.combinations(range(n),3):
        if (A[a][b] and A[b][cc] and A[cc][a]) or (A[a][cc] and A[cc][b] and A[b][a]): c+=1
    return c
def Hpaths(A):
    n=len(A); h=0
    for p in itertools.permutations(range(n)):
        if all(A[p[t]][p[t+1]] for t in range(n-1)): h+=1
    return h
def strongly_connected(A):
    n=len(A)
    # reachability from 0 forward and reverse
    def reach(adj,s):
        seen={s}; st=[s]
        while st:
            u=st.pop()
            for v in range(n):
                if adj(u,v) and v not in seen: seen.add(v); st.append(v)
        return seen
    fwd=reach(lambda u,v:A[u][v],0)
    bwd=reach(lambda u,v:A[v][u],0)
    return len(fwd)==n and len(bwd)==n
def maxgap(E,x):
    ps=sorted(set((e*x)%1 for e in E))
    if len(ps)==1: return F(1)
    g=F(0)
    for i in range(len(ps)):
        nxt=ps[(i+1)%len(ps)]+(F(1) if i+1==len(ps) else F(0)); g=max(g,nxt-ps[i])
    return g
def sectors_hit(E,x): return set(int(((e*x)%1)*7) for e in E)
def breakpoints(E):
    bps=set([F(0),F(1)]); Es=sorted(set(E))
    for e in Es:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
        for m in range(0,2*e+1): bps.add(F(m,2*e))
    diffs=set()
    for a in range(len(Es)):
        for b in range(a+1,len(Es)):
            diffs.add(Es[b]-Es[a]); diffs.add(Es[a]+Es[b])
    for d in diffs:
        if d==0: continue
        for m in range(0,2*d+1): bps.add(F(m,2*d))
    return sorted(x for x in bps if 0<=x<=1)

def full_profile(E):
    k=len(E); bps=breakpoints(E)
    M={"S7":F(0),"N":F(0),"SC":F(0),"hasCondorcet":F(0)}
    A_={"c3":F(0),"H":F(0),"transtriples":F(0)}
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        w=x1-x0; xm=(x0+x1)/2
        A,tf=round_tournament(E,xm); s=scores(A); mg=maxgap(E,xm)
        if len(sectors_hit(E,xm))==7: M["S7"]+=w
        if mg<=F(1,7): M["N"]+=w
        if strongly_connected(A): M["SC"]+=w
        if (k-1) in s: M["hasCondorcet"]+=w
        c3v=c3_count(A); A_["c3"]+=w*c3v; A_["transtriples"]+=w*(comb(k,3)-c3v)
        if k<=8: A_["H"]+=w*Hpaths(A)
    return M,A_

print("="*100)
print("OCF / ODD-CYCLE ANGLE + AP=TRANSITIVE EXTREMALITY  (exact)")
print("="*100)
print(f"{'shape':<16}{'meas(S7)':>10}{'measN':>10}{'measSC':>10}{'P[Cond]':>9}"
      f"{'E[c3]':>9}{'E[trans3]':>11}{'E[H]':>10}")
print("-"*100)
shapesets={
 7:[("consec",list(range(7))),("AP3",[0,3,6,9,12,15,18]),
    ("perforated",[0,1,2,3,4,5,7]),("dissoc",[0,1,3,7,15,31,63]),
    ("Sidon",[0,1,3,7,12,20,30])],
 8:[("consec",list(range(8))),("AP3",[0,3,6,9,12,15,18,21]),
    ("perforated",[0,1,2,3,4,5,6,8]),("dissoc",[0,1,3,7,15,31,63,127]),
    ("Sidon",[0,1,3,7,12,20,30,44]),("random",[0,5,13,27,41,58,79,97])],
}
results={}
for k in (7,8):
    print(f"--- k={k} ---")
    for name,E in shapesets[k]:
        M,A_=full_profile(E)
        results[(k,name)]=(M,A_)
        muN=F(1)-M["N"]
        print(f"{name:<16}{float(M['S7']):>10.5f}{float(M['N']):>10.5f}"
              f"{float(M['SC']):>10.5f}{float(M['hasCondorcet']):>9.5f}"
              f"{float(A_['c3']):>9.4f}{float(A_['transtriples']):>11.4f}"
              f"{float(A_['H']):>10.3f}")

print("\n" + "="*100)
print("AP=TRANSITIVE EXTREMALITY readout:")
for k in (7,8):
    rows=[(name,results[(k,name)][1]['c3'],results[(k,name)][0]['S7'],
           results[(k,name)][0]['N']) for name,_ in shapesets[k]]
    minc3=min(r[1] for r in rows); maxS7=max(r[2] for r in rows); maxN=max(r[3] for r in rows)
    consec_c3=[r[1] for r in rows if r[0]=='consec'][0]
    consec_S7=[r[2] for r in rows if r[0]=='consec'][0]
    consec_N =[r[3] for r in rows if r[0]=='consec'][0]
    print(f"  k={k}: consec is min-E[c3]? {consec_c3==minc3}  (consec E[c3]={float(consec_c3):.4f}, "
          f"min={float(minc3):.4f})")
    print(f"        consec is max-meas(S7)? {consec_S7==maxS7}  max-measN? {consec_N==maxN}")
print("\nINTERPRETATION: if consec minimises E[c3] (=maximises transitive triples) AND")
print("maximises meas(S7)/measN, then the AP cluster is the one whose winding tournament is")
print("the MOST TRANSITIVE on average AND fills the most sectors — a clean parity-side")
print("statement of AP-extremality: 'the dangerous LRC cluster = the most-transitive winder'.")
print("="*100)
print("DONE.")
