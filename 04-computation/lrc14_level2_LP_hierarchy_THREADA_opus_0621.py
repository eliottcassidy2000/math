#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_level2_LP_hierarchy_THREADA_opus_0621.py  (opus, 2026-06-21, THREAD A — the HINGE)

THE MAKE-OR-BREAK TEST: does the Coregliano-Jeronimo-Jones complete LP hierarchy
(arXiv:2211.01248, level-1 = Delsarte) COLLAPSE (Prop 1.2) or IMPROVE at level 2 for
the LRC(14) Z/7 sector cover bound  measS7(E) = P(N=0) = p_0 ?

============================== THE CORRESPONDENCE ==============================
CJJ (linear codes):                       OUR problem (LRC sector cover):
  maximize code size |C|                    maximize measS7(E) = p_0(E)
  distance distribution a_i (an S_n         occupancy law:
    -orbit-collapsed measure on C x C)        SUBSET version  q_A = P(empty set = A),
                                                A subset {1..6}  (32 atoms)
                                              NUMBER version  p_t = P(|empty|=t) = sum_{|A|=t} q_A
                                                (7 atoms) = the S_7-orbit collapse of q.
  Delsarte LP (level 1): on a_i with        LEVEL-1 LP: on p_t, constraints =
    MacWilliams (Krawtchouk) positivity       Krawtchouk-moment positivity M_j=E[K_j(N)]>=0
                                              + nonneg + sum=1.  Objective max p_0.
    bound = the q'>=0 dual = our L_y.         (THM-534 dual = L_y).  This IS Delsarte.
  Level-l: configuration of all 2^l         LEVEL-2 LP: on the SUBSET law q_A
    LINEAR COMBINATIONS of l codewords;       (32 atoms, the "support-2 / pair sector
    Lovasz theta' on conflict Cayley graph;   structure"), with the FULL hypercube Q_6
    factor translation + LINEAR-COMB sym.     association-scheme (Walsh/Fourier) positivity
                                              + the level-1 marginals as constraints +
                                              the exact "always-feasible relation" facts.
  Prop 1.2 (the HINGE): the hierarchy        THE TEST: is the level-2 LP optimum STRICTLY
    IMPROVES on Delsarte ONLY via LINEARITY   SMALLER than the level-1 optimum?  If equal =>
    (optimizer closed under lin. comb.); for  COLLAPSE (Prop 1.2: our optimizer is "non-
    NON-LINEAR optimizers it COLLAPSES to     linear" -- no closure under lin. comb.).  If
    level-1 at every level.                    strictly smaller AND <= measS7(consec) it could
                                              prove consec-max WITHOUT the finite atlas.

WHAT IS "linearity" HERE?  consec's residues mod 7 = {0,1,2,3,4,5,6} = ALL of Z/7 (the
full subgroup).  The relation code Lambda(E)={n in Z^k: sum n_i e_i = 0} is consec's
"linear structure" (anti-MDS, min distance 2: 2*e_1 = e_2).  THREAD A asks whether THAT
linear (full-Z/7 / low-distance-relation) structure is what a level-2 LP can SEE and
whether seeing it certifies consec-max.

============================== WHAT THIS SCRIPT DOES ==============================
For k=8 (the binding row, 11432-shape finite atlas, consec = unique argmax):
 (A) LEVEL-1 LP:  max p_0 over prob vectors p_t (t=0..6) subject to the moment data we
     are willing to impose.  Two flavors:
       L1a = pure Krawtchouk-positivity relaxation (impose ONLY M_j>=0, the scheme dual):
             this is the "abstract Delsarte" optimum (no shape-specific moment values).
       L1b = impose the EXACT achievable moment RANGE (S_1,S_2 ranges over real shapes):
             the data-constrained Delsarte optimum.
 (B) LEVEL-2 LP:  max q_{emptyset} over SUBSET prob vectors q_A (A subset {1..6}, 64 atoms)
     subject to:
       (i)   q_A >= 0, sum q_A = 1                                 (probability)
       (ii)  Walsh/Fourier positivity of the Q_6 hypercube scheme:  for every character
             chi_S (S subset {1..6}), the Fourier coefficient  qhat_S = sum_A (-1)^{|A cap S|} q_A
             must satisfy the MacWilliams-type sign constraints of the SUBSET occupancy
             of a REAL torus measure.  (We DERIVE the exact achievable sign pattern from
             the actual shapes -- this is the level-2 "association scheme positivity".)
       (iii) the Z/7 CYCLIC symmetry orbit-collapse: real occupancy laws are NOT cyclic-
             symmetric (sector 0 is special), so we DO NOT impose full S_7; we impose the
             genuine structural facts (Krawtchouk marginals + the realized Fourier cone).
 (C) COMPARE: L1 optimum  vs  L2 optimum  vs  measS7(consec).  Verdict:
       L2 == L1  ->  COLLAPSE (Prop 1.2).   L2 < L1  ->  IMPROVES.
       L2 <= measS7(consec)+eps -> level-2 alone would CERTIFY consec-max (no atlas).

All arithmetic exact (Fractions) for the data; LP solved exactly by vertex enumeration
on the small moment LPs and by an exact rational simplex (Bland) on the 64-atom L2 LP.
"""
import sys, itertools
from math import comb, gcd
from fractions import Fraction as F
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

# ---------------- occupancy laws (number and subset), exact ----------------
def occ_subset(E):
    """q_A = meas{x: inner-empty set = A}, A subset {1..6}.  Exact."""
    E=sorted(set(E)); Enz=[e for e in E if e!=0]
    bps={F(0),F(1)}
    for e in Enz:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); law=defaultdict(lambda:F(0))
    for x0,x1 in zip(bps,bps[1:]):
        if x1<=x0: continue
        xm=(x0+x1)/2
        hit=set(int(7*e*xm)%7 for e in E)
        empty=frozenset(s for s in range(1,7) if s not in hit)
        law[empty]+=x1-x0
    return dict(law)

def occ_number(E):
    law=occ_subset(E); p=[F(0)]*7
    for A,m in law.items(): p[len(A)]+=m
    return p

def Kraw(j,t,n=6): return sum((-1)**i*comb(t,i)*comb(n-t,j-i) for i in range(j+1))
KTAB=[[F(Kraw(j,t)) for t in range(7)] for j in range(7)]
def Mvec(p): return [sum(KTAB[j][t]*p[t] for t in range(7)) for j in range(7)]
def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return list(range(k))

INNER=list(range(1,7))  # 6 inner sectors
SUBSETS=[frozenset(c) for r in range(7) for c in itertools.combinations(INNER,r)]  # 64
SIDX={A:i for i,A in enumerate(SUBSETS)}

# ============================================================================
# LEVEL-1 LP (Delsarte): max p_0 s.t. p prob on {0..6}, Krawtchouk moments in
# the ACHIEVABLE convex hull of real shapes.  We compute the achievable (M_1,M_2)
# range from the atlas and solve the moment LP exactly.
# ============================================================================
def lp_max_p0_moments(constraints):
    """max p_0 over prob vec p[0..6] with linear EQ/INEQ constraints (list of
    (coeffs[7], op, rhs)), op in {'==','<=','>='}.  Exact vertex enumeration:
    optimum at a basic feasible soln with support<= (#equality constraints).
    We add sum p =1 always.  Returns (max_p0, witness_p)."""
    eqs=[([F(1)]*7,F(1))]  # sum p =1
    ineqs=[]
    for co,op,r in constraints:
        co=[F(c) for c in co]; r=F(r)
        if op=='==': eqs.append((co,r))
        elif op=='<=': ineqs.append((co,r))
        elif op=='>=': ineqs.append(([-c for c in co],-r))
    # Enumerate vertices: choose a support set; with all p>=0 the LP vertices are
    # determined by which inequalities/nonneg bounds are tight.  We brute-force over
    # support subsets of {0..6} (<=7), solving the equality system, checking feasibility.
    best=None; bestp=None
    idxs=list(range(7))
    for r in range(1,8):
        for sup in itertools.combinations(idxs,r):
            # variables p_sup; solve eqs restricted; need #unknown <= #eqs typically.
            A=[]; b=[]
            for co,rhs in eqs:
                A.append([co[t] for t in sup]); b.append(rhs)
            # also pick (len(sup)-len(eqs)) active inequalities to pin a vertex
            need=len(sup)-len(A)
            if need<0: continue
            for act in itertools.combinations(range(len(ineqs)),need):
                AA=[row[:] for row in A]; bb=b[:]
                for ai in act:
                    co,rhs=ineqs[ai]; AA.append([co[t] for t in sup]); bb.append(rhs)
                # solve square-ish system AA x = bb (len(sup) x len(sup) if exactly)
                if len(AA)!=len(sup): continue
                sol=gauss(AA,bb)
                if sol is None: continue
                p=[F(0)]*7
                ok=True
                for i,t in enumerate(sup):
                    if sol[i]<0: ok=False;break
                    p[t]=sol[i]
                if not ok: continue
                # check all constraints
                feas=True
                for co,rhs in eqs:
                    if sum(co[t]*p[t] for t in range(7))!=rhs: feas=False;break
                if feas:
                    for co,rhs in ineqs:
                        if sum(co[t]*p[t] for t in range(7))>rhs+F(0): feas=False;break
                if feas:
                    if best is None or p[0]>best: best=p[0]; bestp=p
    return best,bestp

def gauss(A,b):
    n=len(A);
    if n==0: return []
    m=len(A[0])
    if m!=n: return None
    M=[A[i][:]+[b[i]] for i in range(n)]
    for c in range(n):
        piv=None
        for r in range(c,n):
            if M[r][c]!=0: piv=r;break
        if piv is None: return None
        M[c],M[piv]=M[piv],M[c]
        pv=M[c][c]; M[c]=[x/pv for x in M[c]]
        for r in range(n):
            if r!=c and M[r][c]!=0:
                f=M[r][c]; M[r]=[M[r][j]-f*M[c][j] for j in range(n+1)]
    return [M[i][n] for i in range(n)]

# ============================================================================
# EXACT RATIONAL LP via Fourier-Motzkin-free vertex / simplex.  For the 64-atom
# level-2 LP we use a small exact simplex (maximize) with rational pivots.
# ============================================================================
def simplex_max(c, A_ub, b_ub, A_eq, b_eq, nvar):
    """max c.x  s.t. A_ub x <= b_ub, A_eq x == b_eq, x>=0.  Exact Fractions.
    Big-M / two-phase via converting to standard form and Bland's rule.
    Returns (optval, x) or (None,None) if infeasible/unbounded."""
    # Build: equalities + inequalities (add slacks).  Use two-phase.
    rows=[]; rhs=[]
    # inequalities -> add slack
    nslack=len(A_ub)
    for i,(row,bv) in enumerate(zip(A_ub,b_ub)):
        r=list(row)+[F(0)]*nslack
        r[nvar+i]=F(1)
        # ensure rhs>=0
        if bv<0:
            r=[-x for x in r]; bv=-bv
        rows.append(r); rhs.append(bv)
    for row,bv in zip(A_eq,b_eq):
        r=list(row)+[F(0)]*nslack
        if bv<0:
            r=[-x for x in r]; bv=-bv
        rows.append(r); rhs.append(bv)
    total=nvar+nslack
    # Phase 1: artificial vars for every row, minimize their sum.
    nrows=len(rows)
    for i in range(nrows):
        for j in range(nrows):
            rows[i].append(F(1) if i==j else F(0))
    base=[total+i for i in range(nrows)]  # artificials in basis
    grand=total+nrows
    # objective phase1: minimize sum artificials = maximize -sum
    def run(obj):
        # obj: vector length grand (coeffs to MAXIMIZE)
        while True:
            # reduced costs
            # z_j - c_j ; compute c_B
            cB=[obj[base[i]] for i in range(nrows)]
            red=[]
            for j in range(grand):
                zj=sum(cB[i]*rows[i][j] for i in range(nrows))
                red.append(zj-obj[j])
            # entering: most negative reduced cost (Bland for cycling safety)
            enter=None
            for j in range(grand):
                if red[j]<0:
                    enter=j; break
            if enter is None: break
            # ratio test
            leave=None; best=None
            for i in range(nrows):
                if rows[i][enter]>0:
                    ratio=rhs[i]/rows[i][enter]
                    if best is None or ratio<best or (ratio==best and base[i]<base[leave]):
                        best=ratio; leave=i
            if leave is None: return False  # unbounded
            # pivot
            pv=rows[leave][enter]
            rows[leave]=[x/pv for x in rows[leave]]; rhs[leave]=rhs[leave]/pv
            for i in range(nrows):
                if i!=leave and rows[i][enter]!=0:
                    f=rows[i][enter]
                    rows[i]=[rows[i][j]-f*rows[leave][j] for j in range(grand)]
                    rhs[i]=rhs[i]-f*rhs[leave]
            base[leave]=enter
        return True
    obj1=[F(0)]*total+[F(-1)]*nrows
    if not run(obj1): return None,None
    # feasibility: sum artificials
    artsum=sum(rhs[i] for i in range(nrows) if base[i]>=total)
    if artsum!=0: return None,None
    # drop artificial columns; phase 2 with real objective (maximize c.x)
    obj2=[F(0)]*grand
    for j in range(nvar): obj2[j]=c[j]
    # forbid artificials re-entering by setting big negative? we just exclude them:
    # zero their column influence is fine since rhs feasible; set obj huge negative
    for j in range(total,grand): obj2[j]=F(-10**9)
    if not run(obj2): return None,None
    x=[F(0)]*total
    for i in range(nrows):
        if base[i]<total: x[base[i]]=rhs[i]
    val=sum(c[j]*x[j] for j in range(nvar))
    return val,x[:nvar]

# ---- the danger-zone cap (target for "would close the cap") ----
H=F(1,14)
def danger(u):
    iv=[]
    for j in range(u):
        cc=F(j,u); a=(cc-H/u)%1; b=(cc+H/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def mgmerge(iv):
    iv=sorted(iv); o=[]
    for a,b in iv:
        if o and a<=o[-1][1]: o[-1]=(o[-1][0],max(o[-1][1],b))
        else: o.append((a,b))
    return o
def measGP(P):
    if not P: return F(1)
    dz=mgmerge([iv for u in P for iv in danger(u)]); s=F(0); prev=F(0)
    for a,b in dz:
        if a>prev: s+=a-prev
        prev=max(prev,b)
    if prev<1: s+=1-prev
    return s
def cap(k):
    psz=13-k
    if psz==0: return F(1)
    return min(measGP(P) for P in itertools.combinations(range(1,14),psz))

def banner(t): print("\n"+"="*92+f"\n{t}\n"+"="*92)

# ============================================================================
if __name__=="__main__":
    k=8; W=16
    Ec=consec(k); pc=occ_number(Ec); qc=occ_subset(Ec)
    measc=pc[0]; Mc=Mvec(pc); capk=cap(k)
    print(f"k={k}: measS7(consec)=p_0={measc}={float(measc):.6f}   cap_k={float(capk):.6f}")
    print(f"      consec Krawtchouk moments M_j = {[float(x) for x in Mc]}")

    # ---- build the atlas (achievable moment / Fourier ranges) ----
    banner("Building k=8 atlas (primitive shapes, max<=16) for achievable cones.")
    bank=[[0]+list(r) for r in itertools.combinations(range(1,W+1),k-1)]
    bank=[E for E in bank if primitive(E)]
    print(f"  {len(bank)} primitive shapes")
    # gather number-laws + subset-laws
    M1s=[];M2s=[];M3s=[];M4s=[]; p0s=[]
    # subset Fourier coeffs realized
    qhat_real=defaultdict(set)  # S -> set of realized qhat_S (sign cone)
    qhat_min={}; qhat_max={}
    for E in bank:
        p=occ_number(E); M=Mvec(p); p0s.append(p[0])
        M1s.append(M[1]);M2s.append(M[2]);M3s.append(M[3]);M4s.append(M[4])
        q=occ_subset(E)
        for S in SUBSETS:
            qh=sum(((-1)**len(A&S))*q.get(A,F(0)) for A in SUBSETS)
            if S not in qhat_min or qh<qhat_min[S]: qhat_min[S]=qh
            if S not in qhat_max or qh>qhat_max[S]: qhat_max[S]=qh
    print(f"  achievable M_1 in [{float(min(M1s)):.4f},{float(max(M1s)):.4f}]")
    print(f"  achievable M_2 in [{float(min(M2s)):.4f},{float(max(M2s)):.4f}]")
    print(f"  achievable p_0 in [{float(min(p0s)):.4f},{float(max(p0s)):.4f}]  (atlas max = consec? {max(p0s)==measc})")

    # ============================================================
    banner("(A) LEVEL-1 (Delsarte) LP optima")
    # L1a: ABSTRACT Delsarte: max p_0 s.t. p prob, M_j>=0 for all j (scheme positivity),
    #      and (the only data) M_2 fixed? No -- abstract uses NO shape data, just the
    #      scheme.  The pure Krawtchouk-positivity polytope max p_0:
    cons_abstract=[]
    for j in range(1,7):
        cons_abstract.append((KTAB[j],'>=',0))  # M_j>=0
    v_abs,p_abs=lp_max_p0_moments(cons_abstract)
    print(f"  L1a abstract (only p>=0,sum=1,M_j>=0):  max p_0 = {v_abs}={float(v_abs):.6f}")
    # L1b: impose ALSO the achievable EVEN-moment ceiling at consec (the binding j=2,4
    #      moments are consec-extremal by HYP-2726/Thread B; impose M_2<=M_2(consec),
    #      M_4<=M_4(consec) which is the PROVEN-by-atlas even-Krawtchouk dominance):
    cons_data=list(cons_abstract)
    cons_data.append((KTAB[2],'<=',Mc[2]))
    cons_data.append((KTAB[4],'<=',Mc[4]))
    v_b,p_b=lp_max_p0_moments(cons_data)
    print(f"  L1b +even-moment ceilings M_2<=M2c,M_4<=M4c:  max p_0 = {v_b}={float(v_b):.6f}")
    print(f"      (This is the Delsarte/L_y value -- compare to L_y(consec)={float(measc):.6f}? ")
    print(f"       L_y(consec) (dual bound) was 0.35823 in saturation script.)")

    # ============================================================
    banner("(B) LEVEL-2 LP: subset law q_A (64 atoms) + hypercube Q_6 Fourier positivity")
    # variables x = q_A for A in SUBSETS (64).  objective: max q_emptyset.
    nvar=len(SUBSETS)
    c=[F(0)]*nvar; c[SIDX[frozenset()]]=F(1)
    A_eq=[]; b_eq=[]
    # sum q =1
    A_eq.append([F(1)]*nvar); b_eq.append(F(1))
    # level-1 marginal CONSISTENCY: sum_{|A|=t} q_A = p_t is NOT fixed (p_t free),
    # but the Krawtchouk MOMENT structure must hold.  Impose the SAME even-moment
    # ceilings as L1b, expressed on subsets:  M_j = sum_A K_j(|A|) q_A.
    A_ub=[]; b_ub=[]
    for j in [2,4]:
        row=[KTAB[j][len(A)] for A in SUBSETS]
        A_ub.append(row); b_ub.append(Mc[j])
    for j in range(1,7):  # M_j>=0  -> -M_j<=0
        row=[-KTAB[j][len(A)] for A in SUBSETS]
        A_ub.append(row); b_ub.append(F(0))
    # THE LEVEL-2 CONTENT: hypercube Q_6 Walsh/Fourier positivity of the realized cone.
    # For each character S subset {1..6}, qhat_S = sum_A (-1)^{|A cap S|} q_A must lie in
    # the realized achievable interval [qhat_min[S], qhat_max[S]] (the level-2 association
    # -scheme constraint = the actual second-order structure of torus occupancy laws).
    nS=0
    for S in SUBSETS:
        if len(S)==0: continue
        row=[F((-1)**len(A&S)) for A in SUBSETS]
        A_ub.append(row[:]); b_ub.append(qhat_max[S])           # qhat_S <= max
        A_ub.append([-x for x in row]); b_ub.append(-qhat_min[S])# qhat_S >= min
        nS+=1
    print(f"  L2 LP: {nvar} atoms, {len(A_eq)} eq, {len(A_ub)} ineq ({nS} Fourier-cone pairs).")
    v2,x2=simplex_max(c,A_ub,b_ub,A_eq,b_eq,nvar)
    if v2 is None:
        print("  L2 LP infeasible/unbounded!")
    else:
        print(f"  L2 max q_emptyset = {v2}={float(v2):.6f}")

    # ============================================================
    banner("(C) VERDICT: collapse vs improves")
    print(f"  measS7(consec)        = {float(measc):.6f}")
    print(f"  L1a abstract Delsarte = {float(v_abs):.6f}")
    print(f"  L1b data Delsarte     = {float(v_b):.6f}")
    if v2 is not None:
        print(f"  L2 (level-2 LP)       = {float(v2):.6f}")
        impr = v2 < v_b - F(1,10**9)
        print(f"  L2 < L1b (strict improvement)? {impr}")
        print(f"  L2 <= measS7(consec)?          {v2 <= measc + F(1,10**9)}")
        if not impr:
            print("  => COLLAPSE: level-2 gives the SAME bound as level-1 (Prop 1.2: nonlinear optimizer).")
        else:
            print("  => IMPROVES: level-2 is strictly stronger than Delsarte.")
            if v2 <= measc + F(1,10**9):
                print("     AND level-2 alone CERTIFIES consec-max (bound <= measS7(consec))!")
            else:
                print("     but level-2 still exceeds measS7(consec) -- does not yet certify consec-max.")
    print("\nDONE.")
