#!/usr/bin/env python3
"""
lrc14_cyclespace_bridge_ADVERSARIAL_kps.py
kind-pasteur-2026-06-19  (ADVERSARIAL re-verification, independent code)

Goal: independently CHECK (not assert) the four claims of the cycle-space/
even-graph bridge angle, and hunt for counterexamples.

I deliberately AVOID the original script's closed-form shortcuts:
 - FORM (1): actually compute the integer kernel rank of the relation map via
   exact Hermite/Smith reduction over Z, not 'k-1' by fiat.  Also build the
   matroid (column matroid of [e_i]) and CHECK it is the uniform matroid
   U_{1,k} (rank 1, every singleton independent) and the affine one rank 2.
 - FORM (2): compute the GF(2) cycle-space dim by an independent rank routine.
 - FORM (3): the LOAD-BEARING test.  Verify the THM-501 Mayer expansion
       D(q,S)/q = (1-b)^k + sum_{T nonempty} (1-b)^{k-|T|} (-1)^{|T|}
                   sum_{ relations on T mod q } prod chat(t_v)
   by EXACT brute force of the deficit D(q,S) and EXACT brute force of the RHS
   for small synthetic sets, confirming the (-1)^{|T|} support-parity sign is
   the actual sign that makes the identity hold (this is the OCF-parity claim).
 - FORM (4): doubling order mod 7, QR/NQR split, and the (non)factorization,
   recomputed independently, plus a search over many sets for ANY clean factor.
"""
import itertools, math
from fractions import Fraction as F

# ---------------- exact integer linear algebra ----------------
def hnf_rank_Z(rows):
    """Integer rank of a matrix (list of int rows) via fraction-free Gaussian elim over Q.
    Rank over Q = rank over Z for the image; kernel rank = ncols - rank."""
    M=[[F(x) for x in r] for r in rows]
    if not M: return 0
    nr=len(M); nc=len(M[0]); rank=0
    for c in range(nc):
        piv=None
        for r in range(rank,nr):
            if M[r][c]!=0: piv=r; break
        if piv is None: continue
        M[rank],M[piv]=M[piv],M[rank]
        for r in range(nr):
            if r!=rank and M[r][c]!=0:
                f=M[r][c]/M[rank][c]
                M[r]=[M[r][j]-f*M[rank][j] for j in range(nc)]
        rank+=1
        if rank==nr: break
    return rank

def integer_kernel_rank(E):
    """rank_Z ker( n |-> sum n_i e_i ).  Map Z^k -> Z, matrix is the 1xk row [e_i].
    kernel rank = k - rank(row)."""
    k=len(E)
    r=hnf_rank_Z([[int(e) for e in E]])
    return k-r, r

def integer_affine_kernel_rank(E):
    """ker of the 2xk map n|->(sum n_i, sum n_i e_i)."""
    k=len(E)
    r=hnf_rank_Z([[1]*k,[int(e) for e in E]])
    return k-r, r

def gf2_rank(rows):
    rows=[r[:] for r in rows]; nr=len(rows); nc=len(rows[0]) if rows else 0
    rank=0
    for c in range(nc):
        piv=None
        for r in range(rank,nr):
            if rows[r][c]: piv=r;break
        if piv is None: continue
        rows[rank],rows[piv]=rows[piv],rows[rank]
        for r in range(nr):
            if r!=rank and rows[r][c]:
                rows[r]=[a^b for a,b in zip(rows[r],rows[rank])]
        rank+=1
        if rank==nr: break
    return rank

def gf2_cycle_dim(E, affine=False):
    k=len(E); rows=[[int(e)%2 for e in E]]
    if affine: rows.append([1]*k)
    return k-gf2_rank(rows)

# ---------------- matroid checks ----------------
def is_uniform_U1k(E):
    """column matroid of the 1xk matrix [e_i]: a singleton {i} is independent iff e_i!=0.
    With all e_i != 0 it is U_{1,k}: rank 1, every single column independent, every pair
    dependent."""
    cols=[int(e) for e in E]
    all_single_indep = all(c!=0 for c in cols)
    # every pair dependent in a 1-row matrix (rank<=1) -> automatically true
    rank = hnf_rank_Z([cols])
    return all_single_indep, rank

# ---------------- exact deficit + Mayer RHS (FORM 3, load-bearing) ----------------
def deficit_exact(S, q):
    """D(q,S) = #{a in Z/q : for all v in S, v*a mod q NOT in B_q}, B_q=+-{0..floor(q/14)}.
    Here we use band half-width w=floor(q/14) to mirror THM-501. We use 14 since LRC(14)."""
    w=q//14
    band=set((x%q) for x in range(-w,w+1))
    cnt=0
    for a in range(q):
        ok=True
        for v in S:
            if (v*a)%q in band: ok=False;break
        if ok: cnt+=1
    return cnt

def chat(t,q,w):
    """chat(t) = (1/q) sum_{x=-w..w} e_q(-t x), real Dirichlet kernel. Exact via cos sum."""
    # (1/q) * sum_{x=-w}^{w} exp(-2pi i t x/q) = (1/q)*Dirichlet = sin(pi t (2w+1)/q)/(q sin(pi t/q))
    # compute exactly as a real number sum of cosines (t x term + symmetric) -> use floats only
    import cmath
    s=0.0
    for x in range(-w,w+1):
        s+=math.cos(2*math.pi*t*x/q)
    return s/q

def mayer_rhs(S,q):
    """RHS of THM-501 expansion at modulus q (FINITE form, all t_v in 1..q-1, congruence mod q).
    D(q,S)/q ?= (1-b)^k + sum_{T!=0} (1-b)^{k-|T|} (-1)^{|T|}
                  sum_{(t_v): t_v in 1..q-1, sum t_v v ==0 mod q} prod chat(t_v).
    Returns the float RHS. (Brute over all t-tuples on each subset T; small q only.)"""
    k=len(S); w=q//14; b=(2*w+1)/q
    total=(1-b)**k
    for r in range(1,k+1):
        for T in itertools.combinations(range(k),r):
            sign=(-1)**r
            pref=(1-b)**(k-r)
            # sum over t_v in 1..q-1 with sum t_v S[v] == 0 mod q
            acc=0.0
            for ts in itertools.product(range(1,q),repeat=r):
                if sum(ts[i]*S[T[i]] for i in range(r))%q==0:
                    p=1.0
                    for i in range(r): p*=chat(ts[i],q,w)
                    acc+=p
            total+=pref*sign*acc
    return total

# ---------------- FORM 4 doubling ----------------
def doubling_order_mod7():
    z=2; o=1
    while True:
        if (2**o)%7==1: return o
        o+=1
        if o>10: return None

def main():
    print("="*84)
    print("ADVERSARIAL re-verification of LRC(14) cycle-space / even-graph bridge")
    print("="*84)

    Es={
        "consec {0..8}":      list(range(0,9)),
        "consec {0..13}":     list(range(0,14)),
        "AP {0,2,4,6,8}":     [0,2,4,6,8],
        "Sidon {0,1,4,9,11}": [0,1,4,9,11],
        "two-large {0,1,500}":[0,1,500],
        "no-zero {3,5,7,11}": [3,5,7,11],
        "rand {2,3,7,13,29}": [2,3,7,13,29],
    }

    print("\n--- FORM (1): integer relation-lattice rank, computed (not asserted) ---")
    ok1=True
    for name,E in Es.items():
        k=len(E)
        kr,ir=integer_kernel_rank(E)
        ar,iar=integer_affine_kernel_rank(E)
        # expected: lin kernel = k-1 (img rank 1 since some e!=0); affine = k-2 if distinct & some !=0
        exp_lin = k-(1 if any(e!=0 for e in E) else 0)
        distinct = len(set(E))>1
        exp_aff = k-(2 if distinct else 1)
        single_indep, mrank = is_uniform_U1k(E)
        flag = "OK" if (kr==exp_lin and ar==exp_aff) else "MISMATCH"
        if flag=="MISMATCH": ok1=False
        print(f"  {name:22s} k={k}: ker_lin={kr}(exp {exp_lin}) ker_aff={ar}(exp {exp_aff}) "
              f"matroid_rank={mrank} singles_indep={single_indep}  [{flag}]")
    print(f"  FORM(1) all-match: {ok1}")

    print("\n--- FORM (2): GF(2) cycle-space (even-subgraph) dims ---")
    for name,E in Es.items():
        k=len(E)
        d=gf2_cycle_dim(E); da=gf2_cycle_dim(E,affine=True)
        print(f"  {name:22s} k={k}: gf2_lin={d}  gf2_aff={da}")

    print("\n--- FORM (3) LOAD-BEARING: brute deficit vs Mayer (-1)^|T| expansion ---")
    print("  Verifying THM-501's support-parity sign IS the correct sign (OCF-parity claim).")
    # small synthetic sets and small q so the full t-tuple brute is feasible
    tests=[([1,2,3],14),([1,2,3],21),([1,3,5],14),([1,2,4],28),([2,3,7],21)]
    ok3=True
    for S,q in tests:
        D=deficit_exact(S,q)
        lhs=D/q
        rhs=mayer_rhs(S,q)
        match = abs(lhs-rhs)<1e-9
        if not match: ok3=False
        print(f"  S={S} q={q}: D/q={lhs:.10f}  Mayer-RHS={rhs:.10f}  diff={abs(lhs-rhs):.2e}  "
              f"{'MATCH' if match else 'FAIL'}")
    # CONTROL: flip the sign to (+1)^|T| and show it BREAKS (proves the parity is load-bearing)
    def mayer_rhs_wrongsign(S,q):
        k=len(S); w=q//14; b=(2*w+1)/q; total=(1-b)**k
        for r in range(1,k+1):
            for T in itertools.combinations(range(k),r):
                pref=(1-b)**(k-r)
                acc=0.0
                for ts in itertools.product(range(1,q),repeat=r):
                    if sum(ts[i]*S[T[i]] for i in range(r))%q==0:
                        p=1.0
                        for i in range(r): p*=chat(ts[i],q,w)
                        acc+=p
                total+=pref*(+1)*acc   # WRONG: no (-1)^{|T|}
        return total
    S,q=[1,2,3],14
    print(f"  CONTROL (drop the (-1)^|T| sign): wrong-RHS={mayer_rhs_wrongsign(S,q):.10f} "
          f"vs D/q={deficit_exact(S,q)/q:.10f}  -> {'still matches?!' if abs(mayer_rhs_wrongsign(S,q)-deficit_exact(S,q)/q)<1e-9 else 'breaks (sign is load-bearing)'}")
    print(f"  FORM(3) Mayer-with-parity all-match: {ok3}")

    print("\n--- FORM (4): doubling order mod 7, QR/NQR, factorization hunt ---")
    print(f"  doubling order of 2 mod 7 = {doubling_order_mod7()}  (claim: 3)")
    QR={1,2,4}; NQR={3,5,6}
    print(f"  <2> orbit from 1: {[ (2**i)%7 for i in range(3)]} -> QR={sorted(QR)}")
    print(f"  <2> orbit from 3: {[ (3*2**i)%7 for i in range(3)]} -> NQR={sorted(NQR)}")
    # factorization hunt: does covAll == covQR*covNQR / cov0 ever hold?
    def meas_cover(E,target):
        E=sorted(set(int(e) for e in E));
        bps={F(0),F(1)}
        for e in E:
            if e==0: continue
            for a in range(0,7*e+1): bps.add(F(a,7*e))
        bps=sorted(b for b in bps if 0<=b<=1); tot=F(0)
        for i in range(len(bps)-1):
            lo,hi=bps[i],bps[i+1]
            if hi==lo: continue
            v=(lo+hi)/2
            hit=frozenset(((F(e)*v - (F(e)*v).numerator//(F(e)*v).denominator).numerator*7)//
                          ((F(e)*v - (F(e)*v).numerator//(F(e)*v).denominator).denominator) for e in E)
            if set(target)<=hit: tot+=hi-lo
        return tot
    worst_ratio=0
    for name,E in Es.items():
        if len(E)>8 or max(E)>50: continue
        pall=meas_cover(E,set(range(7)))
        pQR=meas_cover(E,{0}|QR); pNQR=meas_cover(E,{0}|NQR)
        prod=pQR*pNQR
        ratio=float(pall/prod) if prod!=0 else float('nan')
        print(f"  {name:22s}: covAll={float(pall):.4f} covQR={float(pQR):.4f} covNQR={float(pNQR):.4f} "
              f"QR*NQR={float(prod):.4f}  ratio covAll/(QR*NQR)={ratio:.3f}")

if __name__=="__main__":
    main()
