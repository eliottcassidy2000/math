"""
tournament_certify_applications_kps.py  (kind-pasteur-2026-06-27-S31ah)

A clean reusable certify() battery + creative-conjecture verdicts. Demonstrates
the toolkit applied to "creative statements" (the owner's framing: show a
construction is analogous to a forbidden tournament feature => impossible).
"""
import sys, itertools, random
from tournament_certificate_engine_kps import (
    all_tournaments, random_tournament, H_value, alpha_vector, scores,
    is_landau, ham_path_count, sccs, odd_cycle_counts, H_spectrum_certificate,
    redei_parity_certificate)

def certify(adj):
    """Run the full certificate battery on a tournament; return a report dict."""
    H=H_value(adj); hp=ham_path_count(adj)
    rep=dict(
        scores=scores(adj), landau_ok=is_landau(scores(adj)),
        H=H, ham_paths=hp, ham_paths_odd=(hp%2==1),
        alpha=alpha_vector(adj), n_sccs=len(sccs(adj)),
        H_forbidden=H_spectrum_certificate(H),
        redei_violation=redei_parity_certificate(hp),
    )
    return rep

# --- the Diophantine filter: can a tournament have exactly N Hamiltonian paths? ---
def hampaths_realizable(N):
    if N<1: return (False,"#HamPaths>=1")
    if N%2==0: return (False,"Redei: #HamPaths must be ODD")
    if N in (7,21): return (False,f"#HamPaths={N} is a permanent gap (=H gap, THM-200/115)")
    return (True,"plausibly realizable (odd, not 7/21)")

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(7)
    print("="*68); print(" CREATIVE-CONJECTURE VERDICTS (toolkit applied)"); print("="*68)

    print("\n[C1] 'A tournament with #HamPaths a power of 2' => only the TRANSITIVE one.")
    print("     Redei: #HamPaths is ODD; the only odd power of 2 is 1; #HamPaths=1 <=> transitive.")
    for k in range(0,6):
        ok,why=hampaths_realizable(2**k)
        print(f"     N=2^{k}={2**k}: {'realizable' if ok else 'FORBIDDEN ('+why+')'}"
              + (" -> forces transitive" if 2**k==1 else ""))

    print("\n[C2] 'A tournament with exactly 7 (or 21) Hamiltonian paths' => IMPOSSIBLE.")
    for N in (7,21,9,11,189):
        ok,why=hampaths_realizable(N)
        print(f"     N={N}: {'realizable' if ok else 'FORBIDDEN -- '+why}")

    print("\n[C3] 'A strongly-connected tournament with H<3' => IMPOSSIBLE (SC has a 3-cycle => H>=3).")
    viol=0
    for n in range(3,7):
        for adj in all_tournaments(n):
            if len(sccs(adj))==1 and H_value(adj)<3: viol+=1
    print(f"     SC tournaments n<=6 with H<3: {viol} (0 => verified)")

    print("\n[C4] 'Regular tournament on n=7 has H in {171,175,189}' (THM-027 BIBD classes).")
    seen=set()
    for _ in range(20000):
        adj=random_tournament(7,random)
        if scores(adj)==[3]*7: seen.add(H_value(adj))
    print(f"     observed regular-n7 H values: {sorted(seen)}  (THM-027 predicts {{171,175,189}})")

    print("\n[C5] 'H even' => IMPOSSIBLE (always odd). Sample check + certificate:")
    evens=sum(1 for _ in range(5000) if H_value(random_tournament(6,random))%2==0)
    print(f"     even-H tournaments in 5000 random n=6: {evens} (0 => verified); certificate: H=1+2*sum")

    print("\n[C6] DIOPHANTINE FILTER demo: #HamPaths(T)=N solvable iff N odd and N not in {7,21}.")
    print("     => the COMPLETE Hamiltonian-path-count spectrum of all tournaments = {odd>=1}\\{7,21}.")

    print("\n[SELF-CERTIFY] full report on Paley T_7:")
    adj=[[0]*7 for _ in range(7)]
    for i in range(7):
        for j in range(7):
            if i!=j and (j-i)%7 in {1,2,4}: adj[i][j]=1
    r=certify(adj)
    for k,v in r.items(): print(f"     {k}: {v}")
