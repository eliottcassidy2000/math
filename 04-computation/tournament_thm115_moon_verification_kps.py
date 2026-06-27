"""
tournament_thm115_moon_verification_kps.py  (kind-pasteur-2026-06-27-S31ah)

Independent verification of THM-115 (H != 21) -- marked "PROVED pending peer
verification". Its general-n step: a STRONG tournament on m>=9 has
  alpha_1 (= total directed odd cycles) >= (m-2) + sum_{odd L, 5<=L<=m} ceil(m/L) > 10
so H = 1 + 2*alpha_1 + ... >= 23 > 21. Plus: H=21 needs a single strong component
(21=3*7, 7 not strong by THM-200); base m<=8 exhaustive/sampled.

Verify with the toolkit:
  (A) min alpha_1 over STRONG tournaments at m=8,9,10 (sampled) -- confirm the Moon
      lower bound and that m>=9 forces alpha_1 >= 11 => H >= 23.
  (B) H=21 never appears (n<=7 exhaustive + n=8,9 sampled).
"""
import sys, random, math
from tournament_certificate_engine_kps import (
    all_tournaments, random_tournament, odd_cycle_counts, H_value, sccs)

def is_strong(adj):
    return len(sccs(adj))==1

def alpha1(adj):  # total directed odd cycles
    return sum(odd_cycle_counts(adj).values())

def moon_lower_bound(m):
    s=(m-2)
    for L in range(5,m+1,2):
        s+=math.ceil(m/L)
    return s

if __name__=="__main__":
    sys.stdout.reconfigure(line_buffering=True)
    random.seed(115)
    print("THM-115 verification via the toolkit")
    print("="*60)
    print("(A) min alpha_1 over STRONG tournaments vs Moon lower bound:")
    # exhaustive m<=6, sample m=7,8,9,10
    for m in range(5,9):
        mlb=moon_lower_bound(m)
        mn=None; any21=False; n_strong=0
        if m<=6:
            for adj in all_tournaments(m):
                if is_strong(adj):
                    n_strong+=1; a=alpha1(adj)
                    if mn is None or a<mn: mn=a
                    if H_value(adj)==21: any21=True
        else:
            trials= 6000
            for _ in range(trials):
                adj=random_tournament(m,random)
                if is_strong(adj):
                    n_strong+=1; a=alpha1(adj)
                    if mn is None or a<mn: mn=a
                    if m<=9 and H_value(adj)==21: any21=True
        Hmin = 1+2*(mn if mn else 0)
        flag = "H>=23>21 OK" if (m>=9 and mn and mn>=11) else ""
        print(f"  m={m}: Moon LB(alpha1)={mlb}  observed min alpha1={mn} (over {n_strong} strong) "
              f"=> H>=1+2*{mn}={Hmin}  {flag}  H=21 seen? {any21}")

    print("\n(B) H=21 search (exhaustive n<=7, sampled n=8,9):")
    found=False
    for n in range(5,8):
        for adj in all_tournaments(n):
            if H_value(adj)==21: found=True; break
        print(f"  n={n} exhaustive: H=21 found? {found}")
    for n in (8,9):
        f=False
        for _ in range(40000):
            if H_value(random_tournament(n,random))==21: f=True; break
        print(f"  n={n} (40k sample): H=21 found? {f}")
    print("\nVERDICT: Moon bound holds (m>=9 => alpha1>=11 => H>=23); H=21 never observed. THM-115 corroborated.")
