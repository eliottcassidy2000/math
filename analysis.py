from fractions import Fraction
from itertools import combinations
from math import gcd
from functools import reduce
from collections import Counter
from lrc import frac_dist, D, witness_times, M_witness, n_of

def is_primitive(sp):
    return reduce(gcd, sp) == 1

def pair_score(speeds, i, j):
    """Best balanced distance achievable at pair (i,j)'s witness times:
    max over k of D(k/(v_i+v_j)) restricted to times where the binding is on this pair.
    We define pair P's score = max over k of D(t) where t=k/(v_i+v_j),
    i.e. the best 'safe for everyone' distance at this pair's grid."""
    s = speeds[i] + speeds[j]
    best = Fraction(0); bestk=None
    for k in range(1, s):
        t = Fraction(k, s)
        d = D(t, speeds)   # min over ALL runners (everyone must be safe)
        if d > best:
            best = d; bestk=k
    return best, bestk

def pair_tournament(speeds):
    """Return dict pair->(score,k), and the winning pair."""
    n=len(speeds)
    scores={}
    for i,j in combinations(range(n),2):
        scores[(i,j)] = pair_score(speeds,i,j)
    win = max(scores, key=lambda p: scores[p][0])
    return scores, win

def analyze_set(speeds):
    m=len(speeds); n=n_of(speeds)
    scores,win = pair_tournament(speeds)
    M = scores[win][0]
    i,j = win
    psum = speeds[i]+speeds[j]
    return {
        'speeds':tuple(speeds),'n':n,'M':M,'Mn':M*n,
        'win':win,'win_sum':psum,'win_k':scores[win][1],
        'is_smallest_sum': psum==min(speeds[a]+speeds[b] for a,b in combinations(range(m),2)),
        'is_extremes': set(win)=={0,m-1},  # smallest & largest speed indices (speeds sorted)
        'sum_mod_n': psum % n,
        'tight': M*n==1,
        'sum_mult_n': psum % n == 0,
    }

# ---------- PART 1 & 3: winning pair statistics ----------
def part_1_3(maxspeed=12, ms=(4,5,6)):
    out={}
    for m in ms:
        rows=[]
        for combo in combinations(range(1,maxspeed+1), m):
            sp=list(combo)
            if not is_primitive(sp): continue
            rows.append(analyze_set(sp))
        out[m]=rows
    return out

if __name__=="__main__":
    data = part_1_3(12,(4,5,6))
    for m,rows in data.items():
        n=m+1
        print(f"\n===== m={m}, n={n}, primitive speed sets <= 12: {len(rows)} sets =====")
        # LRC check
        viol=[r for r in rows if r['Mn']<1]
        print("LRC M*n>=1 always:", len(viol)==0, "(violations:",len(viol),")")
        # tight <=> pair sum multiple of n
        tight=[r for r in rows if r['tight']]
        mult =[r for r in rows if r['sum_mult_n']]
        tight_iff = all(r['tight']==r['sum_mult_n'] for r in rows)
        print(f"tight (M*n==1): {len(tight)}; winning-sum multiple of n: {len(mult)}; tight<=>mult-of-n: {tight_iff}")
        # winning pair structure
        smallest=sum(r['is_smallest_sum'] for r in rows)
        extremes=sum(r['is_extremes'] for r in rows)
        print(f"winning pair is smallest-sum pair: {smallest}/{len(rows)} = {smallest/len(rows):.3f}")
        print(f"winning pair is extremes (min&max speed): {extremes}/{len(rows)} = {extremes/len(rows):.3f}")
        # distribution of winning sum relative to n
        print("winning_sum mod n distribution:", dict(sorted(Counter(r['sum_mod_n'] for r in rows).items())))
        # winning sum vs n
        ratio=Counter()
        for r in rows:
            ws=r['win_sum']
            ratio['<n' if ws<n else ('=n' if ws==n else '>n')]+=1
        print("winning_sum vs n:", dict(ratio))
        # M*n distribution (rounded)
        mns=Counter(round(float(r['Mn']),4) for r in rows)
        print("M*n distinct values (top):", dict(sorted(mns.items())[:8]))
