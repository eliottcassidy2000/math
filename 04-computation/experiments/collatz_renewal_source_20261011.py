#!/usr/bin/env python3
"""Exact controls for renewal/source separation; no universal stopping claim.

Universe: coefficient first-passage words with total valuation <= 11;
literal cylinder enumeration through total 8; pair concatenations through 9;
exact first-passage mass through length 120; explicit clock/anchor hostiles.
Atomic controls: r=1..4, indices0..31; alias heads at m=0..15 and
q=32,64,128, plus adaptive schedule indices1..100.
All gates stay active under python -O. Standard library only.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import json

checks = 0


def check(ok, message):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(message)


def step(n):
    t = 3*n + 1
    a = (t & -t).bit_length() - 1
    return t >> a, a


def replay(n, word):
    values = [n]
    for a in word:
        n, actual = step(n)
        if actual != a:
            return None
        values.append(n)
    return values


def frame(word):
    p, b, q = 1, 0, 1
    for a in word:
        p, b, q = 3*p, 3*b + q, q << a
    modulus = 2*q
    residue = ((q-b)*pow(p, -1, modulus)) % modulus
    return p, b, q, residue


def code(max_total):
    out = []
    def extend(w, total):
        for a in range(1, max_total-total+1):
            v, t = w+(a,), total+a
            if (1 << t) > 3**len(v):
                out.append(v)
            else:
                extend(v, t)
    extend((), 0)
    return out


def law(depth):
    alive = {0: 1}
    mass = F(0)
    rows = []
    for length in range(1, depth+1):
        amin = (3**length).bit_length()
        count = sum(alive.values())
        exit_mass = F(count, 1 << (amin-1))
        next_alive = defaultdict(int)
        for old, multiplicity in alive.items():
            for a in range(1, amin-old):
                next_alive[old+a] += multiplicity
        survival = sum((F(v, 1 << k) for k, v in next_alive.items()), F(0))
        mass += exit_mass
        check(mass+survival == 1, "exact Kraft prefix budget")
        rows.append((length, amin, count, exit_mass, survival))
        alive = dict(next_alive)
    return rows


def first_descent(n, limit=10000):
    start, word = n, []
    for _ in range(limit):
        n, a = step(n)
        word.append(a)
        if n < start:
            return n, tuple(word)
    raise RuntimeError("finite control exceeded declared horizon")


def main():
    words = code(11)
    for w in words:
        p,b,q,r = frame(w)
        check(p < q and r % 2 == 1, "contracting odd address")
        for shift in (1, 2, 9):
            n = r + 2*q*shift
            path = replay(n, w)
            check(path is not None, "native exact word")
            check(q*path[-1] == p*n+b, "independent affine replay")
            check(all(v > n for v in path[1:-1]), "proper prefixes rise")
            check(path[-1] < n, "cofinal first descent")
        if sum(w) <= 8:
            literal = [n for n in range(1, 2*q, 2) if replay(n,w) is not None]
            check(literal == [r], "exhaustive exact cylinder")
    pairs = 0
    for u,v in product(words, repeat=2):
        if sum(u)+sum(v) > 9:
            continue
        p,b,q,r = frame(u+v)
        n = r+1000*2*q
        middle, observed_u = first_descent(n)
        end, observed_v = first_descent(middle)
        check((observed_u, observed_v) == (u,v), "renewal word order")
        check(q*end == p*n+b, "renewal endpoint")
        pairs += 1
    rows = law(120)
    check([r[3] for r in rows[:5]] == [F(1,2),F(1,8),F(1,8),F(3,64),F(7,128)],
          "independent law controls")
    check(2049**2>2*1444**2,"rational majorant of Chernoff constant")
    schedule=[]
    for i in range(1,101):
        target=F(1,16*(i+2)**2)
        t,bound=0,F(1)
        while bound>target:
            t,bound=t+1,bound*F(19,20)
        check(bound<=target<bound*F(20,19),"minimal rational adaptive cutoff")
        if i<=10:
            schedule.append(t)
    check(sum(w == (1,1,3) for w in words) == 1, "hard word retained")

    # A shortcut coefficient exit can end even; it is not an odd block.
    n = 13
    shortcut = [(3*n+1)//2]
    shortcut.append(shortcut[-1]//2)
    check(shortcut == [20,10] and step(13) == (5,3), "clock hostile 13")
    # Conditioning on fitting an available bit budget changes the law.
    fit_mass = sum((F(1, 1 << a) for a in (2,3)), F(0))
    check(F(1,4)/fit_mass == F(2,3), "conditional fit is renormalized")

    # All finite repetitions have positive native representatives, but their
    # compatible 2-adic limit is 19/5, not an ordinary integer.
    repeated = []
    for k in range(1, 13):
        w = (1,1,3)*k
        p,b,q,r = frame(w)
        check(F(b,q-p) == F(19,5), "common rational anchor")
        check((5*r-19) % (2*q) == 0, "2-adic anchor address")
        n = r+2*q
        path = replay(n,w)
        check(path is not None, "finite repetition has native source")
        check(q*(5*path[-1]-19) == p*(5*n-19), "integer separation identity")
        check(5*path[-1]-19 >= 6, "positive integer gap at anchor")
        repeated.append({"k":k,"residue":r,"modulus":2*q})
    u,v = (1,1,3),(2,)
    uv,vu = frame(u+v),frame(v+u)
    check(uv[0] == vu[0] and uv[2] == vu[2], "commutative slope budget")
    check(uv[1] != vu[1] and uv[3] != vu[3], "ordered carry and address")
    # Word (2) contracts affinely but has the root exception.
    check(replay(1,(2,)) == [1,1], "strict descent needs its threshold")
    for qcap in (1,3,8,32,100):
        length=qcap+1
        sums=[(3**j).bit_length()-1 for j in range(1,length)]
        sums += [(3**length).bit_length()]
        w=tuple(b-a for a,b in zip([0]+sums[:-1],sums))
        check(all(a>0 for a in w), "hard witness positive letters")
        check(all(2**sums[j-1] < 3**j for j in range(1,length)), "hard proper prefixes")
        exponent=qcap*sums[-1]-length
        check(exponent<0 or 2**exponent < 3**(qcap*length), "hard word exceeds cap")
    # A source-address prior, deliberately distinct from ROOT-supported mass.
    def nu(m,r):
        return F(1,(m+1)**r)-F(1,(m+2)**r)
    atomic=[]
    for r in range(1,5):
        total=F(0)
        for m in range(32):
            total+=nu(m,r)
            check(total==1-F(1,(m+2)**r),"rational prior telescopes")
        for m in range(16):
            for q in (32,64,128):
                head=sum((nu(m+k*q,r) for k in range(1,17)),F(0))
                bound=F(r+1,q**(r+1))
                check(head<bound,"finite alias head below proved total bound")
                atom=nu(m,r)
                floor=atom/(atom+bound)
                check(atom/(atom+head)>floor,"source survives prefix refinement")
                if r==2:
                    for target_score in (F(0),F(1,12),F(1,3)):
                        measurement=atom*target_score+head/3
                        lower=(measurement-bound/3)/atom
                        check(lower<=target_score<=lower+bound/(3*atom),
                              "score reconstruction with retained alias bill")
                        if target_score==0:
                            check(lower<0,"certified neighbours cannot create the missing target")
                if r==2 and m in (0,5) and q==128:
                    atomic.append({"index":m,"atom":str(atom),"q":q,
                                   "alias_bound":str(bound),"conditional_floor":str(floor)})
    check(nu(1,2)>F(3,16**3),"q>m is required: source17 has lower alias1 modulo16")
    check(F(3,64**3)==F(3,32**3)/8,"one more bit pays three powers in the alias bound")
    result={"status":"FINITE-EXACT controls; universal claims require the companion proof",
            "checks":checks,"code_words_total_at_most_11":len(words),"pair_controls":pairs,
            "law_head":[{"length":l,"amin":a,"prefix_count":n,"mass":str(m),"survival":str(s)}
                        for l,a,n,m,s in rows[:12]],
            "tail_length_120":str(rows[-1][-1]),
            "adaptive_length_cutoffs_epsilon_1_over_16":schedule,
            "ordered_frames":{"1132":uv,"2113":vu},
            "atomic_source_controls":atomic,
            "repeat_113_addresses":repeated}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
