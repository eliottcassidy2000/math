import random
from collections import Counter
def E(S):
    r = Counter(a+b for a in S for b in S)
    return sum(c*c for c in r.values())
# AP and dilates
print("AP {1..13}:", E(list(range(1,14))))
print("2*AP:", E([2*i for i in range(1,14)]))
print("AP {0..12}:", E(list(range(0,13))))
# adversarial max over 13-distinct-int sets, spread<=60
best=0; bestS=None
rng=random.Random(7)
for _ in range(200000):
    spread=rng.randint(12,60)
    S=sorted(rng.sample(range(0,spread+1),13))
    e=E(S)
    if e>best: best=e; bestS=S
print(f"adversarial max E over 200k random 13-sets (spread<=60): {best}  at {bestS}")
# hill-climb from a random start toward max E
S=sorted(rng.sample(range(0,40),13))
for _ in range(3000):
    i=rng.randrange(13); old=S[i]
    cand=rng.randint(0,40)
    if cand in S: continue
    S2=sorted(set(S)-{old}|{cand})
    if len(S2)==13 and E(S2)>E(S): S=S2
print(f"hill-climb max E: {E(S)}  at {S}  (is it an AP? diffs={[S[i+1]-S[i] for i in range(12)]})")

print("\n--- Freiman sumset connection: max E <=> min |S+S| = 2n-1 (equality iff AP) ---")
def sumset(S): return len(set(a+b for a in S for b in S))
for name,S in [("AP {1..13}", list(range(1,14))), ("2*AP", [2*i for i in range(1,14)]),
               ("near-AP {1..12}+{20}", list(range(1,13))+[20]),
               ("Sidon", [1,2,5,11,22,33,45,60,78,95,110,130,150]),
               ("Fib-ish", [1,2,3,5,8,13,21,34,55,89,144,233,377])]:
    print(f"  {name:24} |S+S|={sumset(S):4}  E={E(S):5}  (2n-1={2*13-1})")
