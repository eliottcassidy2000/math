import random
def E3(S):
    Sset=set(S); return sum(1 for a in S for b in S if a+b in Sset)
print("E3 (ordered Schur triples a+b=c) extremal check among 13-sets:")
print("  AP {1..13}:", E3(list(range(1,14))), " 2*AP:", E3([2*i for i in range(1,14)]), " {0..12}:", E3(list(range(0,13))))
best=0;bestS=None;rng=random.Random(3)
for _ in range(300000):
    sp=rng.randint(12,80); S=sorted(rng.sample(range(0,sp+1),13)); e=E3(S)
    if e>best:best=e;bestS=S
print(f"  adversarial max E3 (300k random, spread<=80): {best} at {bestS}")
# hill-climb
S=sorted(rng.sample(range(0,50),13))
for _ in range(4000):
    i=rng.randrange(13);cand=rng.randint(0,50)
    if cand in S:continue
    S2=sorted(set(S)-{S[i]}|{cand})
    if len(S2)==13 and E3(S2)>E3(S):S=S2
print(f"  hill-climb max E3: {E3(S)} at {S}  diffs={[S[i+1]-S[i] for i in range(12)]}")
# is {0..12} (AP from 0) higher? sum-closure loves 0
print(f"  {{0..12}} check: E3={E3(list(range(0,13)))}  (0 is a Schur-identity a+0=a)")
