"""
LRC transfer test: does the witness-count dispersion CV^2_S(D(q,S)/q) over covering sets VANISH as q->inf
(like the metagraph CV^2(H)~2/n, full concentration) or PLATEAU (a persistent tight tail = the hard part)?
"""
import random, statistics as st
def D_over_q(S,q,n=14):
    cnt=0
    for a in range(1,q):
        if all(min((v*a)%q, q-(v*a)%q)*n>=q for v in S): cnt+=1
    return cnt/q
def is_cov(S): return all(any(s%qq==0 for s in S) for qq in range(2,15))
random.seed(5)
samples=[[1,2,3,4,5,6,7,8,9,10,11,13,84],[1,2,3,4,5,6,7,8,9,10,11,13,28],list(range(2,15))]
cc=0
while len(samples)<50 and cc<100000:
    cc+=1; S=sorted(random.sample(range(1,30),13))
    if is_cov(S): samples.append(S)
print("CV^2_S(D/q) over covering sets vs shell q (metagraph analog: CV^2(H)~2/n -> 0):")
print(f"{'q':>6} {'mean D/q':>10} {'var_S':>10} {'CV^2_S':>9} {'min D/q (tail)':>14}")
for q in [101,503,1009,2003,4001]:
    ds=[D_over_q(S,q) for S in samples]
    m=st.mean(ds); v=st.pvariance(ds); cv=v/m**2 if m>0 else 0
    print(f"{q:>6} {m:>10.5f} {v:>10.6f} {cv:>9.5f} {min(ds):>14.6f}")
print("\nIf CV^2_S PLATEAUS (not ->0): the covering ensemble has a persistent spread = the singular-series")
print("variance L(S); the near-tight sets are a real TAIL the metagraph (which fully concentrates) misses.")
