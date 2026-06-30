"""
WHY the H-maximizer is regular / doubly-regular / Paley -- test the theory.
Prediction: n=19 (prime, 3 mod4) -> Paley wins (doubly-regular exists & is Paley).
            n=17 (prime, 1 mod4) -> NOT Paley (no doubly-regular tournament).
Also: is each maximizer (a) REGULAR, (b) DOUBLY-REGULAR? (the structural forcing).
"""
from itertools import combinations

def Hcirc(n, out):  # out = set of difference 'd' with i->i+d
    adj=[0]*n
    for i in range(n):
        for d in out: adj[i]|=1<<((i+d)%n)
    full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]; dp[1][0]=1
    for mask in range(1<<n):
        if not mask&1: continue
        row=dp[mask]
        for last in range(n):
            c=row[last]
            if not c or not (mask>>last)&1: continue
            avail=adj[last] & ~mask
            while avail:
                nb=avail&-avail; nxt=nb.bit_length()-1
                dp[mask|nb][nxt]+=c; avail^=nb
    return n*sum(dp[full])

def out_from_S(n,S):  # S subset of 1..(n-1)/2 chosen 'low'
    return set(k if k in S else n-k for k in range(1,(n-1)//2+1))

def is_regular(n,out):  # circulant is always regular; trivially true (out-deg=(n-1)/2)
    return len(out)==(n-1)//2

def double_reg_defect(n,out):
    """doubly regular <=> for every nonzero diff d, #{(a,b) in out^2: a-b=d} is constant.
       return (min,max) of that count over nonzero d; constant (defect 0) <=> doubly regular."""
    cnt={d:0 for d in range(1,n)}
    for a in out:
        for b in out:
            d=(a-b)%n
            if d: cnt[d]+=1
    vals=list(cnt.values())
    return min(vals),max(vals)

def QRset(n): return set((x*x)%n for x in range(1,n))

print("=== n=17 full circulant search (prime, 1 mod4): expect NOT Paley ===")
n=17; half=8; best=-1; bo=None
for r in range(half+1):
    for S in combinations(range(1,half+1),r):
        out=out_from_S(n,set(S)); h=Hcirc(n,out)
        if h>best: best=h; bo=out
qr=QRset(17)
print(f"  maxH(circ,n=17)={best}  maximizer out-set={sorted(bo)}")
print(f"  QR mod 17={sorted(qr)}  is maximizer Paley/QR? {bo==qr or bo==set(range(1,17))-qr}")
print(f"  maximizer double-reg defect (min,max over diffs)={double_reg_defect(17,bo)}  (equal=>doubly regular)")
print(f"  H(Paley_17)={Hcirc(17,qr)}")

print("\n=== n=19 (prime, 3 mod4): PREDICTION Paley wins ===")
n=19; qr19=QRset(19)
hP=Hcirc(19,qr19)
hRot=Hcirc(19,set(range(1,10)))  # rotational i beats next 9
# a couple of other circulants for comparison
import random
others=[]
for seed in range(4):
    random.seed(seed); S=set(random.sample(range(1,10), random.randint(0,9)))
    others.append((sorted(out_from_S(19,S)), Hcirc(19,out_from_S(19,S))))
print(f"  H(Paley_19)={hP}   double-reg defect={double_reg_defect(19,qr19)}")
print(f"  H(rotational_19)={hRot}")
print(f"  random circulants H: {[h for _,h in others]}")
print(f"  Paley beats rotational & all sampled? {hP>hRot and all(hP>h for _,h in others)}")

print("\n=== double-regularity of the known maximizers (the structural test) ===")
maxers={7:out_from_S(7,set()),11:out_from_S(11,set()),13:out_from_S(13,{1,2,3,4,5,6}),
        15:out_from_S(15,{1,2,3,4,5,6,7})}
# use the actual maximizer out-sets found earlier; recompute cleanly:
known_out={7:{3,5,6},11:{2,6,7,8,10},13:{7,8,9,10,11,12},15:{8,9,10,11,12,13,14}}
for n,out in known_out.items():
    dr=double_reg_defect(n,out)
    print(f"  n={n} ({n%4} mod4) maximizer: doubly-regular? {dr[0]==dr[1]} (defect {dr})  [n=3mod4 should be DR; else not]")
