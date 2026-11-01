"""A fixed 32-fraction evaluator for the full one-apex tournament family.

Word b[0]...b[m-1]: vertices x,v0,...,v(m-1), vi->vj for i<j,
x->vi iff b[i]=1. E(b)=2**m+sum(b[i]*2**i). Input 2*19*5**E.
H = exponent of 2 at halt, total directed odd cycles = exponent of 3.
"""
import json
import math
from pathlib import Path

DATA = {'H':2, 'C':3, 'N':5, 'Q':7, 'T':11, 'L':13, 'U':17}
STATES = {'read':19,'read_guard':23,'odd':29,'even':31,
          'restore1':37,'restore1_guard':41,'restore0':43,'restore0_guard':47,
          't1':53,'t1_guard':59,'l1':61,'l1_guard':67,'u1':71,'u1_guard':73,
          't0':79,'t0_guard':83,'l0':89,'l0_guard':97,'u0':101,'u0_guard':103,
          'halt':107,'error':109}
# These rows are in first-eligible order. All controls are distinct primes.
FRACTIONS = (
    (7*23, 5*5*19), (19,23), (29,5*19), (31,19),
    (5*37,7*29), (107,29), (5*43,7*31), (109,31),
    (5*41,7*37), (37,41), (53,37),
    (5*47,7*43), (43,47), (79,43),
    (17*17*59,11*53), (53,59), (61,53),
    (17*67,13*61), (61,67), (71,61),
    (11*73,17*71), (71,73), (13*19,71),
    (2*2*3*17*17*83,11*79), (79,83), (89,79),
    (2*2*3*17*97,13*89), (89,97), (101,89),
    (11*103,17*101), (101,103), (19,101),
)
PRIMES = tuple(DATA.values()) + tuple(STATES.values())

def require(c,msg):
    if not c: raise RuntimeError(msg)

def encode(word):
    require(all(type(b) is int and b in (0,1) for b in word),'binary word required')
    return 2**len(word)+sum(b*2**j for j,b in enumerate(word))

def expected(word):
    c=t=l=0
    for b in word:
        if b==0: c+=t+l
        t=2*t+l
        l=b
    return 1+2*c,c

def factor(x):
    out={}
    for p in PRIMES:
        while x%p==0:
            out[p]=out.get(p,0)+1
            x//=p
    require(x==1,'foreign prime')
    return out

ROWS=tuple((factor(a),factor(b)) for a,b in FRACTIONS)

def run_exponents(word, limit=10_000_000, rows=ROWS):
    exps={p:0 for p in PRIMES}
    exps[2]=1;exps[5]=encode(word);exps[19]=1
    steps=0
    while steps<limit:
        for pos,neg in rows:
            if all(exps[p]>=e for p,e in neg.items()):
                for p,e in neg.items(): exps[p]-=e
                for p,e in pos.items(): exps[p]+=e
                break
        else:
            require(exps[107]==1 and exps[109]==0,'not a successful halt')
            require(sum(exps[p] for p in STATES.values())==1,'non-one-hot terminal')
            require(all(exps[p]==0 for p in (5,7,17)),'nonempty disposable work register')
            return {'H':exps[2],'odd_cycles':exps[3],'steps':steps,'T':exps[11],'last':exps[13]}
        require(sum(exps[p] for p in STATES.values())==1,'control invariant violated')
        steps+=1
    raise RuntimeError('step bound reached')

def run_integer(word,limit=10_000_000):
    value=2*19*5**encode(word);steps=0;maxbits=value.bit_length()
    while steps<limit:
        for a,b in FRACTIONS:
            if value%b==0:
                value=value//b*a
                break
        else:
            e=factor(value)
            require(e.get(107,0)==1 and e.get(109,0)==0,'integer unsuccessful halt')
            require(set(e)<=set((2,3,11,13,107)),'integer bad terminal shape')
            return {'H':e.get(2,0),'odd_cycles':e.get(3,0),'steps':steps,'max_integer_bits':maxbits}
        maxbits=max(maxbits,value.bit_length());steps+=1
    raise RuntimeError('literal integer step bound reached')

def export():
    root=Path(__file__).resolve().parent
    for a,b in FRACTIONS: require(math.gcd(a,b)==1,'cancelled control')
    (root/'one_apex.fractions.txt').write_text(''.join(f'{a}/{b}\n' for a,b in FRACTIONS))
    info={'data_primes':DATA,'state_primes':STATES,'fractions':len(FRACTIONS),
          'input':'2*19*5**E(word)','halt_prime':107,'H_prime':2,'odd_cycle_prime':3}
    (root/'ONE_APEX_PROGRAM.json').write_text(json.dumps(info,indent=2,sort_keys=True)+'\n')
    print(json.dumps(info,indent=2,sort_keys=True))

if __name__=='__main__': export()
