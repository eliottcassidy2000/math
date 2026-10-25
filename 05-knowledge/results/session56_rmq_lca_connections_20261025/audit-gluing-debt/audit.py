from fractions import Fraction as F
from itertools import product, combinations
from collections import Counter
from pathlib import Path
import json, hashlib, math, subprocess, sys, os
ROOT=Path(__file__).resolve().parent
G=ROOT.parent/'correlated-gluing'; D=ROOT.parent/'debt-range'
counts=Counter()
def ck(x,label):
    counts[label]+=1
    if not x: raise RuntimeError(label)
for line in (G/'SHA256SUMS').read_text().splitlines():
    digest,name=line.split(None,1)
    ck(hashlib.sha256((G/name).read_bytes()).hexdigest()==digest,'gluing_packet_hash')
for line in (G/'SOURCE_PINS.txt').read_text().splitlines():
    digest,name=line.split(None,1)
    ck(hashlib.sha256((G/name).read_bytes()).hexdigest()==digest,'gluing_source_pin')
for name,digest in json.loads((D/'MANIFEST.json').read_text())['files'].items():
    ck(hashlib.sha256((D/name).read_bytes()).hexdigest()==digest,'debt_packet_hash')
for item in json.loads((D/'SOURCES.json').read_text())['sources']:
    ck(hashlib.sha256((D/item['snapshot']).read_bytes()).hexdigest()==item['sha256'],'debt_source_pin')
# Separate line-equation enumeration, no producer code or parametrization imported.
LINES=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
PAIRS=tuple(combinations(range(9),2))
def arrays(n):
    Q={x*x%n for x in range(n)}; out=set()
    for a,b,c,e in product(Q,repeat=4):
        target=(a+b+c)%n
        g=(target-c-e)%n; i=(target-a-e)%n
        d=(target-a-g)%n; f=(target-c-i)%n; h=(target-b-e)%n
        row=(a,b,c,d,e,f,g,h,i)
        if set(row)<=Q and all(sum(row[k] for k in L)%n==target for L in LINES):out.add(row)
    return sorted(out)
def mask(a):return sum(1<<k for k,(i,j) in enumerate(PAIRS) if a[i]==a[j])
mods=(2,5,3); local=[arrays(q) for q in mods]; full=arrays(30)
ck(list(map(len,local))==[8,7,2],'independent_local_cardinalities')
# CRT independently solved by search table on residues, not inverse formula.
crt={(x%2,x%5,x%3):x for x in range(30)}
joined={tuple(crt[tuple(a[i] for a in choice)] for i in range(9)) for choice in product(*local)}
ck(joined==set(full) and len(full)==112,'full_mod30_bijection')
ck(all(len(set(a))<9 for a in full),'no_genuine_distinct')
left=list(product(local[0],local[1])); right=list(product(local[1],local[2])); fake=0
for a,b in product(left,right):
    if mask(a[0])&mask(a[1])&mask(b[0])&mask(b[1])==0:fake+=1
ck(len(left)*len(right)==784 and fake==128,'fabricated_distinct_128')
ck(sum(c*c for c in Counter(map(mask,local[1])).values())==17,'naive_mask_pairs_17')
cert=json.loads((G/'certificate.json').read_text())
ck(Counter(map(mask,full))=={int(k,16):v for k,v in cert['hostile']['histogram'].items()},'full_histogram')
# Exhaustive binary chain relations, independent full-assignment filtering vs messages.
pairs=list(product(range(2),repeat=2))
rels=[{p for k,p in enumerate(pairs) if bits>>k&1} for bits in range(16)]
w=[[F(-2,3),F(0)],[F(1,2),F(-3)],[F(0),F(4,5)],[F(2),F(-1,7)]]
for A,B,C in product(rels,repeat=3):
    global_rows={x for x in product(range(2),repeat=4) if x[:2] in A and x[1:3] in B and x[2:] in C}
    exact=sum((math.prod(w[i][x[i]] for i in range(4)) for x in global_rows),F(0))
    # Own a at left, b/c at middle, d at right, including signed and zero values.
    ml={b:sum((w[0][a] for a,bb in A if bb==b),F(0)) for b in range(2)}
    mr={c:sum((w[3][d] for cc,d in C if cc==c),F(0)) for c in range(2)}
    got=sum((ml[b]*w[1][b]*w[2][c]*mr[c] for b,c in B),F(0))
    ck(got==exact,'signed_zero_tree_dp')
    reduced=[set(A),set(B),set(C)]
    for source,target in [(0,1),(2,1),(1,0),(1,2)]:
        # Root middle, bottom-up then top-down, coordinate-aware intersections.
        shared=max(source,target)
        si=shared-source; ti=shared-target
        reduced[target]={r for r in reduced[target] if any(s[si]==r[ti] for s in reduced[source])}
    ck(all(reduced[i]=={x[i:i+2] for x in global_rows} for i in range(3)),'semijoin_projection_exact')
for weights in product((-2,-1,0,1,2),repeat=4):
    pushed=Counter()
    for s,v in enumerate(weights):pushed[s&s]+=F(v,3)
    ck(all(pushed[s]==F(weights[s],3) for s in range(4)),'signed_diagonal')
f={0:F(1,2),1:F(1,2)}
ind=Counter()
for a,b in product(f,repeat=2):ind[a&b]+=f[a]*f[b]
ck(ind!=f,'diagonal_not_tensor')
# Lossy projections and cyclic extension obstruction checked by brute assignment.
even={x for x in product(range(2),repeat=3) if sum(x)%2==0}
projections=({x[:2] for x in even},{x[1:] for x in even})
ck(len(even)==4 and sum(x[:2] in projections[0] and x[1:] in projections[1] for x in product(range(2),repeat=3))==8,'lossy_parity')
ck(not any(a==b and b==c and a!=c for a,b,c in product(range(2),repeat=3)),'cyclic_unsatisfiable')
# Independent affine algebra: derive thresholds as preimages, compare literal evaluation.
I=(F(1),F(0),F(0))
def compose(A,B):
    l,d,r=A;m,e,s=B
    return m*l,m*d+e,max(r,(d+s)/l)
def summary(word):
    l=F(1);b=F(0);bounds=[F(0)]
    for scale,debt in word:
        l*=scale;b=scale*b-debt;bounds.append(-b/l)
    return l,-b,max(bounds)
def trajectory(word,x):
    out=[x]
    for a,b in word:x=a*x-b;out.append(x)
    return out
alphabet=list(product((F(1,3),F(1),F(3,2)),(F(-2),F(0),F(3,2))))
for k in range(5):
    for word in product(alphabet,repeat=k):
        s=summary(word);r=s[2]
        for cut in range(k+1):ck(compose(summary(word[:cut]),summary(word[cut:]))==s,'affine_all_cuts')
        for x in (r-F(1,17),r,r+F(1,17),F(math.ceil(r)),F(math.ceil(r)-1)):
            ck((min(trajectory(word,x))>=0)==(x>=r),'affine_prefix_threshold')
        for modulus,residue in ((3,1),(5,4)):
            low=residue+modulus*math.ceil((r-residue)/modulus)
            ck(low>=r and low-modulus<r and low%modulus==residue,'guarded_integer_ceiling')
# Literal actual odd Collatz, full small residue universes and source pullback.
def step(n):
    x=3*n+1;a=0
    while x%2==0:x//=2;a+=1
    return x,a
def native(n,k):
    states=[n];vals=[]
    for _ in range(k):n,a=step(n);states.append(n);vals.append(a)
    return tuple(vals),states
def cylinder(word):
    power=1;carry=0
    for a in word:carry=3*carry+power;power*=2**a
    M=2*power;r=((power-carry)*pow(3**len(word),-1,M))%M
    return r,M
for k in range(1,5):
    for word in product((1,2,3),repeat=k):
        r,M=cylinder(word)
        for n in range(1,M,2):
            vals,states=native(n,k)
            ck((vals==word)==(n==r),'exhaustive_exact_native_cylinder')
        for lift in range(3):
            n=r+lift*M;vals,states=native(n,k)
            for cut in range(k+1):
                ra,ma=cylinder(word[:cut]);rb,mb=cylinder(word[cut:])
                ck(n%ma==ra and states[cut]%mb==rb,'native_guard_pullback')
            l=F(1);P=Dd=F(0)
            for a in word:
                scale=F(3,2**a);p=F(int(a==1));d=F(1)-F(4,2**a) if a>=3 else F(0)
                l*=scale;P=scale*P+p;Dd=scale*Dd+d
            ck(l*(n-1)+P-Dd==states[-1]-1,'native_positive_charge_ledger')
            ck(l*(n-1)-Dd<=states[-1]-1,'native_lower_envelope')
ck(summary(((F(1),F(2)),(F(1),F(-2))))==(F(1),F(0),F(2)),'endpoint_cone_not_prefix')
ck(trajectory(((F(1,2),F(0)),),F(1))[-1]==F(1,2),'integer_input_not_integer_path')
ck(summary(((F(1),F(1)),)*3)[2]==3 and compose(summary(((F(1),F(1)),)*2),summary(((F(1),F(1)),)*2))[2]==4,'overlap_debt_doublecount')
# Replay producer checks in each interpreter mode, without writing in producer dirs.
env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
for folder,script,args in ((G,'verify.py',['--verify']),(D,'check.py',[])):
    normal=subprocess.run([sys.executable,str(folder/script),*args],capture_output=True,text=True,env=env,check=True).stdout
    opt=subprocess.run([sys.executable,'-O',str(folder/script),*args],capture_output=True,text=True,env=env,check=True).stdout
    ck(normal==opt,'producer_mode_agreement')
    expected=(folder/('replay.txt' if folder==G else 'normal.json')).read_text()
    ck(normal==expected,'producer_frozen_receipt_match')
print(json.dumps({'status':'PASS','checks':sum(counts.values()),'counts':dict(sorted(counts.items())),'mod30_arrays':len(full),'independent_pairs':len(left)*len(right),'fabricated_distinct':fake,'scope':'Independent finite audit; elementary all-size proofs reviewed separately. No coverage or novelty claim.'},indent=2,sort_keys=True))
