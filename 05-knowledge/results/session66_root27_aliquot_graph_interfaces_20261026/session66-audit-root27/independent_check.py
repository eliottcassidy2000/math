"""Independent exact audit: native guard is enforced before any state mutation.
No producer code imports. Universe: all residues through 12 bits; all balanced
words through 15 bits (five internal vertices); finite marked and refill runs.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import argparse, hashlib, json, subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
parser=argparse.ArgumentParser(description=__doc__)
producer_default=next((p for p in (ROOT/'session66-root27',ROOT/'root27') if p.is_dir()),ROOT/'session66-root27')
parser.add_argument('--producer',type=Path,default=producer_default)
parser.add_argument('--repo',type=Path,default=ROOT/'session66-repo')
mode=parser.add_mutually_exclusive_group()
mode.add_argument('--offline',action='store_true',help='Verify sealed source bytes and Git blob hashes only (default).')
mode.add_argument('--verify-git-pins',action='store_true',help='Additionally inspect each historical tree in --repo.')
parser.add_argument('--pin-check-only',action='store_true',help='Emit a standalone source-pin verification receipt.')
args=parser.parse_args()
PRODUCER=args.producer.resolve()
REPO=args.repo.resolve()
COUNTS=Counter()

def need(x,why):
    if not x: raise ValueError(why)

def valuation(x):
    x=F(x)
    if x==0: return None
    def vi(a):
        a=abs(a); return (a & -a).bit_length()-1
    return vi(x.numerator)-vi(x.denominator)

def parity(z):
    z=F(z); need(z.denominator%2==1,'odd denominator')
    return z.numerator%2

def step(z,q):
    need(q>=1,'positive credit')
    b=parity(z)
    if q==1 and b==0: return None
    return (3*z+b)/2, q+3*b-1, b

def replay(z,q,n):
    z=F(z); word=''; states=[(z,q)]
    for _ in range(n):
        out=step(z,q); need(out is not None,'native blocked before requested endpoint')
        z,q,b=out; word+=str(b); states.append((z,q))
    return word,states

def raw_bits(z,n):
    z=F(z); bits=[]
    for _ in range(n):
        b=parity(z); bits.append(str(b)); z=(3*z+b)/2
    return ''.join(bits)

@lru_cache(None)
def data(w):
    n=len(w)
    c=sum((1<<i)*3**(n-1-i) for i,b in enumerate(w) if b=='1')
    return n,3**n,1<<n,c,3**n-(1<<n)

def affine(w,z):
    n,a,b,c,d=data(w); return (a*F(z)+c)/b

def center(w):
    n,a,b,c,d=data(w); return -F(c,d)

def depth(w,z): return valuation(F(z)-center(w))

def bitlift(w):
    r=0
    for i,b in enumerate(w):
        if raw_bits(r,i+1)[-1]!=b: r+=1<<i
    need(raw_bits(r,len(w))==w,'bitlift literal parity')
    return r

def balanced(w):
    c=0
    for b in w:
        c+=2 if b=='1' else -1
        if c<0:return False
    return c==0

def prim(w):
    return next(w[:k] for k in range(1,len(w)+1) if len(w)%k==0 and w[:k]*(len(w)//k)==w)

def pull(w,p):
    n,a,b,c,d=data(w); return (b*center(p)-c)/a

def source_pins():
    pins=json.loads((PRODUCER/'SOURCE_PINS.json').read_text())['sources']
    for pin in pins:
        payload=(PRODUCER/pin['snapshot']).read_bytes()
        need(hashlib.sha256(payload).hexdigest()==pin['sha256'],'pin SHA256')
        blob=hashlib.sha1(b'blob '+str(len(payload)).encode()+b'\0'+payload).hexdigest()
        need(blob==pin['git_blob'],'pin git object digest')
        if args.verify_git_pins:
            committed=subprocess.check_output(['git','-C',str(REPO),'show',pin['commit']+':'+pin['repo_path']])
            need(payload==committed,'pin committed bytes')
            named=subprocess.check_output(['git','-C',str(REPO),'rev-parse',pin['commit']+':'+pin['repo_path']],text=True).strip()
            need(named==blob,'pin tree blob identity')
        COUNTS['source_pins_verified']+=1
    for name in (['05-knowledge/results/root27_machine_bridges_20261026/anchor/README.md','05-knowledge/results/lca_ports_reverse_20261025/arithmetic/ANCESTRY.md'] if args.verify_git_pins else []):
        old=subprocess.check_output(['git','-C',str(REPO),'show','63f52936:'+name])
        new=subprocess.check_output(['git','-C',str(REPO),'show','2807ee83:'+name])
        need(old==new,'inherited root anchor/ancestry unchanged in incoming repair')
    return pins

def main():
    source_pins()
    blocks=[]
    # Generate words by literal residue itineraries, not recursive tree grammar.
    for n in range(13):
        seen=set(); carries=set()
        for r in range(1<<n):
            word=raw_bits(r,n); seen.add(word)
            nn,a,b,c,d=data(word); carries.add(c)
            need((a*r+c)%b==0,'literal residue implies carry cylinder')
            need((-c*pow(a,-1,b))%b==r,'carry recovers literal source residue')
            COUNTS['all_residue_word_controls']+=1
        need(len(seen)==len(carries)==1<<n,'carry and itinerary injectivity')
    for n in (3,6,9,12,15):
        for word in map(''.join,product('01',repeat=n)):
            if balanced(word): blocks.append(word)
    need(Counter(map(len,blocks))=={3:1,6:3,9:12,12:55,15:273},'independent balanced census')
    common=0; regimes=Counter()
    for p in blocks:
        lp,ap,bp,cp,dp=data(p)
        for q in blocks:
            lq,aq,bq,cq,dq=data(q)
            delta=cq*dp-cp*dq
            need(affine(q,affine(p,0))-affine(p,affine(q,0))==-F(delta,bp*bq),'ordered commutator constant sign')
            need(affine(q,affine(p,1))-affine(p,affine(q,1))==-F(delta,bp*bq),'ordered commutator linear part')
            eq=delta==0
            need(eq==(center(p)==center(q))==(p+q==q+p)==(prim(p)==prim(q)),'primitive/commutation/center iff')
            common+=eq; COUNTS['balanced_ordered_pairs']+=1
            z=(-cp*pow(ap,-1,bp))%bp+bp
            out,states=replay(z,1,lp); y=states[-1][0]
            need(out==p and states[-1][1]==1,'balanced native execution and q preservation')
            need(bp*dp*(dq*y+cq)==ap*dq*(dp*z+cp)+bp*delta,'signed cocycle')
            u=depth(p,z)-lp; hq=depth(q,y)
            need(depth(p,y)==u,'exact own-fuel consumption')
            if eq: need(hq==u,'zero-defect branch'); regimes['zero']+=1
            else:
                t=valuation(delta)
                if u==t: need(hq>=t+1,'equal-valuation cancellation'); regimes['equal']+=1
                else: need(hq==min(u,t),'unequal valuation minimum'); regimes['unequal']+=1
            # Pullbacks tested with nontrivial arbitrary word; no execution claim here.
            for prefix in ('','111','111010'):
                if delta:
                    need(valuation(pull(prefix,q)-pull(prefix,p))==len(prefix)+valuation(delta),'pair pullback Eq6')
                    COUNTS['pair_pullback_separations']+=1
    word,states=replay(7,1,110)
    need([(int(z),q) for z,q in states[:5]]==[(7,1),(11,3),(17,5),(26,7),(39,6)],'root7 to odd-x root27 mapping')
    for n,(z,credit) in enumerate(states):
        prefix=word[:n]
        for p in blocks:
            need(valuation(F(7)-pull(prefix,p))==n+depth(p,z),'original source Eq5')
            need(pull(prefix+p,p)==pull(prefix,p),'retained center source transport')
            COUNTS['marked_prefix_center_controls']+=1
    z,q0=states[86];y,q1=states[89]
    need((int(z),q0,int(y),q1,word[86:105])==(10952205219944621,47,36963692617313097,47,'1001010001011101000'),'marked macro86 states bits guards')
    need((depth('100',z),depth('100',y),depth('101000',y))==(5,2,9),'marked seven-bit refill')
    np=19*z+9; delta=684
    need((valuation(27*665*np),valuation(8*delta),valuation(27*665*np+8*delta))==(5,5,12),'marked exact equal cancellation')
    # Reconstruct old odd-x returns only by the native scalar macros.
    z,q=states[3]; old=None; returned=0
    def v3(x):
        x=abs(int(x)); need(x>0,'nonzero v3'); k=0
        while x%3==0:x//=3;k+=1
        return k
    for j in range(20):
        x=z+1;a=q-1; count=0
        while True:
            out=step(z,q);need(out is not None,'old native guard'); z,q,b=out;count+=1
            if parity(z+1)==1:break
        returned+=v3(z+1)
        need(q-1==6+2*returned-(j+1),'returned v3 indexing excludes initial 27')
        need(q-1==a+2*count-3,'native odd-x cumulative credit')
        if j==19:old={'index':j,'x':int(x),'a':a,'k':count,'v3_x_minus_1':v3(x-1),'charge_before':a+v3(x-1),'charge_after':q-1+v3(z)}
    need(old==json.loads((PRODUCER/'normal.json').read_text())['old_root27_charge_hostile'],'old hostile independent physical cuts')
    family=[]
    for m in [1,2,3,4,5,8,16,32,64]:
        w='100'+'101000'*m+'0'; r=bitlift(w); z=r+(1<<len(w))
        need(z>=3,'positive physical arithmetic source')
        native,st=replay(z,1,len(w)-1)
        need(native==w[:-1] and st[-1][1]==1 and parity(st[-1][0])==0,'native prefix ends before blocked leaf')
        need(step(*st[-1]) is None,'blocked leaf not executed')
        y=st[3][0]
        need((depth('100',z),depth('101000',y))==(5,6*m),'all-m refill finite controls')
        family.append({'m':m,'seed':str(z),'executed':len(native),'blocked_state':str(st[-1][0]),'credit':st[-1][1]})
        COUNTS['family_native_macros']+=len(native);COUNTS['blocked_leaf_rejections']+=1
    # Same center does not imply same clocked source transform.
    need(center('100')==center('100100'),'same center positive control')
    need(affine('100',7)!=affine('100100',7),'same center different transition hostile')
    need(center('100')!=center('101000'),'two-center finite closure obstruction')
    need(step(F(4),2)==(F(6),1,0) and step(F(6),1) is None,
         'initial-credit-two lifetime-one hostile')
    need(all(depth(p,4)==0 for p in blocks),'q=1 lifetime lemma cannot ignore credit')
    # Lemma B: every finite observed legal prefix has a nonempty balanced completion.
    for seed in range(1,257):
        z=F(seed);q=1;prefix=''; h=None
        for n in range(81):
            if prefix:
                p=prefix+'0'*(q-1)
                need(balanced(p),'legal prefix admits balanced zero completion')
                need(depth(p,seed)>=n,'completion center matches legal prefix')
                COUNTS['legal_prefix_balanced_completions']+=1
            out=step(z,q)
            if out is None:
                h=n;break
            z,q,b=out;prefix+=str(b)
        if h is not None:
            # Finite bank cannot exceed exact legal lifetime.
            vals=[depth(p,seed) for p in blocks]
            need(max(vals)<=h,'periodic balanced centers cannot match blocked symbol')
            if h==0:need(max(vals)==0,'even-source lifetime zero boundary')
            else:need(balanced(prefix) and depth(prefix,seed)==h,'finite lifetime attained by its whole legal block')
            COUNTS['halting_seed_lifetime_controls']+=1
        else:COUNTS['not_halted_by_prefix80_seeds']+=1
    # Rational periodic centers really replay forever at q=1 (three periods controls).
    for p in blocks:
        r=center(p); word3,st=replay(r,1,3*len(p))
        need(word3==p*3 and st[-1]==(r,1),'rational periodic balanced-center legal orbit')
        COUNTS['rational_periodic_center_controls']+=1
    return {'status':'PASS','counts':dict(sorted(COUNTS.items())),'common_center_pairs_five_internal_vertices':common,'valuation_regimes':dict(regimes),'refill_family':family,'lemmas_A_B':'PROVED symbolically in CENTER_BANK_APPENDIX.md; checks are finite controls only','open':'Root27 all-prefix survival and a global source-specific rank remain open.'}

if __name__=='__main__':
    if args.pin_check_only:
        pins=source_pins()
        result={'status':'PASS','verification_mode':'HISTORICAL_GIT_TREE_AND_SEALED_BYTES' if args.verify_git_pins else 'OFFLINE_SEALED_SHA256_AND_GIT_BLOB',
                'sources':pins,'counts':dict(COUNTS),
                'anchor_and_ancestry_unchanged_63f52936_to_2807ee83':True if args.verify_git_pins else 'not checked offline'}
    else:
        result=main()
    print(json.dumps(result,indent=2,sort_keys=True))
