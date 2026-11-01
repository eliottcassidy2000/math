"""Producer tests: exact scopes, hostile controls, and optimization-safe checks."""
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
from build_tm import Builder,delta,validate,expand
from simulate_tm import encode,run
from square_fractran import FRACTIONS,expected_steps,run as run_fraction
HERE=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def predicate(n,word,mode):
    if n<1 or (mode=='cycle' and n<3) or sorted(word)!=list(range(1,n+1)):
        return False
    pairs=list(zip(word,word[1:]))
    if mode=='cycle': pairs.append((word[-1],word[0]))
    return all(math.isqrt(a+b)**2==a+b for a,b in pairs)

def cost(n,word,mode):
    roots=sum(math.isqrt(a+b) for a,b in zip(word,word[1:]))
    if mode=='path': return 5*n*n+14*n+3-2*word[0]-5*word[-1]+6*roots
    roots+=math.isqrt(word[-1]+word[0])
    return 5*n*n+14*n+11+3*word[0]+6*roots

def check():
    receipt={'status':'FINITE-EXACT','fractran':{},'tm':{},'shells':{}}
    for n in range(1025):
        a,steps=run_fraction(n)
        k=math.isqrt(n); square=(k*k==n)
        expected=({31:1} if n==0 else {31:1,5:2*k-1}) if square else {37:1,3:2*k-(n-k*k),5:n-k*k+1}
        expected={p:e for p,e in expected.items() if e}
        need(a==expected,('fraction output',n,a,expected))
        need(steps==expected_steps(n),('fraction steps',n,steps))
        if n<=256: need(run_fraction(n,literal=True)==(a,steps),('literal discrepancy',n))
    need(all(math.gcd(a,b)==1 for a,b in FRACTIONS),'unreduced fraction')
    # A literal first-eligible mutation: placing the initial zero fallback first
    # accepts nonsquare input2 without consuming its exponent.
    mutant=FRACTIONS[:]; mutant[0],mutant[1]=mutant[1],mutant[0]
    value=7*2**2
    for a,b in mutant:
        if value%b==0: value=value//b*a; break
    need(value==31*2**2 and all(value%b for _,b in mutant),'priority hostile')
    receipt['fractran']={'exponent_inputs':'all n=0..1024','literal_inputs':'all n=0..256','checks':'exact terminal exponent vector and exact step formula','fraction_count':len(FRACTIONS),'priority_swap_n2_final_integer':124}
    tables={mode:json.loads((HERE/(mode+'.tm.json')).read_text()) for mode in ['path','cycle']}
    frozen={item['mode']:item for item in json.loads((HERE/'TM_BUILD.json').read_text())}
    for mode,machine in tables.items():
        need(machine==Builder(mode).build(),('builder mismatch',mode)); validate(machine)
        if not (HERE/(mode+'.tm.tsv.gz')).exists(): expand(machine,HERE/(mode+'.tm.tsv.gz'))
        h=hashlib.sha256(); count=0; seen=set()
        with gzip.open(HERE/(mode+'.tm.tsv.gz'),'rt') as f:
            for line in f:
                h.update(line.encode()); q,r,nxt,w,m=line.rstrip('\n').split('\t')
                need((q,r) not in seen,('duplicate expanded key',mode,q,r)); seen.add((q,r))
                need(delta(machine,q,r)==(nxt,tuple(w),tuple(m)),('expanded row mismatch',mode,q,r))
                count+=1
        need(count==len(machine['rules'])*math.prod(map(len,machine['alphabets'])),('expanded row count',mode))
        need(h.hexdigest()==frozen[mode]['expanded_sha256'],('frozen expanded hash',mode))
        need(count==frozen[mode]['concrete_delta_rows'],('frozen row count',mode))
        receipt['tm'][mode]={'states':len(machine['rules']),'explicit_rows':count,'expanded_sha256':h.hexdigest()}
    checks=0; accepts=0
    # Every bounded unary word, not just permutations. Zero is intentionally
    # encoded as an empty field to exercise a malformed-token rejection.
    for n in range(5):
        for length in range(n+2):
            for word in itertools.product(range(n+2),repeat=length):
                for mode in ['path','cycle']:
                    result=run(encode(n,word),mode,machine=tables[mode])
                    yes=predicate(n,word,mode)
                    need((result['state']=='ACCEPT')==yes,('bounded predicate',n,word,mode,result))
                    if yes:
                        need(result['steps']==cost(n,word,mode),'accepted cost')
                        accepts+=1
                    checks+=1
    malformed=0
    for length in range(7):
        for chars in itertools.product('^1#$',repeat=length):
            text=''.join(chars)
            # No valid positive path except the explicit length5 singleton can
            # fit this range; parse independently rather than assume rejection.
            yes=False
            if text.startswith('^') and text.endswith('$') and text.count('^')==text.count('$')==1:
                fields=text[1:-1].split('#')
                if len(fields)>=2 and all(field and set(field)=={'1'} for field in fields):
                    yes=predicate(len(fields[0]),[len(x) for x in fields[1:]],'path')
            for mode in ['path','cycle']:
                result=run(text,mode,machine=tables[mode])
                need((result['state']=='ACCEPT')==(yes and mode=='path'),('raw syntax',text,mode))
                malformed+=1
    w=json.loads((HERE/'witnesses.json').read_text())
    cases=[('path',15,w['path15'])]+[('cycle',c['n'],c['cycle']) for c in w['cycles']]
    positive=[]; mutations=0
    for mode,n,word in cases:
        need(predicate(n,word,mode),'bad supplied witness')
        result=run(encode(n,word),mode,machine=tables[mode])
        need(result['state']=='ACCEPT' and result['steps']==cost(n,word,mode),('positive',mode,n,result))
        roots=sum(math.isqrt(word[i]+word[(i+1)%n]) for i in range(n)) if mode=='cycle' else sum(math.isqrt(a+b) for a,b in zip(word,word[1:]))
        positive.append({'mode':mode,'n':n,'steps':result['steps'],'root_weight':roots})
        # Every cyclic gauge of each cycle checks the asserted start-coordinate
        # cost, not only one favorable normalization.
        if mode=='cycle':
            for reverse in [False,True]:
                oriented=word[::-1] if reverse else word
                for shift in range(n):
                    rotated=oriented[shift:]+oriented[:shift]
                    r=run(encode(n,rotated),mode,machine=tables[mode])
                    need(r['state']=='ACCEPT' and r['steps']==cost(n,rotated,mode),('cycle gauge',n,shift,reverse))
            need(roots%2==0,'root-weight parity')
        bads=[word[:-1],word+[word[0]],[word[0]]+word[:-1],word[:-1]+[n+1],[]]
        for bad in bads:
            need(run(encode(n,bad),mode,machine=tables[mode])['state']=='REJECT',('hostile',mode,n,bad)); mutations+=1
    hostile=[1,8]*4
    need(len(hostile)==8 and sum(hostile)==36 and all(math.isqrt(hostile[i]+hostile[(i+1)%8])**2==hostile[i]+hostile[(i+1)%8] for i in range(8)),'hostile premise')
    need(run(encode(8,hostile),'cycle',machine=tables['cycle'])['state']=='REJECT','named-mark hostile')
    need(run(encode(15,w['path15']),'cycle',machine=tables['cycle'])['state']=='REJECT','path/cycle closing hostile')
    receipt['tm'].update({'bounded_words':{'n':'0..4','length':'0..n+1','digits':'0..n+1','runs':checks,'accepts':accepts},'raw_strings':{'alphabet':'^1#$','length':'0..6','runs':malformed},'positive_witnesses':positive,'cycle_gauge_runs':2*sum(c['n'] for c in w['cycles']),'mutations':mutations,'global_sidecar_hostile':hostile})
    total=0
    for n in range(1,1001):
        direct={v for v in range(1,n) if math.isqrt(n+v)**2==n+v}
        shell={k*k-n for k in range(math.isqrt(n)+1,math.isqrt(2*n-1)+1)}
        need(direct==shell and len(shell)==math.isqrt(2*n-1)-math.isqrt(n),('shell',n))
        total+=len(shell)
    receipt['shells']={'arrivals':'all n=1..1000','edge_total_G1000':total}
    return receipt

if __name__=='__main__': print(json.dumps(check(),sort_keys=True,indent=2))
