#!/usr/bin/env python3
"""Literal odd-word payments and a separate primitive-P annulus ledger.

First-exit discovery is bounded and labelled. ROOT suffixes are supplied data.
The annulus scan retains original sources, compatible residues, and full clocks.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO/'04-computation/experiments'))
import native_cancellation_105_20261018g as native

CHECKS = 0
CAP = 10000
W63 = (1,1,1,1,1,4,1,2,1,1,2,1,1,1,2,3,1,1,2,1,2,1,1,1,1,1,3,1,1,1,4,2,2,4)
ROOT_REGISTERS = (71,78,1,8,36,43)


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def integer(x, minimum=0):
    if type(x) is not int or x < minimum:
        raise ValueError('exact integer outside domain')
    return x


def odd(n):
    integer(n, 1)
    if not n & 1:
        raise ValueError('positive odd source required')
    return n


def word(w):
    if type(w) is not tuple:
        raise ValueError('immutable valuation word required')
    for a in w:
        integer(a, 1)
    return w


def v2(n):
    integer(n, 1)
    return (n & -n).bit_length()-1


def state(w):
    """Forward shifted slope, positive emission, reset debt, and raw A,Q,B."""
    word(w)
    lam, pos, debt = F(1), F(0), F(0)
    A, Q, B = 1, 1, 0
    for a in w:
        mu = F(3, 2**a)
        lam *= mu
        pos = mu*pos + (a == 1)
        debt = mu*debt + (1-F(4,2**a) if a >= 3 else 0)
        A, Q, B = 3*A, Q*2**a, 3*B+Q
    return lam, pos, debt, A, Q, B


def replay(n, w, root=False):
    odd(n); word(w)
    if type(root) is not bool:
        raise ValueError('exact ROOT flag required')
    x, states = n, [n]
    for a in w:
        if x == 1:
            raise ValueError('ROOT padding rejected')
        y = 3*x+1
        if v2(y) != a:
            raise ValueError('source is outside the native word')
        x = y >> a
        states.append(x)
    if root and x != 1:
        raise ValueError('supplied suffix does not reach ROOT')
    return tuple(states)


def first_exit(n):
    """Explicit bounded discovery; no cap failure is treated as a theorem."""
    odd(n)
    if n == 1:
        raise ValueError('ROOT has no strict first exit')
    x, w = n, []
    for _ in range(CAP):
        y=3*x+1; a=v2(y); x=y >> a; w.append(a)
        if x < n:
            return tuple(w)
    raise RuntimeError(('UNRESOLVED first exit', n))


def terminal_language(w):
    """{1,2}* or {1,2}* a 2*, a>=3, always a nonempty word."""
    word(w)
    if not w:
        return False
    resets = [i for i,a in enumerate(w) if a >= 3]
    return not resets or len(resets) == 1 and all(a == 2 for a in w[resets[0]+1:])


def cone_cut(w):
    word(w)
    for k in range(len(w)-1, -1, -1):
        lam, _, debt, *_ = state(w[:k])
        if lam >= 1 and debt <= lam-1 and terminal_language(w[k:]):
            return k
    return None


def paid(n, w, k):
    """Authenticate a supplied source, cut and actual descending word."""
    odd(n); word(w); integer(k)
    if not 0 <= k < len(w) or not terminal_language(w[k:]):
        raise ValueError('bad terminal-language cut')
    trajectory = replay(n,w)
    h=trajectory[-1]
    lam, _, debt, *_ = state(w[:k])
    if lam < 1 or debt > lam-1 or h >= n:
        raise ValueError('cone or source descent premise failed')
    mu, _, dv, *_ = state(w[k:])
    if not (mu < 1 and dv < 1):
        raise RuntimeError('proved terminal bounds failed')
    total_debt=state(w)[2]
    if not total_debt < n-h:
        raise RuntimeError('source payment failed')
    return h, total_debt


def discharge(n,w,k,suffix):
    h,_=paid(n,w,k)
    replay(h,suffix,root=True)
    whole=w+suffix
    replay(n,whole,root=True)
    return whole


def coarse63(n):
    """Infinite native ray; variable actual final reset, no ROOT conclusion."""
    odd(n)
    if n < 63 or n % (1 << 54) != 63:
        raise ValueError('outside the coarse63 ray')
    x=replay(n,W63[:-1])[-1]
    a=v2(3*x+1)
    if a < 4:
        raise RuntimeError('coarse final valuation failed')
    w=W63[:-1]+(a,)
    h,debt=paid(n,w,32)
    return {'source':n,'word':w,'endpoint':h,'last_valuation':a,
            'debt':str(debt),'gap':n-h}


def power4_phase(target, precision):
    """Unique J modulo3^(precision-1) with4^J=target modulo3^precision."""
    integer(target,1);integer(precision,1)
    if target%3 != 1:
        raise ValueError('target outside the principal 3-adic units')
    j,period=0,1
    for k in range(1,precision):
        modulus=3**(k+1)
        candidates=[j+i*period for i in range(3)
                    if pow(4,j+i*period,modulus)==target%modulus]
        if len(candidates)!=1:
            raise RuntimeError('unique ternary phase lift failed')
        j=candidates[0];period*=3
    return j,period


def root63_lift(index):
    """A source formula plus its actual first-ROOT word; no huge expansion."""
    integer(index)
    j0,period=power4_phase(61,34)
    # The canonical solution already exceeds the positivity cut J>=3.
    if j0 < 3:
        raise RuntimeError('canonical positive-height cut changed')
    j=j0+period*index
    return {'index':index,'J':j,'period':period,'head':W63[:-1],
            'word':W63[:-1]+(4+2*j,),
            'source_formula':'63 + 2^54 * (4^J - 61) / 3^34',
            'endpoint':1}


def root63_residue(index, modulus):
    integer(index);integer(modulus,1)
    j=root63_lift(index)['J'];denominator=3**34
    residue=pow(4,j,denominator*modulus)
    if (residue-61)%denominator:
        raise RuntimeError('full-denominator quotient is not integral')
    return (63+(1<<54)*((residue-61)//denominator))%modulus


def root63_registers(index):
    integer(index)
    r=ROOT_REGISTERS[index%6]
    return r,r%15,r%35


def rejects(fn):
    try:
        fn()
    except (ValueError,RuntimeError):
        check(True,'rejected hostile')
    else:
        check(False,'accepted hostile')


def annulus(n,routes,retain=False):
    """All sources and all finite arrival horizons, fixed inherited frontier."""
    integer(n,1)
    if type(retain) is not bool:
        raise ValueError('exact retention flag required')
    events=defaultdict(list); rows=[]; hasher=sha256(); old_mass=new_mass=0; maxd=0
    for s in range(1,256*n,2):
        if s%8 == 5:
            continue
        A=native.independent_charge(s,64*n)
        u,t,d=routes.route(s)
        maxd=max(maxd,d)
        norm=native.inherited.norm(A)
        row=[s,s%15,s%35,u,t,d,*A,norm]
        hasher.update((json.dumps(row,separators=(',',':'))+'\n').encode())
        if retain:rows.append(row)
        if s<4*n:
            old_mass+=norm
            check(A==native.independent_charge(s,n),'same old charge at64N')
        else:
            new_mass+=norm
            check(norm==native.outer_band_norm(s,n),'source-specific annular band norm')
        if norm:
            v=A if t%2 == 0 else native.negate(A)
            events[d].append((s,u,v,norm))
    check(new_mass==153*n,'complete annular source mass')
    old=defaultdict(lambda:(0,0));new=defaultdict(lambda:(0,0));ledger=[];arrived=0
    for d in range(maxd+1):
        for s,u,v,norm in events[d]:
            target=old if s<4*n else new
            target[u]=native.add(target[u],v)
            if s>=4*n:arrived+=norm
        us=old.keys()|new.keys()
        eo=sum(native.inherited.norm(v) for v in old.values())
        eh=sum(native.inherited.norm(v) for v in new.values())
        cross=sum(native.cross(old[u],new[u]) for u in us)
        combined=sum(native.inherited.norm(native.add(old[u],new[u])) for u in us)
        check(combined-eo==eh+cross,'actual finite-horizon polarization')
        ledger.append([d,eo,eh,cross,eh+cross,new_mass-arrived])
    out={'N':n,'source_count':sum(1 for s in range(1,256*n,2) if s%8!=5),
         'source_rows_sha256':hasher.hexdigest(),'old_norm':old_mass,'annulus_norm':new_mass,
         'ledger':ledger,'max_depth':maxd,'peak_increment':max(ledger,key=lambda r:r[4]),
         'peak_annulus_energy':max(ledger,key=lambda r:r[2]),'peak_cross':max(ledger,key=lambda r:r[3]),
         'unit_rate_passes_all_horizons':all(r[4]<=new_mass for r in ledger)}
    if retain:out['source_rows']=rows
    return out


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=json.loads((HERE/'terminals.json').read_text(encoding='utf-8'))
    if type(data) is not dict or type(data.get('rows')) is not list:
        raise ValueError('supplied terminal rows required')
    for r in data['rows']:
        if type(r) is not dict:raise ValueError('terminal row required')
        for key in ('source','child','cut','root_rank','root_cost'):integer(r[key])
        odd(r['source']);odd(r['child'])
        for key in ('first_exit_word','root_suffix'):
            if type(r[key]) is not list:raise ValueError('stored word list required')
            word(tuple(r[key]))
    supplied={r['source']:r for r in data['rows']}
    check(len(supplied)==76,'76 separately supplied child certificates')
    hist=Counter();gains=[];pending=[]
    for n in range(3,10002,2):
        w=first_exit(n);xs=replay(n,w);h=xs[-1]
        check(all(x>=n for x in xs[:-1]) and h<n,'literal first exit')
        lam,_,debt,*_=state(w[:-1]);baseline=lam>=1 and debt<=lam-1
        k=cone_cut(w)
        category='baseline' if baseline else 'gain' if k is not None else 'pending'
        hist[category]+=1
        check(state(w)[2]<n-h,'finite actual first-exit payment')
        if category=='gain':
            r=supplied[n]
            check(r['child']==h and tuple(r['first_exit_word'])==w and r['cut']==k,'immutable new source certificate')
            suffix=tuple(r['root_suffix']);whole=discharge(n,w,k,suffix)
            check(r['root_rank']==len(whole) and r['root_cost']==sum(whole),'grounded clock metadata')
            gains.append([n,h,len(w),k,len(suffix),sum(suffix)])
        elif category=='pending':pending.append(n)
    check(dict(hist)=={'baseline':4749,'gain':76,'pending':175},'declared finite controller comparison')
    check(set(supplied)=={r[0] for r in gains},'no unmatched supplied ROOT rows')
    # Independent suffix formula on a complete small word universe.
    language_counts=Counter()
    for length in range(1,7):
        for w in product(range(1,5),repeat=length):
            lam,pos,debt,A,Q,B=state(w)
            d=F(0);p=F(0)
            for i,a in enumerate(w):
                suffix=F(3**(length-i-1),2**sum(w[i+1:]))
                if a==1:p+=suffix
                if a>=3:d+=(1-F(4,2**a))*suffix
            check(d==debt and p==pos and lam==F(A,Q),'independent suffix debt')
            if terminal_language(w):
                check(debt<1,'all short grammar debt bounds');language_counts[length]+=1
    # Whole coarse ray has growing proper prefixes and uniform final payment.
    st=state(W63);check(st[3]==3**34 and st[4]==2**54,'63 word exact clocks')
    check(st[3]*63+st[5]==61*st[4],'63 affine endpoint')
    for i in range(1,len(W63)):
        l,_,_,_,_,b=state(W63[:i]);check(l>1 and b>0,'every proper prefix grows at every positive source')
    check(state(W63[:32])[2]<state(W63[:32])[0]-1,'new checkpoint is in cone')
    check(state(W63[:33])[2]>state(W63[:33])[0]-1,'old last-step cone fails on entire ray')
    lifts=[]
    for t in (*range(65),10**30):
        n=63+(1<<54)*t;r=coarse63(n);w=tuple(r['word'])
        check(r['endpoint']<=61+3**34*t and r['debt']==str(state(w)[2]),'coarse payment and clock retention')
        check(all(x>n for x in replay(n,w)[1:-1]),'literal first exit on sampled ray')
        lifts.append([t,r['last_valuation'],r['endpoint']])
    for n in (True,63.0,61,65):rejects(lambda n=n:coarse63(n))
    rejects(lambda:paid(63,W63,True));rejects(lambda:paid(63,W63,33))
    r=supplied[63];suffix=tuple(r['root_suffix'])
    rejects(lambda:discharge(63,W63,32,suffix+(2,)))
    rejects(lambda:discharge(65,W63,32,suffix))
    rejects(lambda:replay(63,list(W63)))
    # Exhaust all principal-unit targets through precision6 against a separate
    # enumeration of powers, then query the actual huge grounded family.
    for precision in range(1,7):
        modulus=3**precision;period=3**(precision-1)
        explicit={pow(4,j,modulus):j for j in range(period)}
        check(len(explicit)==period,'independent exact order of4')
        for target,j in explicit.items():
            check(power4_phase(target,precision)==(j,period),'all small principal-unit phases')
    j0,period=power4_phase(61,34)
    check(pow(4,j0,3**34)==61 and period==3**33 and 3<=j0<period,'actual exact ROOT phase')
    grounded_lifts=[]
    for index in (0,1,2,10**12):
        lift=root63_lift(index);j=lift['J'];residues=[]
        for modulus in (1,3,5,7,15,35,105,1<<54,3**38,105*(1<<60)):
            residue=root63_residue(index,modulus)
            check((st[3]*residue+st[5])%modulus==pow(2,54+2*j,modulus),'actual affine ROOT identity modulo full query')
            residues.append([modulus,residue])
        check(root63_residue(index,1<<54)==63,'grounded family stays in the new ray')
        check(root63_residue(index,105)%15==root63_residue(index,15) and
              root63_residue(index,105)%35==root63_residue(index,35),'compatible source15/35 registers')
        check(lift['word'][:-1]==W63[:-1] and lift['word'][-1]>=4,'strict ROOT terminal and marked head')
        grounded_lifts.append({'index':index,'J':j,'last_valuation':lift['word'][-1],'residues':residues})
    for bad in (True,1.0,-1):rejects(lambda bad=bad:root63_lift(bad))
    rejects(lambda:root63_residue(0,False));rejects(lambda:power4_phase(2,3))
    for index in range(60):
        joint,r15,r35=root63_registers(index)
        check((joint,r15,r35)==tuple(root63_residue(index,m) for m in (105,15,35)),
              'actual six-phase compatible15/35 register table')
    check(len(set(ROOT_REGISTERS))==6,'all six register phases distinct')
    check((pow(4,j0,105)-61)%(3**34)!=0,'denominator-blind reduction is not an integer quotient')
    rejects(lambda:root63_registers(True))
    # Source-indexed canonical-P first arrivals use the inherited fixed frontier.
    ranks={int(line.split()[0]):int(line.split()[1]) for line in native.FRONTIER.read_text().splitlines()
           if line.strip() and not line.startswith('#')}
    V=frozenset(ranks);routes=native.Routes(V)
    check(ranks[1]==0 and all(ranks[native.inherited.edge(s)[0]]<ranks[s] for s in V if s!=1),
          'supplied frontier has a strict ROOT rank')
    small=[annulus(n,routes,retain=n in (7,15)) for n in (*range(1,8),15)]
    check(all(r['unit_rate_passes_all_horizons'] for r in small[:6]),'N1..6 finite control')
    hostile=small[6]['ledger'][8]
    check(hostile==[8,19,1086,-2,1084,17] and small[6]['annulus_norm']==1071,'N7,d8 exact hostile')
    check(small[-1]['peak_increment']==[43,36,2625,8,2633,0],'N15 completed hostile')
    inherited=annulus(4096,routes)
    check(inherited['unit_rate_passes_all_horizons'] and inherited['annulus_norm']==626688,'N4096 positive renewal')
    check(inherited['peak_increment']==[149,9629,433053,-8039,425014,0],'N4096 exact maximum')
    check(inherited['peak_cross'][3]==0,'all actual N4096 cross terms nonpositive')
    # A small complete route certificate for retained hostile sources only.
    closure=set()
    for item in small:
        if 'source_rows' not in item:continue
        for row in item['source_rows']:
            s=row[0]
            while s not in V:
                closure.add(s);s=routes.edges[s][0]
    dag=[]
    for s in sorted(closure):
        nxt,dt=routes.edges[s];u,t,d=routes.memo[s];uu,tt,dd=routes.memo[nxt]
        check((nxt,dt)==native.inherited.edge(s) and (u,t,d)==(uu,tt+dt,dd+1),'small actual first-arrival DAG')
        dag.append([s,nxt,dt,u,t,d])
    payload={'status':'PROVED scoped grammar/ray; FINITE-EXACT grounding/renewal; universal payment/entry OPEN',
             'checks':CHECKS,'finite_sources':[3,10001,2],'controller_counts':dict(hist),'new_grounded_sources':gains,
             'pending_sources':pending,'terminal_data_sha256':sha256((HERE/'terminals.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
             'terminal_language_small_counts':dict(language_counts),'coarse63':{'modulus':2**54,'residue':63,
             'word_at_63':W63,'carrier':[st[3],st[4],st[5]],'debt':str(st[2]),'sample_lifts':lifts},
             'grounded_ray':{'phase_J':j0,'period_J':period,'source_formula':root63_lift(0)['source_formula'],
                             'sample_modular_readers':grounded_lifts,
                             'register_period':6,'register_table':[[i,*root63_registers(i)] for i in range(6)]},
             'annulus_small':small,'annulus_inherited':inherited,'small_route_dag':dag,
             'source_row_schema':['s','mod15','mod35','terminal','clock','depth','a','b','norm'],
             'ledger_schema':['depth','old_energy','annulus_energy','cross','increment','moving_annulus_norm'],
             'frontier':sorted(V),'frontier_ranks':sorted(ranks.items()),'frontier_sha256':native.digest(sorted(V)),
             'native_dependency_sha256':sha256((REPO/'04-computation/experiments/native_cancellation_105_20261018g.py').read_bytes().replace(b'\r\n',b'\n')).hexdigest()}
    text=json.dumps(payload,sort_keys=True,separators=(',',':'))+'\n';out=HERE/'receipt.json'
    if args.verify:
        if out.read_text(encoding='utf-8').replace('\r\n','\n')!=text:raise RuntimeError('saved receipt differs')
    else:out.write_text(text,encoding='utf-8',newline='\n')
    print(json.dumps({'checks':CHECKS,'counts':dict(hist),'annulus_hostile':hostile,
                     'inherited_peak':inherited['peak_increment'],'LF_SHA256':sha256(text.encode()).hexdigest()},sort_keys=True))


if __name__=='__main__':main()
