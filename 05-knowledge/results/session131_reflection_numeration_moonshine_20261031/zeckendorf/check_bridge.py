#!/usr/bin/env python3
"""Portable exact checks; no imports of producer code or third-party libraries."""
from collections import deque
from itertools import product, combinations
from math import isqrt
from random import Random
from pathlib import Path
import hashlib
import json
import sys

CHECKS = 0

def req(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(message)


def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def weights(kind, size):
    out = [1, 3] if kind == 'cg' else [1, 2]
    while len(out) < size:
        out.append(3*out[-1]-out[-2] if kind == 'cg' else out[-1]+out[-2])
    return out[:size]


def value(digits, kind='cg'):
    return sum(a*w for a, w in zip(digits, weights(kind, len(digits))))


def encode(n, kind='cg'):
    req(isinstance(n, int) and n >= 0, 'encode domain')
    ws = weights(kind, 2)
    while ws[-1] <= n:
        ws.append(3*ws[-1]-ws[-2] if kind == 'cg' else ws[-1]+ws[-2])
    result = []
    for w in reversed(ws):
        q, n = divmod(n, w)
        result.append(q)
    return trim(reversed(result))


def legal(digits, kind='cg'):
    state = 0
    for a in reversed(digits):
        if kind == 'cg':
            if a not in (0, 1, 2) or state and a == 2:
                return False
            state = int(a == 2 or state and a == 1)
        else:
            if a not in (0, 1) or state and a == 1:
                return False
            state = a
    return True


def beatty(n):
    req(n >= 0, 'positive Beatty domain')
    return 0 if n == 0 else (3*n+isqrt(5*n*n))//2+1


def signed_shift(n):
    return beatty(n) if n >= 0 else -beatty(-n)


def signed_z_shift(n):
    return (1 if n >= 0 else -1)*value([0]+encode(abs(n),'z'),'z')


def raw_shift(digits):
    return value([0]+list(digits))


def normalize(coefficients, high_first=False):
    a = trim(coefficients)
    req(all(isinstance(d, int) and d >= 0 for d in a), 'normalization domain')
    before, shifted_before = value(a), raw_shift(a)
    bound, steps, bottom = sum(a), 0, 0
    while True:
        big = [i for i, d in enumerate(a) if d >= 3]
        if big:
            i = big[-1] if high_first else big[0]
            if i+1 == len(a):
                a.append(0)
            a[i] -= 3
            a[i+1] += 1
            if i:
                a[i-1] += 1
            else:
                bottom += 1
        else:
            forbidden = []
            for lo, d in enumerate(a):
                if d != 2:
                    continue
                hi = lo+1
                while hi < len(a) and a[hi] == 1:
                    hi += 1
                if hi < len(a) and a[hi] == 2:
                    forbidden.append((lo, hi))
            if not forbidden:
                break
            lo, hi = forbidden[-1] if high_first else forbidden[0]
            if hi+1 == len(a):
                a.append(0)
            a[lo] -= 2
            a[hi] -= 2
            for j in range(lo+1, hi):
                a[j] -= 1
            a[hi+1] += 1
            if lo:
                a[lo-1] += 1
            else:
                bottom += 1
        steps += 1
        req(steps <= bound, 'decreasing sum termination')
    a = trim(a)
    req(value(a) == before, 'normalization value')
    req(legal(a), 'normalization legality')
    req(shifted_before-raw_shift(a) == bottom, 'exact zero-boundary carry count')
    return a, bottom, steps


def z_to_signed_even(z):
    req(legal(z, 'z'), 'Zeckendorf input')
    e = list(z[0::2]); o = list(z[1::2])
    length = max(len(e), len(o)+1)
    e += [0]*(length-len(e)); o += [0]*(length-len(o))
    return trim(e[i]+(o[i-1] if i else 0)-o[i] for i in range(length))


def signed_even_to_z(c):
    previous = 0
    out = []
    for digit in c:
        if previous == 0:
            possibilities = {0:(0,0), 1:(1,0), -1:(0,1)}
        else:
            possibilities = {0:(0,1), 1:(0,0)}
        if digit not in possibilities:
            raise ValueError('illegal signed even digit')
        e, o = possibilities[digit]
        out.extend((e,o))
        previous = o
    if previous:
        raise ValueError('unclosed signed even carry')
    return trim(out)


SINK = (-1,-1,-1,-1)
START = (0,0,0,0)

def transition(state, ab, kind):
    if state == SINK:
        return SINK
    l, r, n, p = state
    a, b = ab
    t = (a+b) % 2
    if kind == 'cg':
        if l and a == 2 or r and b == 2:
            return SINK
        return (int(a==2 or l and a==1), int(b==2 or r and b==1), (p+t)%2, (n+p+t)%2)
    if l and a or r and b:
        return SINK
    return (a,b,(p+t)%2,(n+p)%2)


def accepts(state):
    return state != SINK and state[2] == 0


def run_dfa(pairs, kind, state=START):
    for ab in pairs:
        state = transition(state, ab, kind)
    return state


def make_dfa_certificate(kind):
    alphabet = list(product(range(3 if kind == 'cg' else 2), repeat=2))
    reach = {START:()}
    queue = deque([START])
    while queue:
        state = queue.popleft()
        for ab in alphabet:
            nxt = transition(state, ab, kind)
            if nxt not in reach:
                reach[nxt] = reach[state]+(ab,)
                queue.append(nxt)
    states = sorted(reach)
    req(len(states) == 17, kind+' reachable state count')
    for state, word in reach.items():
        req(run_dfa(word,kind) == state, 'reachable-state witness')
    distinguishes = []
    for left, right in combinations(states,2):
        queue = deque([(left,right,())])
        visited = {(left,right)}
        suffix = None
        while queue:
            a,b,w = queue.popleft()
            if accepts(a) != accepts(b):
                suffix = w
                break
            for symbol in alphabet:
                aa,bb = transition(a,symbol,kind),transition(b,symbol,kind)
                if (aa,bb) not in visited:
                    visited.add((aa,bb))
                    queue.append((aa,bb,w+(symbol,)))
        req(suffix is not None, 'DFA minimality witness exists')
        req(accepts(run_dfa(suffix,kind,left)) != accepts(run_dfa(suffix,kind,right)), 'distinguishing suffix')
        distinguishes.append({'left':left,'right':right,'suffix':suffix})
    return {'kind':kind,'alphabet':alphabet,'start':START,'sink':SINK,
            'states':[{'state':s,'accept':accepts(s),'reach':reach[s],
                       'next':[transition(s,a,kind) for a in alphabet]} for s in states],
            'distinguish':distinguishes}


BLUE = [(0,5),(3,8),(4,6),(5,6),(5,8),(6,9),(7,9),(8,9)]
BLACK = [(1,4),(2,4),(2,5),(2,6),(2,7),(2,8),(4,5),(4,8),(5,7),(6,6),(6,8),(8,8)]
Q = BLUE+BLACK
G = BLUE+[e for e in BLACK for _ in (0,1)]


def boundary(edges, z):
    out = [0]*10
    for (u,v),a in zip(edges,z):
        out[u] -= a; out[v] += a
    return out


def fundamental_cycles(edges):
    # Independent edge-scan tree construction with signed unique tree paths.
    parent = list(range(10))
    adjacency = [[] for _ in range(10)]
    def root(v):
        while parent[v] != v:
            v = parent[v]
        return v
    chords = []
    for i,(u,v) in enumerate(edges):
        a,b = root(u),root(v)
        if a==b:
            chords.append(i)
        else:
            parent[a] = b
            adjacency[u].append((v,i,1));adjacency[v].append((u,i,-1))
    basis = []
    for i in chords:
        u,v = edges[i]
        queue = deque([(v,[])])
        seen = {v}
        path = None
        while queue:
            w,p = queue.popleft()
            if w==u:
                path=p;break
            for nxt,j,sign in adjacency[w]:
                if nxt not in seen:
                    seen.add(nxt);queue.append((nxt,p+[(j,sign)]))
        req(path is not None,'tree path')
        z = [0]*len(edges);z[i]=1
        for j,sign in path:z[j]=sign
        req(boundary(edges,z)==[0]*10,'fundamental cycle')
        basis.append(z)
    return basis


def glue_encode(z):
    b=z[:8];u=[];v=[];epsilon=[]
    for x,y in zip(z[8::2],z[9::2]):
        s,d=x+y,x-y;e=s%2
        u.append((s-e)//2);v.append((d-e)//2);epsilon.append(e)
    return b,u,v,epsilon


def glue_decode(code):
    b,u,v,e=code
    return list(b)+[c for a,d,t in zip(u,v,e) for c in (a+d+t,a-d)]


def glue_add(left,right):
    b,u,v,e=left;bb,uu,vv,ee=right
    return ([a+c for a,c in zip(b,bb)],
            [a+c+t*w for a,c,t,w in zip(u,uu,e,ee)],
            [a+c+t*w for a,c,t,w in zip(v,vv,e,ee)],
            [t^w for t,w in zip(e,ee)])


def glue_reflect(code):
    b,u,v,e=code
    return b,u,[-d-t for d,t in zip(v,e)],e


def check_lattice(rng):
    qbasis=fundamental_cycles(Q);gbasis=fundamental_cycles(G)
    req((len(qbasis),len(gbasis))==(11,23),'cycle ranks')
    image={0}
    for z in qbasis:
        bit=sum((a%2)<<i for i,a in enumerate(z[8:]))
        image |= {x^bit for x in list(image)}
    target={x for x in range(4096) if x&1==0 and sum((x>>i)&1 for i in range(1,6))%2==0}
    req(image==target and len(image)==1024,'realized gluing code')
    test_cycles=gbasis+[[rng.randrange(-25,26) for _ in range(23)] for _ in range(1000)]
    # The first 23 inputs are already edge cycles; later inputs are basis coordinates.
    prior=[0]*32
    for ix,raw in enumerate(test_cycles):
        z=raw if ix<23 else [sum(c*b[j] for c,b in zip(raw,gbasis)) for j in range(32)]
        req(boundary(G,z)==[0]*10,'cycle sample')
        code=glue_encode(z)
        req(glue_decode(code)==z,'gluing inverse')
        b,u,v,e=code
        req(boundary(Q,b+[2*a+t for a,t in zip(u,e)])==[0]*10,'gluing boundary')
        req(sum(t<<i for i,t in enumerate(e)) in image,'gluing residue image')
        reflected=z[:8]+[c for x,y in zip(z[8::2],z[9::2]) for c in (y,x)]
        req(glue_decode(glue_reflect(code))==reflected,'reflection formula')
        req(glue_decode(glue_add(code,glue_encode(prior)))==[a+b for a,b in zip(z,prior)],'carry addition')
        prior=z
    # A literal loop is a genuine cycle, so the arithmetic hostile is graph-realized.
    z=[0]*32;z[8+2*9]=2;z[9+2*9]=1
    req(boundary(G,z)==[0]*10,'positive loop hostile is cycle')
    req((beatty(3),beatty(1))==(8,3),'independent-shift hostile')
    req((beatty(2)+beatty(1),beatty(2)-beatty(1))==(9,3),'corrected-shift hostile')
    return {'cycle_samples':len(test_cycles),'epsilon_words':len(image),
            'epsilon_constraints':['epsilon_0=0','epsilon_1+...+epsilon_5=0 mod 2'],
            'image_sha256':hashlib.sha256(json.dumps(sorted(image)).encode()).hexdigest()}


def check_numerals(rng):
    for n in range(50001):
        cg,z=encode(n),encode(n,'z')
        req(legal(cg) and legal(z,'z'),'canonical legality')
        req(value(cg)==value(z,'z')==n,'canonical reconstruction')
        req(raw_shift(cg)==beatty(n),'shift ceiling identity')
        c=z_to_signed_even(z)
        req(value(c)==n and signed_even_to_z(c)==z,'interwoven reversible map')
    for _ in range(512):
        n=rng.getrandbits(512)
        for kind in ('cg','z'):
            a=encode(n,kind)
            req(legal(a,kind) and value(a,kind)==n,'large canonical reconstruction')
        z=encode(n,'z');c=z_to_signed_even(z)
        req(signed_even_to_z(c)==z and value(c)==n,'large interwoven inverse')
        req(raw_shift(encode(n))==beatty(n),'large exact shift')
    # Enumerate all signed-even strings through length 8, including invalid ones.
    valid_signed=0
    for length in range(9):
        for c in product((-1,0,1),repeat=length):
            try:z=signed_even_to_z(c)
            except ValueError:continue
            req(legal(z,'z'),'signed-language nonconsecutive output')
            req(z_to_signed_even(z)==trim(c),'signed-language reverse')
            req(value(c)==value(z,'z'),'signed-language value')
            valid_signed+=1
    return {'integers_through':50000,'random_512_bit_integers':512,
            'signed_even_words_lengths':'0..8','accepted_signed_even_words':valid_signed}


def check_dfa(kind,certificate,rng):
    alphabet=certificate['alphabet']
    length_limit=4 if kind=='cg' else 7
    total=0
    for length in range(length_limit+1):
        for word in product(alphabet,repeat=length):
            a=[x[0] for x in reversed(word)];b=[x[1] for x in reversed(word)]
            expected=legal(a,kind) and legal(b,kind) and (value(a,kind)-value(b,kind))%2==0
            req(accepts(run_dfa(word,kind))==expected,'DFA exhaustive evaluation')
            total+=1
    for _ in range(1024):
        a=encode(rng.getrandbits(256),kind);b=encode(rng.getrandbits(256),kind)
        length=max(len(a),len(b));a +=[0]*(length-len(a));b +=[0]*(length-len(b))
        word=list(zip(reversed(a),reversed(b)))
        req(accepts(run_dfa(word,kind))==((value(a,kind)-value(b,kind))%2==0),'large DFA parity')
    return {'states':17,'distinguished_pairs':136,'exhaustive_padded_pair_words':total,
            'exhaustive_lengths':f'0..{length_limit}','random_256_bit_pairs':1024}


def check_carries(rng):
    total=0
    for a in range(129):
        da=encode(a)
        for b in range(129):
            db=encode(b);length=max(len(da),len(db))
            raw=[(da[i] if i<len(da) else 0)+(db[i] if i<len(db) else 0) for i in range(length)]
            expected=encode(a+b);defect=beatty(a)+beatty(b)-beatty(a+b)
            req(defect in (0,1),'two-summand carry range')
            for high in (False,True):
                result,bottom,steps=normalize(raw,high)
                req(result==expected and bottom==defect,'strategy-independent boundary carry')
            total+=1
    for _ in range(256):
        summands=[rng.getrandbits(256) for _ in range(rng.randrange(2,9))]
        words=[encode(n) for n in summands];length=max(map(len,words))
        raw=[sum(w[i] for w in words if i<len(w)) for i in range(length)]
        defect=sum(beatty(n) for n in summands)-beatty(sum(summands))
        req(0<=defect<len(summands),'many-summand carry bound')
        for high in (False,True):
            result,bottom,steps=normalize(raw,high)
            req(result==encode(sum(summands)) and bottom==defect,'large normalization carry')
    for a in range(-100,101):
        for b in range(-100,101):
            defect=signed_shift(a)+signed_shift(b)-signed_shift(a+b)
            req(defect in (-1,0,1),'signed carry range')
            s,d=a+b,a-b
            ss,dd=signed_shift(a)+signed_shift(b),signed_shift(a)-signed_shift(b)
            req((ss-dd)%2==0,'repaired eigen gluing')
            req(abs(ss-signed_shift(s))<=1 and abs(dd-signed_shift(d))<=1,'bounded eigen repair')
            req((ss+dd)//2==signed_shift(a) and (ss-dd)//2==signed_shift(b),'repair inverse')
    for a in range(-100,101):
        one=signed_shift(a);two=signed_shift(one);three=signed_shift(two)
        req(two==3*one-a and three==8*one-3*a,'iterated shift identity')
        req((three-a)%2==0,'three shifts preserve scalar parity')
    for a in range(-40,41):
        for b in range(-40,41):
            if (a-b)%2:continue
            aa,bb=a,b
            second=(signed_shift(a)-signed_shift(b))%2
            for k in range(10):
                req(((aa-bb)%2==0)==(not second or k%3==0),'exact gluing return period')
                aa,bb=signed_shift(aa),signed_shift(bb)
    for a in range(-100,101):
        one=signed_z_shift(a);two=signed_z_shift(one);three=signed_z_shift(two)
        req(two==one+a and three==2*one+a,'Zeckendorf iterated shift identity')
        req((three-a)%2==0,'Zeckendorf three shifts preserve scalar parity')
    for a in range(-20,21):
        for b in range(-20,21):
            if (a-b)%2:continue
            aa,bb=a,b
            second=(signed_z_shift(a)-signed_z_shift(b))%2
            for k in range(10):
                req(((aa-bb)%2==0)==(not second or k%3==0),'Zeckendorf exact gluing return period')
                aa,bb=signed_z_shift(aa),signed_z_shift(bb)
    # Explicit non-equivalences are load-bearing hostile controls.
    req(legal([1],'z') and legal([0,1],'z') and (1-2)%2!=0,'nonconsecutive does not imply gluing')
    req(not legal([1,1],'z') and value([1,1],'z')==3,'naive interleaving forbidden')
    req((beatty(3)-beatty(1))%2!=0,'independent shift loses gluing')
    req(normalize([3])[1]==1,'sharp bottom carry')
    return {'all_nonnegative_pairs':'0..128 x 0..128','pairs':total,
            'strategies':['lowest available rewrite','highest available rewrite'],
            'random_multi_sums_256_bit':256,'signed_pair_box':'[-100,100]^2',
            'shift_gluing_period_box':'[-40,40]^2, equal parity only; iterations 0..9',
            'zeckendorf_shift_period_box':'[-20,20]^2, equal parity only; iterations 0..9'}


def main():
    if '--test-failing-gate' in sys.argv:
        req(False,'intentional hostile failing gate')
    rng=Random(13120261031)
    certificates={kind:make_dfa_certificate(kind) for kind in ('cg','z')}
    numeral_report=check_numerals(rng)
    dfa_reports={kind:check_dfa(kind,certificates[kind],rng) for kind in certificates}
    carry_report=check_carries(rng)
    lattice_report=check_lattice(rng)
    report={'status':'PASS','explicit_checks':CHECKS,'numerals':numeral_report,'dfa':dfa_reports,
            'carries':carry_report,'lattice':lattice_report,'minimality_certificates':certificates,
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
