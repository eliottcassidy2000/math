"""Exact small-module contraction and primitive protected-port certificates."""
from itertools import product
from math import gcd, isqrt, prod
from pathlib import Path
import json

INHERITED_REF = "087740676b20b54adf90171584b8833855805b3c"
A = (2,28,25,81,97,79,7)
B = (2,80,3,7,173,143,189,15)
Q4 = (2,4,9,3)
P = 8983331801295497479938493370497


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def tt(n):
    return tuple(((1 << n)-1) ^ ((1 << (i+1))-1) for i in range(n))


def ordinal(*blocks):
    n, offset, rows = sum(map(len, blocks)), 0, []
    for block in blocks:
        later = ((1 << n)-1) ^ ((1 << (offset+len(block)))-1)
        rows.extend((row << offset) | later for row in block)
        offset += len(block)
    return tuple(rows)


D = ordinal(A, tt(1))


def cyclic(blocks):
    need(len(blocks) == 3 and all(blocks), 'three nonempty cyclic blocks')
    starts = (0,len(blocks[0]),len(blocks[0])+len(blocks[1]))
    rows = []
    for i, block in enumerate(blocks):
        nxt = (i+1) % 3
        outgoing = ((1 << len(blocks[nxt]))-1) << starts[nxt]
        rows.extend((row << starts[i]) | outgoing for row in block)
    return tuple(rows)


def delete(rows, root):
    retained = [i for i in range(len(rows)) if i != root]
    return tuple(sum(((rows[i] >> j) & 1) << k for k,j in enumerate(retained))
                 for i in retained)


def inventory(rows):
    n = len(rows)
    ends = [[0]*n for _ in range(1 << n)]
    H = [1]+[0]*((1 << n)-1)
    for mask in range(1,1 << n):
        for v in range(n):
            if mask >> v & 1:
                previous = mask ^ (1 << v)
                ends[mask][v] = 1 if not previous else sum(
                    ends[previous][u] for u in range(n) if rows[u] >> v & 1)
        H[mask] = sum(ends[mask])
    return H, ends


def covers(rows):
    H, _ = inventory(rows)
    n = len(rows)
    previous = [1]+[0]*((1 << n)-1)
    result = [0]
    for _ in range(1,n+1):
        current = [0]*(1 << n)
        for mask in range(1,1 << n):
            first = mask
            while first:
                current[mask] += H[first]*previous[mask ^ first]
                first = (first-1) & mask
        result.append(current[-1])
        previous = current
    return result


def cyclic_slots(profiles, slots):
    """All cyclic words in each labelled path slot, including empty words."""
    caps = tuple(len(p)-1 for p in profiles)
    words = {(0,0,0): 1}
    for length in range(1,sum(caps)+1):
        for start in range(3):
            visits = [length//3]*3
            for offset in range(length % 3):
                visits[(start+offset) % 3] += 1
            if all(visits[i] <= caps[i] for i in range(3)):
                key = tuple(visits)
                words[key] = words.get(key,0)+1
    coefficients = {(0,0,0): 1}
    for _ in range(slots):
        following = {}
        for old,a in coefficients.items():
            for word,b in words.items():
                new = tuple(old[i]+word[i] for i in range(3))
                if all(new[i] <= caps[i] for i in range(3)):
                    following[new] = following.get(new,0)+a*b
        coefficients = following
    return sum(c*prod(profiles[i][visits[i]] for i in range(3))
               for visits,c in coefficients.items())


def state(rows, root):
    H, E = inventory(rows)
    full = (1 << len(rows))-1
    opposite = tuple(full ^ (1 << i) ^ row for i,row in enumerate(rows))
    _, opposite_E = inventory(opposite)
    return (H[-1],opposite_E[full][root],H[full ^ (1 << root)],
            sum(E[mask][root]*H[full ^ mask] for mask in range(1,full+1)))


def birth(rows, root, command):
    new = len(rows)
    output = list(rows)
    output[root] |= 1 << new
    output.append(((1 << new)-1) ^ (1 << root))
    return tuple(output),root if command == 'R' else new


def trial_prime(n):
    return n >= 2 and all(n % d for d in range(2,isqrt(n)+1))


def character_certificate(aa,bb,cc,dd):
    small_primes = (2,3,17,29,31,83,733,17497757)
    F = 2**7*prod(small_primes[1:])
    need(all(trial_prime(q) for q in small_primes), 'small prime factors by trial division')
    need((P-1) % F == 0 and F*F > P, 'large known order divisor')
    witness = 5
    need(pow(witness,P-1,P) == 1, 'Fermat residue')
    residues = {q:gcd(pow(witness,(P-1)//q,P)-1,P) for q in small_primes}
    need(all(v == 1 for v in residues.values()), 'every prime divisor has order multiple F')
    need(abs(aa*dd-bb*cc) == 4*9*7*P, 'complete resultant factorization')
    bad = (-aa*pow(bb,-1,P)) % P
    need((cc+dd*bad) % P == 0, 'same bad projective ratio')
    need(pow(2,(P-1)//2,P) == 1, 'two has positive quadratic character')
    need(pow(bad,(P-1)//2,P) == P-1, 'bad ratio has negative character')
    return dict(prime=P, known_factor=F, unfactored_quotient=(P-1)//F,
                small_primes=small_primes, witness=witness,
                order_gcds=residues, bad_ratio=bad,
                character_exponent=(P-1)//2, character_two=1, character_bad=-1)


def literal_family(n,k,broken):
    m = 4*n+3+6*k
    host = ordinal(*([Q4]*(2*n)),tt(m))
    inside = cyclic((B,D if broken else B,B))
    graph = cyclic((ordinal(host,inside),tt(1),tt(1)))
    root = len(host)  # vertex0 in the first, unchanged B
    return graph,root


def main():
    cb, cd, cb0 = covers(B), covers(D), covers(delete(B,0))
    need(cb == cd, 'inherited complete neutral profile in literal names')
    aa = cyclic_slots((cb,cb,cb),1)
    bb = cyclic_slots((cb,cb,cb),2)
    cc = cyclic_slots((cb0,cb,cb),1)
    dd = cyclic_slots((cb0,cb,cb),2)
    need((aa,bb,cc,dd) == (127749360032069607,5829968289840209568,
                           10335194072037451,453936161252305532), 'core contraction constants')
    need(gcd(aa,cc) == 1, 'unwrapped protected primitive port')
    certificate = character_certificate(aa,bb,cc,dd)
    small_controls = []
    for blocks in ((tt(1),tt(1),tt(1)),(tt(2),(2,4,1),tt(1)),(B,tt(1),tt(1))):
        profiles = tuple(covers(block) for block in blocks)
        graph = cyclic(blocks)
        H,_ = inventory(graph)
        full = len(H)-1
        actual = (H[-1],sum(H[x]*H[full ^ x] for x in range(full+1)))
        predicted = tuple(cyclic_slots(profiles,j) for j in (1,2))
        need(actual == predicted, 'cyclic word contraction versus full subset DP')
        small_controls.append(dict(block_orders=list(map(len,blocks)), values=actual))
    # Literal protected external-port controls of the module congruence.
    word_controls = []
    for word in ('','R','N','RR','NR'):
        observations = []
        for block in (B,D):
            graph = cyclic((block,tt(1),tt(1)))
            root = len(block)
            for command in word:
                graph,root = birth(graph,root,command)
            observations.append(state(graph,root))
        need(observations[0] == observations[1], 'protected-root program congruence')
        word_controls.append(dict(word=word,state=observations[0]))
    H,s,q,f = word_controls[0]['state']
    observed_totals = {item['word']: item['state'][0] for item in word_controls}
    need(observed_totals['R'] == H+f-s and observed_totals['RR'] == H+3*f-s-2*q
         and observed_totals['NR'] == H+2*f, 'universal root conventions')
    need(gcd(gcd(H,observed_totals['R']),gcd(observed_totals['RR'],observed_totals['NR']))
         == gcd(gcd(H,s),gcd(q,f)), 'four-observer content identity')
    guard_controls = []
    for n in range(5):
        for m in range(12):
            h,p = 25**n,2**(10*n+m)
            H,q = aa*h+bb*p,cc*h+dd*p
            expected = m % 2 == 1 and m % 3 != (n+1) % 3
            need((gcd(H,q) == 1) == expected, 'all-m guard hostile and positive controls')
            guard_controls.append(dict(n=n,m=m,gcd=gcd(H,q),primitive=expected))
    examples = []
    for n,k in ((0,0),(0,1),(1,0)):
        graph,root = literal_family(n,k,False)
        changed,_ = literal_family(n,k,True)
        h,p = 25**n,2**(14*n+3+6*k)
        H,q = aa*h+bb*p,cc*h+dd*p
        need(gcd(H,q) == 1, 'selected infinite family')
        protected = set(range(len(graph))) - set(range(root+8,root+16))
        need(all(((graph[i] >> j) & 1) == ((changed[i] >> j) & 1)
                 for i in protected for j in range(len(graph)) if i != j),
             'all incidences outside changed second leaf retained')
        examples.append(dict(n=n,k=k,order=len(graph),root=root,
                             H=H,q=q,rows=graph,broken_rows=changed,
                             counting_method='proved small-profile contraction, not large subset DP'))
    # Induced host growth n -> n+1 retains U and every old named vertex.
    before,old_root = literal_family(0,0,True)
    after,new_root = literal_family(1,0,True)
    mapping = [i+8 if i < 3 else i+12 for i in range(len(before))]
    need(mapping[old_root] == new_root, 'retained root under changing host')
    need(all(((before[i] >> j)&1) == ((after[mapping[i]] >> mapping[j])&1)
             for i in range(len(before)) for j in range(len(before))), 'literal induced host growth')
    reservoir_controls = []
    for ell in range(12):
        h,p = 95**ell,8*1276**ell
        content = gcd(aa*h+bb*p,cc*h+dd*p)
        need((content == 1) == (ell % 6 in (0,4)), 'right reservoir primitive checkpoint')
        reservoir_controls.append(dict(ell=ell,gcd=content))
    output = dict(inherited_ref=INHERITED_REF, constants=dict(A=aa,B=bb,C=cc,D=dd),
                  B_profile=cb,B_delete0_profile=cb0,prime_character_certificate=certificate,
                  contraction_controls=small_controls,program_controls=word_controls,
                  guard_controls=guard_controls,examples=examples,
                  host_growth_embedding=mapping,
                  right_reservoir_controls=reservoir_controls,
                  scope='small module endpoint/partition DP; universal cyclic itinerary contraction; finite order/character certificate; no order29+ path census')
    Path(__file__).with_name('receipt.json').write_text(json.dumps(output,indent=2)+'\n')
    print('PASS: protected primitive family, exact character certificate and literal host growth')


if __name__ == '__main__':
    main()
