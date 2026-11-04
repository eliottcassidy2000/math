#!/usr/bin/env python3
"""Independent scalar audit: rebuild every atom by permutations, not source DP.

No producer module is imported. Full-signature sets replace producer literal
word lists. A deliberately larger normalized h-envelope is used. Read-only
against the producer receipt; prints a deterministic audit record.
"""
from collections import Counter
from functools import lru_cache
from itertools import permutations
from pathlib import Path
import hashlib
import json

GATES = 0


def require(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)


def rows_at(n, mask):
    result = [0] * n
    k = 0
    for i in range(n):
        for j in range(i+1,n):
            a,b = (i,j) if (mask >> k) & 1 else (j,i)
            result[a] |= 1 << b
            k += 1
    return tuple(result)


def strong(rows):
    r = [v | (1 << i) for i,v in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if (r[i] >> k) & 1:
                r[i] |= r[k]
    return all(v == (1 << len(rows))-1 for v in r)


def all_signatures(rows):
    n = len(rows)
    paths = [[] for _ in range(1 << n)]
    paths[0] = [()]
    for size in range(1,n+1):
        for path in permutations(range(n),size):
            if all((rows[a] >> b) & 1 for a,b in zip(path,path[1:])):
                paths[sum(1 << v for v in path)].append(path)
    full = (1 << n)-1
    h,p = len(paths[full]),0
    fv,gv,ev,sv = ([0]*n for _ in range(4))
    joint = [[0]*n for _ in range(n)]
    for path in paths[full]:
        ev[path[-1]] += 1
        sv[path[0]] += 1
    for mask,left in enumerate(paths):
        for a in left:
            for b in paths[full ^ mask]:
                p += 1
                if a: fv[a[-1]] += 1
                if b: gv[b[0]] += 1
                if a and b: joint[a[-1]][b[0]] += 1
    for mask in range(1 << n):
        W = [i for i in range(n) if (mask >> i) & 1]
        yield (h,p,sum(fv[i] for i in W),sum(gv[i] for i in W),
               sum(joint[i][j] for i in W for j in W),
               sum(ev[i] for i in W),sum(sv[i] for i in W))


def multiply(a,b):
    h,p,f,g,j,e,s = a
    H,P,F,G,J,E,S = b
    return h*H,p*P,f*H+p*F,g*P+h*G,j*H+h*J+e*S+g*F,h*E,s*H


def odd(n):
    while n % 2 == 0:
        n //= 2
    return n


def ordered_union(a,b):
    shift = len(a)
    new = ((1 << len(b))-1) << shift
    return tuple(row | new for row in a) + tuple(row << shift for row in b)


def literal_signature(rows,ports):
    """End-point removal recurrence, using converse for start counts."""
    n = len(rows)
    incoming = tuple(sum(1 << j for j in range(n) if (rows[j] >> i) & 1)
                     for i in range(n))
    def end_table(preds):
        table = [[0]*n for _ in range(1 << n)]
        for v in range(n): table[1 << v][v] = 1
        for mask in range(1,1 << n):
            if mask & (mask-1) == 0: continue
            vertices = mask
            while vertices:
                bit = vertices & -vertices; vertices -= bit
                v = bit.bit_length()-1
                prev = mask ^ bit
                choices = prev & preds[v]
                value = 0
                while choices:
                    old = choices & -choices; choices -= old
                    value += table[prev][old.bit_length()-1]
                table[mask][v] = value
        return table
    ends = end_table(incoming)
    reverse_ends = end_table(rows)
    H = [sum(v) for v in ends]; H[0] = 1
    E = [sum(v[i] for i in range(n) if (ports >> i) & 1) for v in ends]
    S = [sum(v[i] for i in range(n) if (ports >> i) & 1) for v in reverse_ends]
    full = (1 << n)-1
    return H[full],sum(H[m]*H[full^m] for m in range(full+1)),\
        sum(E[m]*H[full^m] for m in range(full+1)),\
        sum(H[m]*S[full^m] for m in range(full+1)),\
        sum(E[m]*S[full^m] for m in range(full+1)),E[full],S[full]


def attach(rows,ports,t):
    prefix = tuple(sum(1 << j for j in range(i+1,t)) for i in range(t))
    out = list(ordered_union(prefix,rows))
    root = len(out)
    root_row = (1 << t)-1
    for i in range(len(rows)):
        if (ports >> i) & 1:
            out[t+i] |= 1 << root
        else:
            root_row |= 1 << (t+i)
    return tuple(out+[root_row])


def main():
    receipt_path = Path(__file__).with_name('shifted_returns.json')
    raw = receipt_path.read_bytes()
    receipt = json.loads(raw)
    bank = set()
    shape_counts = Counter()
    for n in range(1,6):
        for mask in range(1 << (n*(n-1)//2)):
            rows = rows_at(n,mask)
            if not strong(rows):
                continue
            shape_counts[n] += 1
            for sig in all_signatures(rows):
                h,p,f,g,j,e,s = sig
                require(h % 2 == 1 and p % 2 == 0 and p >= 2*h, 'atom parity')
                require(j % 2 == 0 and (f == 0) == (g == 0), 'actual marked atom')
                if f:
                    require(f >= h and g >= h, 'both marked cut injections')
                bank.add(sig)
    atoms = tuple(sorted(bank))
    require(len(atoms) == 170, 'full independently rebuilt bank')
    require(dict(shape_counts) == {1:1,3:2,4:24,5:544}, 'strong labeled counts')
    bank_receipt = json.loads(Path(__file__).with_name('component_inverse_five.json').read_text())
    require({tuple(a['sigma']) for a in bank_receipt['atoms']} == bank,
            'producer and permutation-derived bank sets agree exactly')
    for a in bank_receipt['atoms']:
        require(strong(tuple(a['rows'])) and
                literal_signature(tuple(a['rows']),a['ports']) == tuple(a['sigma']),
                'each producer atom has its stated actual representative')
    marked = tuple(a for a in atoms if a[2])
    zero = tuple(a for a in atoms if not a[2])
    M = max(a[0] for a in atoms)
    factors = sorted({a[0] for a in atoms})

    @lru_cache(None)
    def signatures(h,f,initial_mark):
        # Distinct signatures only; no producer word IDs or representatives.
        if h < 1 or f < h:
            return frozenset()
        out = set()
        for a in marked if initial_mark else atoms:
            H,P,F = a[:3]
            if F and (h,f) == (H,F):
                out.add(a)
            if h % H:
                continue
            tail_h = h // H
            numerator = f-F*tail_h
            if numerator <= 0 or numerator % P:
                continue
            tail_f = numerator // P
            for b in signatures(tail_h,tail_f,False):
                out.add(multiply(a,b))
        return frozenset(out)

    def h_values(bound):
        found = {1}
        queue = [1]
        for h in queue:
            for H in factors:
                q = h*H
                if q <= bound and q not in found:
                    found.add(q); queue.append(q)
        return sorted(found)

    def prefixes(c,R):
        # Values are literal (h,p) for one actual all-unmarked prefix; the
        # identity is explicitly empty, avoiding a false empty full signature.
        out = {(1,1): (1,1,True)}
        queue = [(1,1)]
        for h,o in queue:
            old_h,old_p,_ = out[(h,o)]
            for a in zero:
                H,P = a[:2]
                target = h*H,o*odd(P)
                if abs(c) % target[0] or R % target[1] or target in out:
                    continue
                out[target] = old_h*H,old_p*P,False
                queue.append(target)
        return out

    decisions = []
    normalizations = 0
    wider_searches = 0
    for test in receipt['tests']:
        c,R = test['c'],test['R']
        require(c != 0 and c % 2 and R > 0 and R % 2, 'target type')
        answers = set()
        prefix_states = prefixes(c,R)
        for (hp,op),(literal_h,literal_p,empty) in prefix_states.items():
            cc,rr = c//hp,R//op
            B = M*M+abs(cc)  # strictly no smaller than the producer envelope
            cap = max(max(a[2] for a in marked),
                      max(((a[2]+a[1])*B-a[1]*cc)//a[0] for a in marked))
            hs = h_values(B)
            f = rr
            while f <= cap:
                wider_searches += 1
                for h in hs:
                    if h > f: break
                    for sig in signatures(h,f,True):
                        if sig[0]-sig[4] != cc:
                            continue
                        full = sig if empty else (literal_h*sig[0],literal_p*sig[1],
                            literal_p*sig[2],literal_h*sig[3],literal_h*sig[4],
                            literal_h*sig[5],0)
                        require(full[0]-full[4] == c and odd(full[2]) == R,
                                'actual prefix reconstruction')
                        answers.add(full)
                f *= 2
        expected = {tuple(w['sigma']) for w in test['witnesses']}
        require(answers == expected, 'independent complete decision '+test['address'])
        for search in test['searches']:
            ph,po = search['prefix_h'],search['prefix_odd_p']
            require((ph,po) in prefix_states, 'producer prefix accepted independently')
            cc = c//ph
            tight_B = max(M*M,abs(cc))
            tight_cap = max(max(a[2] for a in marked),
                max(((a[2]+a[1])*tight_B-a[1]*cc)//a[0] for a in marked))
            require(search['normalized_c'] == cc and search['proved_f_bound'] == tight_cap,
                    'literal normalization receipt bound')
            normalizations += 1
        decisions.append((test['address'],len(answers)))

    accepted = ['',*(address for address,n in decisions if n)]
    require(accepted == ['', '2','4','22','42','422','4222'], 'entire returned tree')
    require(receipt['complete_tree'] is True and not receipt['unexpanded_states'],
            'complete exhaustion flag')
    declared = {s['address']:(s['c'],s['R']) for s in receipt['states']}
    require(set(declared) == set(accepted), 'returned state list')
    require(declared[''] == (-5,7) and literal_signature((6,4,0),7) == (1,8,7,7,6,1,1),
            'actual specified initial all-marked TT3')
    for address,(c,R) in declared.items():
        for digit in (2,4):
            t = [t for t in receipt['tests'] if t['address'] == address+str(digit)]
            require(len(t) == 1 and (t[0]['c'],t[0]['R']) == (3*c+digit,3*R),
                    'every state has both successor decisions')

    witnesses = 0
    for test in receipt['tests']:
        for witness in test['witnesses']:
            rows,ports = tuple(witness['rows']),witness['ports']
            sig = tuple(witness['sigma'])
            require(literal_signature(rows,ports) == sig, 'independent literal witness')
            for clock in (1,2):
                actual = attach(rows,ports,clock)
                root_observer = literal_signature(actual,1 << (len(actual)-1))
                require(strong(actual), 'actual full-graph witness strongness')
                require(root_observer[0] == sig[0]-sig[4]+(1 << clock)*sig[2],
                        'actual full-graph family count')
                require(root_observer[2]-root_observer[0] == sig[4],
                        'actual physical-root defect')
            witnesses += 1

    # Independent literal test of the new normalized valuation hostile.
    A = (24,1,3,6,14)
    atom = tuple(all_signatures(A))[16]
    require(atom == (11,100,12,27,0,1,4), 'valuation atom by permutations')
    require(strong(A), 'valuation atom strong')
    for m in range(1,8):
        n = m.bit_length()
        chain = tuple(sum(1 << j for j in range(i+1,n)) for i in range(n))
        chain_sig = tuple(all_signatures(chain))[m]
        sig = literal_signature(ordered_union(A,chain),16 | (m << 5))
        require(sig == multiply(atom,chain_sig), 'valuation literal product')
        require(sig[2] == 12+100*m and sig[4] == 27*m+11*chain_sig[4]+(m % 2),
                'valuation explicit formula')
        require(sig[0]-sig[4] <= 11-27*m, 'valuation changing intercept')
    units = {a for a in range(1,25) if a % 5}
    require({pow(2,a,25) for a in range(20)} == units, 'complete mod25 unit orbit')
    for level in range(21):
        R = 7*3**level
        a0 = next(a for a in range(20) if pow(2,a,25)*(R % 25) % 25 == 3)
        for period in (1,2,3):
            a = a0+20*period
            m = (2**a*R-3)//25
            require(m > 0 and 25*m == 2**a*R-3 and 12+100*m == 2**(a+2)*R,
                    'arbitrarily repeatable valuation arithmetic control')
    print('Independent all-clock-shift audit PASS')
    print('Strong labeled shape counts:',dict(sorted(shape_counts.items())))
    print('Permutation-derived full-signature bank:',len(atoms))
    print('Unmarked (h,p) types:',[(a[0],a[1]) for a in zero])
    print('Decisions / distinct full-signature witnesses:',decisions)
    print('Producer normalization searches:',normalizations)
    print('Independent searches under larger h-envelope:',wider_searches)
    print('Literal witnesses / full physical graphs:',witnesses,2*witnesses)
    print('Valuation hostile: seven literal masks; complete unit orbit; 63 arithmetic controls')
    print('Independent inverse cache:',signatures.cache_info())
    print('Receipt SHA256:',hashlib.sha256(raw).hexdigest())
    print('Always-active audit gates:',GATES)


if __name__ == '__main__':
    main()
