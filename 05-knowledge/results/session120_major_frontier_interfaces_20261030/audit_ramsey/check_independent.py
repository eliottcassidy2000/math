"""Independent audit: no imports from producer; truth tables -> literal arcs -> correlations."""
import itertools
import json
import math
from collections import Counter
from pathlib import Path


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def primes(n):
    return [p for p in range(3, n+1, 2) if n % p == 0 and all(p % q for q in range(2, math.isqrt(p)+1))]


def response_tables(ps):
    # Independent enumeration by the integer binary cube, antipode is XOR all ones.
    size = 1 << len(ps)
    out = []
    for values in itertools.product((-1, 1), repeat=size//2):
        table = list(values) + [-v for v in reversed(values)]
        # A signed character is characterized by multiplicativity after normalization.
        character = all(table[a ^ b] * table[0] == table[a]*table[b]
                        for a in range(size) for b in range(size))
        out.append((table, character))
    return out


def options(n):
    result = []
    for m in range(2, n+1):
        if n % m:
            continue
        ps = primes(m)
        rows = []
        for table, character in response_tables(ps):
            mask = 0
            for t in range(1, m):
                if math.gcd(t,m) != 1:
                    continue
                # bit=1 denotes Legendre -1; order does not match producer enumeration.
                cube = sum((pow(t, (p-1)//2, p) == p-1) << i for i,p in enumerate(ps))
                if table[cube] == 1:
                    mask |= 1 << ((n//m)*t)
            rows.append((mask, character))
        result.append(rows)
    return result


def correlation_moment(n, mask):
    full = (1 << n)-1
    require(mask.bit_count() == (n-1)//2 and mask & 1 == 0, 'invalid tournament support')
    total = (n-1)**2
    for h in range(1,n):
        rotated = ((mask << h) & full) | (mask >> (n-h))
        c = 4*(mask & rotated).bit_count() - n + 2
        total += c*c
    return n*total


def scalar_moment(n, mask):
    s = [0] + [1 if mask >> i & 1 else -1 for i in range(1,n)]
    require(all(s[i] == -s[-i] for i in range(1,n)), 'antipodality failed')
    return n*sum(sum(s[x]*s[(x+h)%n] for x in range(n))**2 for h in range(n))


def census(n):
    hist = {'character': Counter(), 'nonlinear': Counter()}
    witnesses = {}
    for rows in itertools.product(*options(n)):
        mask = 0
        for part, _ in rows:
            require(mask & part == 0, 'strata overlap')
            mask |= part
        kind = 'character' if all(c for _,c in rows) else 'nonlinear'
        moment = correlation_moment(n,mask)
        if n < 50:
            require(scalar_moment(n,mask) == moment, f'all-small scalar mismatch N={n}')
        hist[kind][moment] += 1
        if kind not in witnesses or moment < witnesses[kind][0]:
            witnesses[kind] = (moment,mask)
    for moment,mask in witnesses.values():
        require(scalar_moment(n,mask) == moment,'literal scalar minimum mismatch')
    return hist


def main():
    base = Path(__file__).resolve().parent.parent
    receipt = json.loads((base/'ramsey'/'receipt.json').read_text())
    h = census(231)
    result = {'status':'PASS', 'producer_imports':False, 'N':231, 'summaries':{}, 'controls':[]}
    for kind, counts in h.items():
        actual = {str(k):v for k,v in sorted(counts.items())}
        require(actual == receipt['census']['summaries'][kind]['histogram'],kind+' full histogram differs')
        result['summaries'][kind] = {'count':sum(counts.values()),'minimum':min(counts), 'minimum_multiplicity':counts[min(counts)],'maximum':max(counts),'distinct_moments':len(counts),'histogram':actual}
    require(sum(sum(c.values()) for c in h.values()) == 8192,'family count')
    for control in receipt['controls']:
        n = control['N']
        counts = census(n)
        joined = counts['character'] + counts['nonlinear']
        actual = {str(k):v for k,v in sorted(joined.items())}
        require(actual == control['moments'],f'control mismatch N={n}')
        result['controls'].append({'N':n,'moments':actual,'count':sum(joined.values())})
    print(json.dumps(result, sort_keys=True, indent=2))

if __name__ == '__main__':
    main()
