"""Exact finite affine representations and the positivity-refinement boundary.

Standard library only. No output files are modified by this program.
"""
from collections import deque
from fractions import Fraction as F
from math import gcd
import json


def natural(n, minimum=0):
    if type(n) is not int or n < minimum:
        raise ValueError('exact integer outside domain')
    return n


def modulus(m):
    natural(m, 2)
    if gcd(m, 6) != 1:
        raise ValueError('modulus must be coprime to 6')
    return m


def odd(n):
    natural(n, 1)
    if not n % 2:
        raise ValueError('positive odd source required')
    return n


def step(n):
    odd(n)
    raw = 3*n+1
    a = (raw & -raw).bit_length()-1
    return raw >> a, a


def order2(m):
    modulus(m)
    x, h = 2 % m, 1
    while x != 1:
        x = 2*x % m
        h += 1
    return h


def affine_letter(m, a):
    modulus(m)
    natural(a, 1)
    inv = pow(pow(2, a, m), -1, m)
    return 3*inv % m, inv


def compose(g, h, m):
    """g after h; internal affine pairs have already been constructed."""
    return g[0]*h[0] % m, (g[0]*h[1]+g[1]) % m


def affine_group(m):
    modulus(m)
    generators = (affine_letter(m, 1), affine_letter(m, 2))
    seen, queue = {(1, 0)}, deque([(1, 0)])
    while queue:
        g = queue.popleft()
        for a in generators:
            nxt = compose(a, g, m)
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return seen


def multiplier_group(m):
    modulus(m)
    seen, queue = {1}, deque([1])
    while queue:
        a = queue.popleft()
        for p in (2, 3):
            nxt = a*p % m
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return seen


def word_data(word):
    if type(word) is not tuple:
        raise ValueError('exact tuple of positive valuations required')
    P, Q, B = 1, 1, 0
    for a in word:
        natural(a, 1)
        P, Q, B = 3*P, (1 << a)*Q, 3*B+Q
    return P, Q, B


def lift_word(m, residue, word, height=1):
    """Actual source with given residue and exact word; every node > height.

    Empty words are allowed. This constructs an ancestor, never identifies
    a requested source or silently supplies a ROOT certificate.
    """
    modulus(m)
    natural(residue)
    natural(height, 1)
    if residue >= m:
        raise ValueError('canonical residue required')
    P, Q, B = word_data(word)
    binary_modulus = 2*Q
    binary_residue = (Q-B)*pow(P, -1, binary_modulus) % binary_modulus
    # Final oddness gives the complete exact valuation word.
    t = (residue-binary_residue)*pow(binary_modulus, -1, m) % m
    n = binary_residue+binary_modulus*t
    period = m*binary_modulus
    threshold = height*Q+1
    if n < threshold:
        n += ((threshold-n+period-1)//period)*period
    return n, period


def replay_prefix(n, word):
    odd(n)
    word_data(word)
    path = [n]
    for a in word:
        n, actual = step(n)
        if actual != a:
            raise ValueError('false actual valuation')
        path.append(n)
    return tuple(path)


def root_word(n, limit=512):
    """Independent finite control only; used on the explicitly bounded head."""
    odd(n)
    natural(limit)
    out = []
    for _ in range(limit+1):
        if n == 1:
            return tuple(out)
        if len(out) == limit:
            break
        n, a = step(n)
        out.append(a)
    raise ValueError('finite control did not terminate within its declared cap')


def counter_weight(word):
    word_data(word)
    if not word:
        return F(1)
    L = len(word)-1
    K = sum((a-1)//2 for a in word)
    # 2 B(K+1,L+2), computed by an elementary product.
    out = F(2, (K+1)*(K+2))
    for j in range(L):
        out *= F(j+2, K+j+3)
    return out


def hostile_residue_mass(m, r):
    """q(3)=0; q(n)=2^-n for other positive odd n. Exact full residue sum."""
    natural(m, 3)
    if not m % 2:
        raise ValueError('odd modulus required')
    natural(r)
    if r >= m:
        raise ValueError('canonical residue required')
    a = r if r % 2 else r+m
    value = F(1, 1 << a)/(1-F(1, 1 << (2*m)))
    if r == 3 % m:
        value -= F(1, 8)
    return value


def main():
    checks = 0
    def check(ok, label):
        nonlocal checks
        checks += 1
        if not ok:
            raise ValueError(label)

    moduli = (5, 7, 11, 19, 35, 95)
    report, aliases = {}, []
    for m in moduli:
        group = affine_group(m)
        slopes = multiplier_group(m)
        check(group == {(a,b) for a in slopes for b in range(m)}, 'correct global affine group')
        matrices = {tuple((a*r+b) % m for r in range(m)) for a,b in group}
        check(len(matrices) == len(group), 'faithful permutation representation')
        for a,b in group:
            check(gcd(a,m) == 1, 'unit affine slope')
        remaining, frequency_orbits = set(range(m)), []
        while remaining:
            freq = min(remaining)
            orbit = {a*freq % m for a in slopes}
            check(orbit <= remaining, 'disjoint multiplier-frequency orbit')
            remaining -= orbit
            frequency_orbits.append(tuple(sorted(orbit)))
        h = order2(m)
        check(affine_letter(m,1) == affine_letter(m,h+1), 'same full finite matrix')
        check(F(3,2)>1>F(3,2**(h+1)), 'opposite real drift')
        for residue in range(m):
            left,_ = lift_word(m,residue,(1,))
            right,_ = lift_word(m,residue,(h+1,))
            lp = replay_prefix(left,(1,)); rp = replay_prefix(right,(h+1,))
            check(left % m == right % m == residue, 'actual same source residue')
            check(lp[-1] % m == rp[-1] % m, 'same finite affine image')
            check(left != right, 'actual disjoint valuation cells')
            children = []
            for a in (1,2):
                source,_ = lift_word(m,residue,(a,))
                children.append(step(source)[0] % m)
            check((children[0] == children[1]) == ((3*residue+1)%m == 0), 'exact residue closure boundary')
        for a in range(1,h+1):
            slope,carry = affine_letter(m,a)
            # Fourier character identity, tested as exact root-of-unity exponents.
            for freq in range(m):
                for residue in (0,1,m-1):
                    check(freq*(slope*residue+carry)%m == (freq*carry+(freq*slope % m)*residue)%m,
                          'Fourier phase and frequency transport')
        aliases.append({'m':m,'order2':h,'valuation_pair':[1,h+1],
                        'slopes':[str(F(3,2)),str(F(3,2**(h+1)))]})
        report[str(m)] = {'slopes':len(slopes),'affine_group':len(group),
                          'irreducible_permutation_dimensions':sorted(map(len,frequency_orbits)),
                          'residue_fibres_with_ambiguous_next_residue':m-1}

    check(len(multiplier_group(95)) == 36, 'inherited correlated CRT hostile')
    check(len(multiplier_group(5))*len(multiplier_group(19)) == 72, 'local product strictly larger')
    check(report['11']['irreducible_permutation_dimensions'] == [1,10], 'level11 finite representation dimensions')
    check(report['95']['irreducible_permutation_dimensions'] == [1,4,18,36,36], 'correlated CRT representation dimensions')

    words = ((),(1,),(2,),(1,2),(2,1),(1,1,2),(3,2,4),(1,2,1,1,3))
    lifts = 0
    for m in moduli:
        for word in words:
            P,Q,B = word_data(word)
            for residue in range(m):
                source,period = lift_word(m,residue,word,17)
                for offset in (0,1,3):
                    n = source+offset*period
                    path = replay_prefix(n,word)
                    check(n % m == residue and min(path)>17, 'arbitrarily high exact path lift')
                    check(Q*path[-1] == P*n+B, 'ordered carry retained')
                    lifts += 1

    # Exact color-resolved pushforward of a declared finite measure.
    # Actual root is included here; killed-root operators need a separate flag.
    flux_controls = 0
    for m in moduli:
        h = order2(m)
        actual = [F(0) for _ in range(m)]
        colored = {}
        for n in range(1,256,2):
            target,a = step(n)
            mass = F(1,n*(n+2))
            actual[target % m] += mass
            key = (n % m,a % h)
            colored[key] = colored.get(key,F(0))+mass
        reconstructed = [F(0) for _ in range(m)]
        for (r,c),mass in colored.items():
            a,b = affine_letter(m,c or h)
            reconstructed[(a*r+b)%m] += mass
        check(actual == reconstructed, 'color-resolved exact one-step mass transport')
        flux_controls += 1

    # Grounded controls: every residue in the finite modulus list has an
    # explicitly replayed ROOT representative. The all-modulus result is proved
    # by the inherited guarded finite-group lifting theorem, not this census.
    grounded = 0
    for m in moduli:
        masses = [F(0) for _ in range(m)]
        for n in range(1,2*m,2):
            word = root_word(n)
            check(replay_prefix(n,word)[-1] == 1, 'literal finite ROOT control')
            masses[n % m] += counter_weight(word)
        check(all(x>0 for x in masses), 'positive grounded mass in every residue')
        grounded += m

    hostile_levels = (5,7,11,19,35,95,191)
    for m in hostile_levels:
        masses = [hostile_residue_mass(m,r) for r in range(m)]
        check(all(x>0 for x in masses), 'every finite quotient strictly positive')
        check(sum(masses,F(0)) == F(13,24), 'exact hostile total mass')
        check(masses[3] == F(1,8*((1 << (2*m))-1)), 'positive residue mass tends to zero at missing atom3')

    bad = [lambda: modulus(True),lambda: modulus(5.0),lambda: modulus(9),lambda: order2(2),
           lambda: lift_word(5,0,(True,)),lambda: lift_word(5,0,[1]),
           lambda: lift_word(5,5,(1,)),lambda: lift_word(5,0,(1,),0),
           lambda: replay_prefix(True,()),lambda: replay_prefix(1.0,()),lambda: replay_prefix(2,()),
           lambda: hostile_residue_mass(6,0),lambda: hostile_residue_mass(5,False)]
    for job in bad:
        try:
            job()
        except ValueError:
            check(True,'invalid exact API input rejected')
        else:
            raise ValueError('invalid input accepted')

    print(json.dumps({'status':'PROVED algebra and scopes; FINITE-EXACT controls',
        'groups':report,'same_matrix_opposite_drift':aliases,
        'actual_prefix_lifts':lifts,'color_resolved_flux_controls':flux_controls,
        'grounded_residue_representatives':grounded,
        'hostile_missing_atom':3,'hostile_total_mass':'13/24',
        'hostile_moduli':list(hostile_levels),'invalid_inputs':len(bad),
        'checks':checks},indent=2,sort_keys=True))
    print('PASS: finite representation is faithful to affine residue maps; source positivity still needs an unbounded refinement lower bound.')


if __name__ == '__main__':
    main()
