"""Exact width-two Collatz bridge and a recursively certified 41-tail family.

No floating arithmetic decides a barrier, residue, or orbit certificate.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import json


def require(value, label):
    if not value:
        raise RuntimeError(label)


def word_data(word):
    a, c, positive = 0, 0, True
    for j, bit in enumerate(word):
        if bit:
            a += 1
            c = 3*c + 2**j
        if 3**a <= 2**(j+1):
            positive = False
    modulus = 2**len(word)
    residue = (-c*pow(3**a, -1, modulus)) % modulus
    return a, c, residue, positive


def barrier(j):
    i = 0
    while 3**i <= 2**(i+j):
        i += 1
    return i


def poset_words(a, b):
    """Independent minimal-element recursion for two chains and cross edges."""
    thresholds = [0]+[barrier(j) for j in range(1, b+1)]
    def walk(i, j, prefix):
        if i == a and j == b:
            yield prefix
        if i < a:
            yield from walk(i+1, j, prefix+(1,))
        if j < b and i >= thresholds[j+1]:
            yield from walk(i, j+1, prefix+(0,))
    return list(walk(0, 0, ()))


def cross_counts(words, a, b):
    counts = [[0]*b for _ in range(a)]
    for word in words:
        ones, zeros = 0, 0
        for bit in word:
            if bit:
                ones += 1
            else:
                for i in range(ones):
                    counts[i][zeros] += 1
                zeros += 1
    return counts


def balance(words, a, b):
    total = len(words)
    return max((F(min(c, total-c), total) for row in cross_counts(words, a, b)
                for c in row), default=F(0))


def root_bridge():
    universes, swaps = 0, 0
    totals = []
    for length in range(1, 15):
        good_total = 0
        for a in range(length+1):
            b = length-a
            independent = []
            for positions in combinations(range(length), a):
                chosen = set(positions)
                word = tuple(int(j in chosen) for j in range(length))
                if word_data(word)[3]:
                    independent.append(word)
            extensions = poset_words(a, b)
            require(set(independent) == set(extensions), 'barrier/linear-extension bijection')
            good_total += len(extensions)
            universes += 1
            if len(extensions) > 1:
                require(balance(extensions, a, b) >= F(1, 3), 'width-two balancing control')
            for word in extensions:
                _, c, residue, _ = word_data(word)
                n = residue
                actual = []
                for _ in word:
                    bit = n % 2
                    actual.append(bit)
                    n = (3*n+1)//2 if bit else n//2
                require(tuple(actual) == word, 'actual residue realization')
                for j in range(length-1):
                    if word[j:j+2] != (1, 0):
                        continue
                    other = word[:j]+(0, 1)+word[j+2:]
                    if not word_data(other)[3]:
                        continue
                    _, other_c, other_residue, _ = word_data(other)
                    suffix_ones = sum(word[j+2:])
                    require(other_c-c == 2**j*3**suffix_ones, 'exact swap carry')
                    difference = (other_residue-residue) % (2**length)
                    require((difference & -difference) == 2**j, 'swap destroys fixed source at bit j')
                    swaps += 1
        totals.append(good_total)
    print(json.dumps({'fixed_count_universes': universes, 'max_length': 14,
                      'exact_swap_controls': swaps, 'W_1_through_14': totals}))
    words = poset_words(5, 1)
    records = sorted((''.join(map(str, w)), word_data(w)[2]) for w in words)
    chosen = [w for w in words if 1 <= word_data(w)[2] <= 31]
    require(set(word_data(w)[2] for w in chosen) == {27, 31}, '27 height-cut hostile')
    counts = cross_counts(chosen, 5, 1)
    # Every unanimous cross relation holds for all four full extensions.
    for w in words:
        own = cross_counts([w], 5, 1)
        require(all(counts[i][0] not in (0, len(chosen))
                    or own[i][0] == counts[i][0]//len(chosen) for i in range(5)),
                'intersection-order closure strictly enlarges height-cut law')
    print(json.dumps({'six_bit_words_and_least_residues': records,
                      'cutoff31_selected_words': [''.join(map(str,w)) for w in chosen],
                      'selected_law_is_full_LE_law_of_no_poset_on_these_labels': True}))


def lift_clock(depth):
    residues, r, modulus, period = [], 0, 3, 2
    require((41*pow(2, r, modulus)+1) % modulus == 0, 'first clock phase')
    for k in range(1, depth+1):
        require(0 <= r < period, 'least phase range')
        require((41*pow(2, r, modulus)+1) % modulus == 0, 'phase congruence')
        residues.append(r)
        candidates = [r+period*d for d in range(3)]
        good = [s for s in candidates if (41*pow(2, s, 3*modulus)+1) % (3*modulus) == 0]
        require(len(good) == 1, 'unique ternary lift')
        r = good[0]
        modulus *= 3
        period *= 3
    return residues


def family_count(X):
    """All k>=1 and r>=0 with source<=X, by integer logarithms only."""
    answer, distinct = 0, set()
    # In this finite probe retain sources too, independently checking collisions.
    phases = lift_clock((X+1).bit_length())
    for k, phase in enumerate(phases, 1):
        if 2**k > X+1:
            break
        period = 2*3**(k-1)
        cap = (3**k*(X+1)//2**k-1)//41
        if cap < 1:
            continue
        rmax = cap.bit_length()-1
        if phase > rmax:
            continue
        count = 1+(rmax-phase)//period
        answer += count
        for r in range(phase, rmax+1, period):
            u, rem = divmod(41*2**r+1, 3**k)
            require(rem == 0, 'family integrality')
            n = 2**k*u-1
            require(0 < n <= X, 'family cutoff')
            require(n not in distinct, 'unique family parameters')
            distinct.add(n)
    return answer, distinct


def root_family():
    phases = lift_clock(200)
    actual_rows = []
    for k, r in enumerate(phases[:9], 1):
        u = (41*2**r+1)//3**k
        n = 2**k*u-1
        x = n
        for j in range(k):
            require(x % 2 == 1, 'prescribed odd run')
            x = (3*x+1)//2
            require(x == 3**(j+1)*2**(k-j-1)*u-1 and x > n, 'strict growth formula')
        require(x == 41*2**r, 'exact 41 landing')
        require(x >> r == 41 and x % 2**r == 0, 'certified halving tail')
        period = 2*3**(k-1)
        if k <= 5:
            B = 2**period
            intercept = (B-1)*(3**k-2**k)//3**k
            next_n = 2**k*(41*2**(r+period)+1)//3**k-1
            require(next_n == B*n+intercept, 'exact within-family affine recursion')
        actual_rows.append({'k':k,'r_k':r,'n':str(n) if n.bit_length()<150 else 'bits='+str(n.bit_length()),
                            'total_shortcut_time':69+k+r})
    x, tail = 41, 0
    while x != 1:
        x = (3*x+1)//2 if x % 2 else x//2
        tail += 1
    require(tail == 69, 'known 41 tail length')
    require(actual_rows[0]['n']=='27' and actual_rows[1]['n']=='291'
            and actual_rows[2]['n']=='12439', 'initial long-growth family')
    print(json.dumps({'ternary_clock_lifts':200,'first_phases':phases[:12],
                      'actual_growth_family':actual_rows}))
    results = []
    for power in (10,20,40,80,160,320):
        total, sources = family_count(2**power)
        require(total == len(sources), 'union floor count')
        results.append([power, total])
    # Independent brute characterization of the rise-then-halve-to41 words.
    exact, sources = family_count(100000)
    brute = set()
    for n in range(1,100001,2):
        x, k = n, 0
        while x % 2:
            x = (3*x+1)//2
            k += 1
            if x == 41:
                brute.add(n)
            if x % 2 == 0:
                y = x
                while y % 2 == 0:
                    y //= 2
                if y == 41:
                    brute.add(n)
                break
    require(sources == brute, 'all sources<=100000 by independent orbit test')
    print(json.dumps({'counts_at_X_2_power':results,'brute_sources_up_to100000':sorted(brute)}))
    covers = []
    for J in (1,2,4,8,16,32,64,128,256):
        modulus = 2**J
        residues = {modulus-1}
        cover_phases = lift_clock(J)
        for k in range(1,J):
            inverse = pow(3**k,-1,modulus)
            q = (2**k*inverse-1) % modulus
            residues.add(q)
            for r in range(cover_phases[k-1],J-k,2*3**(k-1)):
                residues.add((q+41*pow(2,k+r,modulus)*inverse) % modulus)
        require(4*len(residues) <= 11*J, 'linear dyadic closure cover')
        covers.append([J,len(residues)])
    print(json.dumps({'closure_residue_counts_mod_2_power':covers}))
    print('PASS: exact bridge, carry loss, and certified convergent long-growth family; no Collatz proof.')


if __name__ == '__main__':
    root_bridge()
    root_family()
