#!/usr/bin/env python3
"""Independent exact controls. No imports from the submitted verifier.

Proofs of universal claims remain in the review; this program supplies
independently implemented arithmetic, cylinder and rational-floor controls.
"""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
from math import gcd
from pathlib import Path
import json

checks = 0

def check(test, label):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(label)

def atanh_log(z, terms):
    """Bounds for log((1+z)/(1-z)), with iterative powers."""
    partial, power = Q(0), z
    square = z*z
    for j in range(terms):
        partial += 2*power/(2*j+1)
        power *= square
    return partial, partial + 2*power/((2*terms+1)*(1-square))

# Deliberately different logarithmic basis and truncation from submission:
# log 2 from z=1/3, log(3/2) from z=1/5, each with 192 terms.
l2, u2 = atanh_log(Q(1,3), 192)
l15, u15 = atanh_log(Q(1,5), 192)
check(l2 > u15, 'positive denominator for direct logarithmic ratio')
tlo = (77*l15-45*u2)/(u2-l15)
thi = (77*u15-45*l2)/(l2-u15)
check(Q(1,10) < tlo < thi < Q(1,9), 'independent rotation bracket')

K = '110110110110110101011011101010'
words = {'K':K, 'V':K*4+'10', 'F':'10', 'H':'110101000'}

def coefficients(w):
    # Closed sum b = sum_{odd position j} 2^j 3^(number of later odd bits).
    power, carry = 1, 0
    for position in range(len(w)-1, -1, -1):
        if w[position] == '1':
            carry += (1 << position)*power
            power *= 3
    return power, carry, 1 << len(w)

def valuation3(z):
    if not z:
        return None
    z, answer = abs(z), 0
    while z % 3 == 0:
        z //= 3
        answer += 1
    return answer

def replay(x, w):
    for bit in w:
        check(x % 2 == int(bit), 'literal parity authenticated by actual orbit')
        x = (3*x+1)//2 if x % 2 else x//2
    return x

bank = {name:coefficients(w) for name,w in words.items()}
bank_report = {}
for name, expected in [('K',3), ('V',4), ('H',1)]:
    a,b,M = bank[name]
    h = words[name].count('1')
    check(valuation3(b-M) == expected < h, 'universal output valuation algebra')
    source = -b*pow(a,-1,M) % M
    for lift in range(6):
        x = source + lift*M
        if x == 0:
            x = M
        y = replay(x, words[name])
        check(M*y == a*x+b, 'closed-sum branch agrees with actual orbit')
        check(valuation3(y-1) == expected, 'lawful output type control')
        check((M*y-b) % (3**h) == 0, 'full image congruence control')
    bank_report[name] = {'odd_bits':h, 'type':expected, 'a':a, 'b':b, 'M':M}

check(bank['F'] == (3,1,4), 'exact F branch and type increment identity')
for j in range(160):
    x = 1+4*j
    y = replay(x, words['F'])
    check(4*(y-1) == 3*(x-1), 'F increment identity')
    old, new = valuation3(x-1), valuation3(y-1)
    check(new == (None if old is None else old+1), 'F finite and infinite type controls')

modulus = 243
a,b,M = bank['V']
rho_v = b*pow(M,-1,modulus) % modulus
avf,bvf,Mvf = coefficients(words['V']+words['F'])
rho_vf = bvf*pow(Mvf,-1,modulus) % modulus
check((b % 243, M % 243, pow(M,-1,243)) == (103,22,232), 'V residue constituents')
check((rho_v,rho_vf) == (82,1), 'independent V versus VF residues')
check(valuation3(bvf-Mvf) == 5, 'direct VF type agrees with increment proof')
for r in range(163):
    residue = (rho_v-rho_vf)*pow(1 << r,-1,243) % 243
    check(residue in (81,162), 'full modular power cycle for factor 81')
    check(min(residue,243-residue) == 81, 'minimum CRT distance for image pair')

# Generic CRT controls include negative k, rather than only positive differences.
for s in range(1,6):
    modulus = 3**s
    for residue in range(1,modulus):
        distance = min(residue,modulus-residue)
        candidates = [residue+j*modulus for j in range(-3,4)]
        check(min(map(abs,candidates)) == distance, 'nearest signed representative')
        for r in (0,1,2,7,31):
            check(all(abs((1 << r)*k) >= (1 << r)*distance for k in candidates),
                  'signed CRT lower bound')

small = (1003,9885,1116187)
large = (8116910408926339217391896275524023351,
         79995672446647253540458867489349070728,
         15001586057656515185428419394326237781535)

def certify(p,q,m):
    delta_lo, delta_hi = Q(p,q)-thi, Q(p,q)-tlo
    check(0 < p < q and q >= 2 and gcd(p,q) == 1, 'reduced rational approximation')
    check((m*p-1) % q == 0 and m > 2*q, 'inverse congruence and room')
    check(Q(1,m*q) < delta_lo < delta_hi < Q(1,(m-q)*q), 'strict logarithmic bracket')
    check(m*delta_hi < Q(2,q), 'residue separation bound')
    floor_lo, floor_hi = (m*tlo).__floor__(), (m*thi).__floor__()
    check(floor_lo == floor_hi, 'independent exact final floor')
    check(floor_lo == (m*p-1)//q-1, 'exceptional residue floor at m')
    L = 122*m
    N = L + 2*(floor_lo+1)
    P = 122*q + 2*p
    r = (L+99)//100
    cycles, remainder = divmod(m-1,q)
    H = lambda n:n*p//q-int(n%q == 0)
    remainder_ones = H(remainder+1)-H(1)
    check(122+cycles*P+122*remainder+2*remainder_ones == N, 'compressed exact word length')
    check(122+P+r < N, 'both matching factors fit')
    check(p > q*thi, 'negative periodic multiplier logarithm')
    return {'p':p,'q':q,'m':m,'L':L,'N':N,'P':P,'r':r,
            'full_tail_periods':cycles,'remainder_images':remainder,
            'last_floor':floor_lo,'all_long_factors_hostile':r >= 122+P}

controls = [certify(*small), certify(*large)]
check(controls[0]['last_floor'] == 113255, 'small displayed floor')
check(controls[0]['N'] == 136401326, 'small displayed length')
check(controls[0]['N'] < 2**52, 'small inside previous short range')
check(controls[1]['last_floor'] == 1522163965844627930660358762978517428452,
      'enormous displayed floor')
check(controls[1]['N'] == 1833237826965784108483587883633758044204176,
      'enormous displayed length')
check(controls[1]['N'] > 2**52, 'enormous beyond previous short range')
check(all(c['all_long_factors_hostile'] for c in controls), 'both whole-future hostiles')
check(controls[1]['full_tail_periods'] == 187, 'enormous full period count')

# Independent derivation of 39 upper convergents and their next complete-quotient digits.
u,v = tlo,thi
pm,pn,qm,qn = 0,1,1,0
pending = None
convergents = []
for index in range(79):
    digit = u.__floor__()
    check(digit == v.__floor__(), 'certified continued-fraction digit')
    if pending is not None:
        p,q,prior_q,prior_index = pending
        m = (digit+1)*q+prior_q
        record = certify(p,q,m)
        record['index'] = prior_index
        convergents.append(record)
        pending = None
    p,q = digit*pn+pm, digit*qn+qm
    if index % 2:
        check(p*qn-pn*q == 1, 'upper determinant sign')
        check(Q(p,q) > thi, 'upper orientation')
        check(qn > 0, 'positive preceding denominator ensures m greater than 2q')
        pending = p,q,qn,index
    pm,pn,qm,qn = pn,p,qn,q
    u,v = 1/(v-digit),1/(u-digit)
check(len(convergents) == 39, 'upper convergent controls count')
check(tuple(convergents[3][k] for k in ('p','q','m')) == small, 'small convergent recovered')
check(tuple(convergents[-1][k] for k in ('p','q','m')) == large, 'large convergent recovered')
check(any(not c['all_long_factors_hostile'] for c in convergents),
      'unbounded descent family does not silently assert universal hostile ratio')

# Exhaust rational instances of the floor argument, independently of theta.
# A rational interior value has exactly the same required floor geometry.
rational_instances = 0
for q in range(2,31):
    for p in range(1,q):
        if gcd(p,q) != 1:
            continue
        inverse = pow(p,-1,q)
        for t in range(2,7):
            m = t*q+inverse
            delta = (Q(1,m*q)+Q(1,(m-q)*q))/2
            theta = Q(p,q)-delta
            H = lambda n:n*p//q-int(n%q == 0)
            model = [H(n)-H(n-1) for n in range(1,q+1)]
            check(model[0] == 1 and model[-1] == 0, 'rational model boundary pair')
            check(set(model) <= {0,1} and sum(model) == p, 'binary period and one count')
            prior = 0
            actual = []
            for n in range(1,m+1):
                floor = (n*theta).__floor__()
                check(floor == H(n)-int(n == m), 'complete rational floor identity')
                actual.append(floor-prior)
                prior = floor
            periodic = [model[n%q] for n in range(1,m)]
            check(actual[0] == actual[-1] == 0, 'initial and final reset losses')
            restored = actual[1:-1]+[actual[-1]+1]
            check(restored == periodic, 'entire final-reset-corrected periodic tail')
            rational_instances += 1

# Exhaust short binary preperiod/period choices for the universal continuation
# lemma and its exact frequency conclusion, including nonprimitive periods.
periodic_extension_pairs = 0
periodic_extension_examples = 0
for prelength in range(4):
    for pretuple in product('01',repeat=prelength):
        A = ''.join(pretuple)
        for periodlength in range(1,6):
            for periodtuple in product('01',repeat=periodlength):
                B = ''.join(periodtuple)
                bit = lambda pos:A[pos] if pos < len(A) else B[(pos-len(A))%len(B)]
                r = len(A)+len(B)
                for i in range(len(A)+2*len(B)+1):
                    for j in range(i+1,len(A)+2*len(B)+1):
                        if any(bit(i+t) != bit(j+t) for t in range(r)):
                            continue
                        check(all(bit(i+t) == bit(j+t) for t in range(r+3*len(B))),
                              'matched periodic-tail factors have identical futures')
                        distance = j-i
                        separator_ones = sum(bit(i+t) == '1' for t in range(distance))
                        check(separator_ones*len(B) == distance*B.count('1'),
                              'separator has the exact tail frequency')
                        periodic_extension_pairs += 1
                periodic_extension_examples += 1

short_word = '0'+'01'*5
check(short_word[0] == short_word[1] and short_word[1] != short_word[2],
      'below-threshold equal-factor hostile retains a later mismatch')

# Whole-word lawful-cylinder controls for the first two convergent family members.
whole_word_controls = []
for control in convergents[:2]:
    p,q,m = (control[k] for k in ('p','q','m'))
    floors = [0]
    for n in range(1,m+1):
        lo,hi = (n*tlo).__floor__(),(n*thi).__floor__()
        check(lo == hi, 'all floors of short true-theta word certified')
        floors.append(lo)
    reset = [floors[n]-floors[n-1] for n in range(1,m+1)]
    word = ''.join(words['V']+words['F']*bit for bit in reset)+words['F']
    H = lambda n:n*p//q-int(n%q == 0)
    B = ''.join(words['V']+words['F']*(H(n)-H(n-1)) for n in range(2,q+2))
    check(B.endswith(words['V']+words['F']), 'period ends at VF image type')
    aB,bB,MB = coefficients(B)
    check(aB < MB and bB*pow(MB,-1,243) % 243 == 1,
          'exact short contracting period has VF output residue')
    check(word[122:] == (B*((len(word)-122)//len(B)+1))[:len(word)-122],
          'literal final forced-reset periodicity')
    r,P = control['r'],control['P']
    check(word[122:122+r] == word[122+P:122+P+r], 'literal chosen future match')
    a,b,M = coefficients(word)
    check(a < M, 'whole greedy word contracts')
    source = -b*pow(a,-1,M) % M
    for lift in range(3):
        x = source+lift*M
        check(x > 0, 'positive cylinder lift')
        y = replay(x,word)
        check(M*y == a*x+b, 'whole closed-sum affine formula authenticates replay')
        check(y < x, 'actual endpoint descent in whole-word controls')
        yi = replay(x,word[:122])
        yj = replay(x,word[:122+P])
        check(valuation3(yi-1) == 4 and valuation3(yj-1) == 5,
              'actual selected checkpoint types differ')
        check((yi-yj) % (1 << r) == 0 and abs(yi-yj) >= 81*(1 << r),
              'actual equal-future CRT separation')
    whole_word_controls.append({'m':m,'N':len(word),'three_lifts_descend':True})

base = Path(__file__).resolve().parent.parent
manifest = json.loads((base/'MANIFEST.json').read_text())
hashes = {}
for entry in manifest['files']:
    data = (base/entry['path']).read_bytes()
    digest = sha256(data).hexdigest()
    check(len(data) == entry['bytes'] and digest == entry['sha256'], 'frozen manifest verification')
    hashes[entry['path']] = digest
check(hashes['IMAGE_GUARD_DESCENT.md'] ==
      '8e32d46d23eceefee6a45b2bc7c00548d5409cc14dda99171298d0d56dce3c84',
      'assigned primary frozen hash')

print(json.dumps({'status':'PASS independent image-guard controls',
    'checks':checks,'log_terms':192,'log_basis':['log(2)','log(3/2)'],
    'bank':bank_report,'image_residues_mod_243':{'V':rho_v,'VF':rho_vf},
    'crt_factor':81,'exact_compressed_controls':controls,
    'upper_convergent_controls':len(convergents),
    'upper_controls_meeting_stronger_hostile_ratio':sum(c['all_long_factors_hostile'] for c in convergents),
    'rational_floor_instances':rational_instances,'whole_word_controls':whole_word_controls,
    'periodic_extension_examples':periodic_extension_examples,
    'matched_periodic_extension_pairs':periodic_extension_pairs,
    'frozen_source_hashes':hashes,'giant_word_materialized':False,
    'universal_proofs_supplied_in_review':True},indent=2,sort_keys=True))
