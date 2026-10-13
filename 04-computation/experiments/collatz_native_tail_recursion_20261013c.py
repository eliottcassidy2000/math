"""Native tail payment and supplied-terminal recursion; no ROOT discovery API."""
from dataclasses import dataclass, replace
import hashlib
import collatz_adaptive_return_budget_20261013a as previous

old = previous.old


def need(ok, why):
    if not ok:
        raise ValueError(why)


@dataclass(frozen=True)
class Half:
    source: int
    word: tuple


def audit_half(receipt):
    need(type(receipt) is Half, 'exact half receipt')
    old.odd(receipt.source), old.word(receipt.word)
    need(receipt.source > 1 and bool(receipt.word), 'nonroot source and nonempty word')
    h = old.replay(receipt.source, receipt.word)
    need(h % 3 != 0, 'every nonempty actual endpoint is a ternary unit')
    need(2*(h+1) <= receipt.source+1, 'actual half-band endpoint')
    p, q, b = old.carrier(receipt.word)
    need(q > 2*p and p*receipt.source+b == q*h, 'retained positive carry and slope')
    return h, p, q, b


def lift(receipt, parameter):
    """The whole coarse native AP, with its actual final valuation retained."""
    h, p, q, _ = audit_half(receipt)
    old.natural(parameter)
    source = receipt.source+q*parameter
    tail = h+p*parameter
    extra = old.valuation(tail)
    word = receipt.word[:-1]+(receipt.word[-1]+extra,)
    endpoint = tail >> extra
    need(old.replay(source, word) == endpoint, 'exact coarse native lift')
    need(2*(endpoint+1) <= source+1, 'half-band monotonicity')
    charge = sum(a % 2 for a in word)
    if parameter:
        need(charge <= previous.E(source)**2-previous.E(endpoint)**2,
             'strict native tail pays coefficient one')
    return Half(source, word), endpoint, charge


def surplus(source, word, rate=3):
    """Exact telescoping credit for any actual segment, not only half exits."""
    old.odd(source), old.word(word), previous.rate_value(rate)
    endpoint = old.replay(source, word)
    return rate*(previous.E(source)**2-previous.E(endpoint)**2)-sum(a % 2 for a in word)


def splice(first, continuation):
    """Keep the original source and all debt through a supplied continuation."""
    h, _, _, _ = audit_half(first)
    old.word(continuation)
    old.replay(h, continuation)
    combined = Half(first.source, first.word+continuation)
    audit_half(combined)
    return combined


@dataclass(frozen=True)
class Grounded:
    receipt: Half
    child_word: tuple


def audit_family(family):
    need(type(family) is Grounded, 'exact grounded family')
    h, p, q, b = audit_half(family.receipt)
    old.word(family.child_word)
    need(old.replay(h, family.child_word) == 1, 'supplied exact-child first-hit ROOT')
    depth = len(family.receipt.word)
    period = 2*3**(depth-1)
    return h, p, q, b, depth, period


def member_word(family, index):
    """A symbolic valuation certificate: the large final letter is one integer."""
    *_, period = audit_family(family)
    old.natural(index)
    w = family.receipt.word
    return w[:-1]+(w[-1]+period*index,)+family.child_word


def source_mod(family, index, modulus):
    """Read the source using the full division precision P*modulus."""
    h, p, q, b, _, period = audit_family(family)
    old.natural(index)
    need(type(modulus) is int and modulus >= 1, 'positive exact modulus')
    denominator_modulus = p*modulus
    numerator = (q*h*pow(2, period*index, denominator_modulus)-b) % denominator_modulus
    need(numerator % p == 0, 'symbolic source integrality')
    return (numerator//p) % modulus


def materialize(family, index, bit_cap=8192):
    """Optional bounded expansion, separate from the exact symbolic certificate."""
    h, p, q, b, _, period = audit_family(family)
    old.natural(index), old.natural(bit_cap)
    exponent = period*index
    need(sum(family.receipt.word)+h.bit_length()+exponent <= bit_cap,
         'declared conservative materialization cap')
    numerator = (q*h << exponent)-b
    need(numerator % p == 0, 'exact materialized division')
    source = numerator//p
    word = member_word(family, index)
    need(old.replay(source, word) == 1, 'strict materialized ROOT certificate')
    return source, word


def recognize(family, source):
    """Total membership for a supplied finite integer; no orbit search."""
    h, p, q, b, _, period = audit_family(family)
    old.odd(source)
    numerator, denominator = p*source+b, q*h
    if numerator % denominator:
        return None
    ratio = numerator//denominator
    if ratio < 1 or ratio & (ratio-1):
        return None
    exponent = ratio.bit_length()-1
    return exponent//period if exponent % period == 0 else None


def extend_family(family, source, prefix):
    """Prepend an authenticated strict descent; its child address scales by 3^j."""
    audit_family(family)
    old.odd(source), old.word(prefix)
    need(bool(prefix) and source > family.receipt.source, 'strict source-ranked extension')
    need(old.replay(source,prefix) == family.receipt.source, 'exact supplied descent to factory seed')
    result = Grounded(Half(source,prefix+family.receipt.word),family.child_word)
    audit_family(result)
    return result


SHARED_CHILD = 16504781
# Discovered once at the exact named endpoint, with a 10,000-odd-step cap.
# Production only authenticates this frozen 80-letter first-hit certificate.
SHARED_WORD = (3,3,2,1,4,1,1,1,2,1,2,1,1,3,2,1,1,2,2,1,3,5,1,4,2,4,2,1,3,1,1,4,1,1,1,1,2,1,2,2,5,3,1,1,1,2,2,1,2,1,1,2,1,1,1,2,3,1,1,2,1,2,1,1,1,1,1,3,1,1,1,4,2,2,4,3,1,1,5,4)


def hostile_families():
    """The eleven explicitly declared sources, paid prefix discovery only."""
    sources = []
    for n in previous.record_sources():
        r = previous.search(n)
        _, _, first = previous.audit_return(r)
        if not previous.paid(n, first[1], first[2], 3):
            sources.append(n)
    sources.extend((96883183, 129177577))
    result = []
    for n in sorted(sources):
        r = previous.search(n)
        h, _, _ = previous.audit_return(r)
        suffix = SHARED_WORD if h == SHARED_CHILD else previous.root_suffix(h)
        family = Grounded(Half(n, r.word), suffix)
        audit_family(family)
        result.append(family)
    need(len(result) == 11, 'declared eleven-source universe')
    return tuple(result)


def main():
    checks = 0
    def check(ok, why):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(why)
    def reject(fn):
        try:
            fn()
        except ValueError:
            check(True, 'hostile rejected')
        else:
            check(False, 'hostile accepted')

    # Whole native tails from actual first-half receipts, not a convergence oracle.
    cell_count = 0
    for n in range(3, 1024, 2):
        paid = previous.search(n)
        _, _, first = previous.audit_return(paid)
        receipt = Half(n, paid.word[:first[0]])
        h, p, q, b = audit_half(receipt)
        for t in (0,1,2,3,8,17,2**30+1):
            lifted, endpoint, charge = lift(receipt, t)
            expected_numerator = p*lifted.source+b
            check(expected_numerator == (q*h+p*q*t), 'independent full carrier')
            check(endpoint == old.replay(lifted.source, lifted.word), 'strict native control')
            check(not t or surplus(lifted.source, lifted.word, 1) >= 0, 'coefficient-one invoice')
            cell_count += 1
    # The original source is a genuine exception, not licensed by its native tail.
    paid = previous.search(56924955)
    h, _, first = previous.audit_return(paid)
    first_half = Half(paid.source, paid.word[:first[0]])
    continuation = paid.word[first[0]:]
    intermediate, _, _, _ = audit_half(first_half)
    check(surplus(first_half.source,first_half.word,3) == -4, 'source deficit retained')
    check(surplus(intermediate,continuation,3) == 141, 'continuation pays retained deficit')
    complete = splice(first_half,continuation)
    check(complete.word == paid.word and surplus(complete.source,complete.word,3) == 137,
          'surplus cocycle exact')
    check(not previous.paid(first_half.source,intermediate,sum(a%2 for a in first_half.word),1),
          'no backwards transfer from native tail')
    for t in (1,2,17):
        r, _, _ = lift(first_half,t)
        check(surplus(r.source,r.word,1) >= 0, 'unpaid seed still has paid tail')

    # Bounded materializations independently replay the symbolic recursion.
    small = (Grounded(Half(3,(1,4)),()), Grounded(Half(5,(4,)),()),
             Grounded(Half(7,(1,1,2,3,4)),()), Grounded(Half(11,(1,2,3)),(4,)))
    for family in small:
        original_charge = sum(a%2 for a in member_word(family,0))
        for k in range(4):
            source, word = materialize(family,k,8192)
            check(recognize(family,source) == k, 'canonical recursive address')
            check(sum(a%2 for a in word) == original_charge, 'parity-preserving inflation')
            for modulus in (1,2,3,9,27,64,135,1024,2187):
                check(source_mod(family,k,modulus) == source % modulus, 'full denominator reader')
            if k:
                _, _, _, _, _, period = audit_family(family)
                prefix = family.receipt.word[:-1]+(family.receipt.word[-1]+period*k,)
                shifted = Grounded(Half(source,prefix),family.child_word)
                next_source, next_word = materialize(shifted,1,8192)
                direct_source, direct_word = materialize(family,k+1,8192)
                check((next_source,next_word) == (direct_source,direct_word), 'exact additive address law')
            for neighbor in (source-2,source+2):
                if neighbor > 0 and recognize(family,neighbor) is not None:
                    index = recognize(family,neighbor)
                    check(materialize(family,index,8192)[0] == neighbor, 'membership never guesses adjacent source')
    extension = extend_family(small[1],13,(3,))
    for k in range(4):
        n, _ = materialize(extension,k)
        child, _ = materialize(small[1],3*k)
        check(old.replay(n,(3,)) == child < n, 'index grows while actual source decreases')

    rows = []
    families = hostile_families()
    for family in families:
        h,p,q,b,depth,period = audit_family(family)
        n = family.receipt.source
        source, root_word = materialize(family,0)
        check(source == n, 'original point remains exact')
        # The k=1 source has astronomical bit length; all following are symbolic.
        for index in (0,1,2,17,2**64+1):
            word = member_word(family,index)
            check(len(word) == len(root_word) and sum(a%2 for a in word) == sum(a%2 for a in root_word),
                  'all-index exact rank and odd-charge invariants')
            for modulus in (3,9,19,32,729,1024,3**20,5**9):
                r = source_mod(family,index,modulus)
                # Independent modular division by splitting source into seed plus increment.
                lifted_modulus = p*modulus
                increment = q*h*((pow(2,period*index,lifted_modulus)-1) % lifted_modulus)
                check(increment % p == 0, 'independent increment integral')
                check(r == (n+increment//p) % modulus, 'independent symbolic source reader')
            check(source_mod(family,index,q) == n % q, 'same coarse native cylinder')
        reject(lambda family=family: materialize(family,1,8192))
        rows.append((n,h,len(family.receipt.word),sum(family.receipt.word),depth,
                     len(root_word),sum(a%2 for a in root_word),
                     source_mod(family,1,3**10),source_mod(family,1,2**20)))
    check(old.replay(SHARED_CHILD,SHARED_WORD) == 1 and (len(SHARED_WORD),sum(SHARED_WORD)) == (80,151),
          'new finite shared-child ROOT certificate')
    check(len(rows) == 11, 'eleven grounded source factories')
    first = next(f for f in families if f.receipt.source == 96883183)
    second = next(f for f in families if f.receipt.source == 129177577)
    check(extend_family(first,129177577,(2,)) == second, 'literal shared-child recursive factory')
    for k in (0,1,2,17,2**64+1):
        for modulus in (3,9,19,128,729,5**9):
            nmod = source_mod(second,k,4*modulus)
            check((3*nmod+1) % 4 == 0, 'exact prefix division precision')
            check(((3*nmod+1)//4) % modulus == source_mod(first,3*k,modulus),
                  'all-height recursive address bridge modular reader')
    g = small[0]
    for fn in (lambda: audit_half(replace(first_half,source=True)),
               lambda: audit_half(replace(first_half,word=list(first_half.word))),
               lambda: audit_half(Half(1,(2,))),lambda: audit_half(Half(9,(2,1,1,2,3))),
               lambda: lift(first_half,True),lambda: lift(first_half,-1),
               lambda: audit_family(Grounded(g.receipt,(2,))),
               lambda: source_mod(g,1,3.0),lambda: source_mod(g,True,3),
               lambda: materialize(g,1,-1),lambda: recognize(g,3.0),
               lambda: surplus(3,(1,4),True),lambda: splice(first_half,(1,)),
               lambda: extend_family(small[1],3,(1,)),lambda: extend_family(small[1],13,[3]),
               lambda: extend_family(small[1],13.0,(3,))):
        reject(fn)

    print('Native-tail recursion: every strict coarse lift of an actual half exit pays coefficient1')
    print('Declared finite control universe:511 odd sources3..1023;7 parameters each;cells=',cell_count,sep='')
    print('Original-source hostile56924955: deficit4; continuation(2,2) surplus141; total surplus137')
    print('Small supplied-terminal factories:4; indices0..3, literal ROOT and recursive address controls')
    for row in rows:
        print('seed,child,prefix_rank,cost,ternary_depth,ROOT_rank,odd_charge,next_mod59049,next_mod1048576:',row)
    print('Shared child16504781 frozen ROOT rank80,cost151,charge53,SHA256='+hashlib.sha256(bytes(SHARED_WORD)).hexdigest())
    print('Recursive bridge: B_k --(2)--> A_(3k), seeds129177577/96883183; actual source strictly decreases')
    print('Eleven infinite grounded subfamilies; k>0 local coefficient1; all-index ROOT rank and odd charge unchanged')
    print('No symbolic family is claimed to contain every native parameter, or to cover every positive odd source')
    print('checks=',checks,sep='')


if __name__ == '__main__':
    main()
