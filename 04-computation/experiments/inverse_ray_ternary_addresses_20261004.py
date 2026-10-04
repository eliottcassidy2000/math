"""Proof-carrying inverse rays with row, block-height, and ternary addresses.

Exact integers only. Normal and -O execute identical checks. Large symbolic
certificates are queried modulo prime powers without expanding their sources.
The companion note separates this codec from global Collatz coverage.
"""
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import re


def need(condition, message):
    if not condition:
        raise ValueError(message)


def valuation(n, prime):
    need(n != 0 and prime in (2, 3), 'nonzero valuation argument')
    n = abs(n)
    count = 0
    while n % prime == 0:
        n //= prime
        count += 1
    return count


def odd_step(n):
    need(type(n) is int and n > 0 and n % 2 == 1, 'positive odd source')
    z = 3*n+1
    k = valuation(z, 2)
    return z >> k, k


def kappa(target_mod9, row):
    need(type(row) is int and row in (0, 1, 2), 'source row modulo three')
    need(target_mod9 % 3 != 0, 'inverse target must be a ternary unit')
    options = [k for k in range(1, 7)
               if pow(2, k, 9)*(target_mod9 % 9) % 9 == 1+3*row]
    need(len(options) == 1, 'one exponent in a six-block')
    return options[0]


def ray_source(target, row, block):
    need(type(target) is int and target > 0 and target % 2 == 1 and target % 3 != 0,
         'positive odd certified-hub candidate')
    need(type(block) is int and block >= 0, 'nonnegative integer block height')
    k = kappa(target % 9, row)+6*block
    numerator = (target << k)-1
    need(numerator % 3 == 0, 'inverse integrality')
    return numerator//3, k


def ray_residue(target, row, block, depth):
    need(type(depth) is int and depth >= 1, 'positive ternary precision')
    need(type(block) is int and block >= 0, 'nonnegative block height')
    k = kappa(target % 9, row)+6*block
    modulus = 3**(depth+1)
    numerator = (pow(2, k, modulus)*(target % modulus)-1) % modulus
    need(numerator % 3 == 0, 'modular inverse integrality')
    return numerator//3


def address_block(target, residue, depth):
    """One carry subtraction per digit; no discrete-log search over a fibre."""
    need(type(depth) is int and depth >= 1, 'positive ternary precision')
    need(type(target) is int and target > 0 and target % 2 and target % 3,
         'positive odd ternary-unit target')
    need(type(residue) is int and 0 <= residue < 3**depth, 'canonical ternary residue')
    row, block = residue % 3, 0
    for a in range(1, depth):
        current = ray_residue(target, row, block, a+1)
        need((residue-current) % 3**a == 0, 'retained lower address digits')
        digit = ((residue-current)//3**a) % 3
        block += digit*3**(a-1)
    need(block < 3**(depth-1), 'canonical block address')
    need(ray_residue(target, row, block, depth) == residue, 'requested source address')
    return row, block


def address_family(target, residue, depth):
    row, block = address_block(target, residue, depth)
    period = 3**(depth-1)
    first_parameter = int(target == 1 and row == 1 and block == 0)
    return row, block, period, first_parameter


@dataclass(frozen=True)
class Certificate:
    parent: object = None
    row: int = -1
    block: int = 0


ROOT = Certificate()


@lru_cache(None)
def mod3(cert, depth):
    need(type(depth) is int and depth >= 0, 'nonnegative ternary precision')
    if depth == 0:
        return 0
    if cert.parent is None:
        need(cert == ROOT, 'the unique root certificate')
        return 1
    parent_residue = mod3(cert.parent, depth+1)
    k = kappa(parent_residue % 9, cert.row)+6*cert.block
    modulus = 3**(depth+1)
    numerator = (pow(2, k, modulus)*parent_residue-1) % modulus
    need(numerator % 3 == 0, 'symbolic ternary division')
    return numerator//3


@lru_cache(None)
def exponent(cert):
    need(cert.parent is not None, 'root has no incoming inverse edge')
    return kappa(mod3(cert.parent, 2), cert.row)+6*cert.block


@lru_cache(None)
def mod2(cert, depth):
    need(type(depth) is int and depth >= 0, 'nonnegative binary precision')
    if depth == 0:
        return 0
    if cert.parent is None:
        need(cert == ROOT, 'the unique root certificate')
        return 1
    modulus, k = 1 << depth, exponent(cert)
    factor = 0 if k >= depth else 1 << k
    return (factor*mod2(cert.parent, depth)-1)*pow(3, -1, modulus) % modulus


def extend(parent, row, block):
    audit_certificate(parent)
    need(type(block) is int and block >= 0, 'nonnegative block height')
    need(type(row) is int and row in (0, 1, 2), 'source row')
    need(mod3(parent, 1) != 0, 'a row-zero leaf has no odd predecessor')
    need(not (parent == ROOT and row == 1 and block == 0), 'exclude root self-return')
    child = Certificate(parent, row, block)
    need(mod3(child, 1) == row and mod2(child, 1) == 1, 'typed odd child')
    return child


def chain(cert):
    output = []
    need(isinstance(cert, Certificate), 'typed certificate root or node')
    while cert.parent is not None:
        output.append(cert)
        cert = cert.parent
        need(isinstance(cert, Certificate), 'typed parent pointer')
    need(type(cert.row) is int and type(cert.block) is int and cert == ROOT,
         'retained typed terminal witness')
    return output


def audit_certificate(cert):
    """Validate supplied structural data without expanding any source integer."""
    for node in reversed(chain(cert)):
        need(type(node.row) is int and node.row in (0, 1, 2), 'stored source row')
        need(type(node.block) is int and node.block >= 0, 'stored nonnegative block height')
        need(mod3(node.parent, 1) != 0, 'stored parent cannot be a ternary leaf')
        need(not (node.parent == ROOT and node.row == 1 and node.block == 0),
             'stored root self-return is not first-hit data')
    return True


def ranks(cert):
    nodes = chain(cert)
    return len(nodes), sum(exponent(node)+1 for node in nodes)


def bit_bounds(cert):
    lower = upper = 1
    for node in reversed(chain(cert)):
        k = exponent(node)
        lower, upper = max(1, k+lower-2), k+upper-1
    return lower, upper


def expand(cert, bit_cap=10000):
    need(bit_bounds(cert)[1] <= bit_cap, 'declared expansion bit cap')
    value = 1
    for node in reversed(chain(cert)):
        k = exponent(node)
        numerator = (value << k)-1
        need(numerator % 3 == 0, 'expanded integer division')
        value = numerator//3
    return value


def encode_source(source, step_cap=10000):
    """Canonical encoder on a supplied finite first-hit route; may hit its cap."""
    need(type(source) is int and source > 0 and source % 2 == 1, 'positive odd integer source')
    need(type(step_cap) is int and step_cap >= 0, 'nonnegative integer forward cap')
    current, trail = source, []
    while current != 1:
        need(len(trail) < step_cap, 'declared forward certificate search cap')
        target, k = odd_step(current)
        trail.append((current, target, k))
        current = target
    cert = ROOT
    for source, target, k in reversed(trail):
        row = source % 3
        least = kappa(target % 9, row)
        need(k >= least and (k-least) % 6 == 0, 'canonical row/block decoder')
        cert = extend(cert, row, (k-least)//6)
    return cert


def literal_certificate_check(cert):
    audit_certificate(cert)
    source = expand(cert)
    current = source
    ordinary = 0
    for node in chain(cert):
        need(current != 1, 'first arrival has not occurred early')
        target, k = odd_step(current)
        need(k == exponent(node) and current % 3 == node.row, 'independent literal edge')
        current, ordinary = target, ordinary+k+1
    need(current == 1 and ranks(cert) == (len(chain(cert)), ordinary), 'literal first-hit ranks')
    need(encode_source(source) == cert, 'canonical source/code round trip')
    need(bit_bounds(cert)[0] <= source.bit_length() <= bit_bounds(cert)[1], 'bit-length interval')
    for a in range(1, 8):
        need(mod3(cert, a) == source % 3**a, 'symbolic/expanded ternary projection')
    for a in range(1, 17):
        need(mod2(cert, a) == source % 2**a, 'symbolic/expanded binary projection')
    return source


def main():
    print('INVERSE RAY TERNARY ADDRESSES: PROVED codec; FINITE-EXACT declared controls')
    print('kappa table: target mod9 -> source rows0,1,2')
    for residue in (1, 2, 4, 5, 7, 8):
        print(residue, [kappa(residue, row) for row in range(3)])
    print('root nonself first sources by rows0,2,1:',
          [ray_source(1, 0, 0)[0], ray_source(1, 2, 0)[0], ray_source(1, 1, 1)[0]])

    count = 0
    for depth in range(1, 5):
        modulus, period = 3**depth, 3**(depth-1)
        for target in range(1, 2*3**(depth+1), 2):
            if target % 3 == 0:
                continue
            for row in range(3):
                # Independent orbit iteration from one modular source, not powers or logs.
                current = ray_residue(target, row, 0, depth)
                observed = []
                for block in range(period):
                    need(address_block(target, current, depth) == (row, block), 'complete clock-address inverse')
                    need(ray_residue(target, row, block, depth) == current, 'power/affine orbit agreement')
                    observed.append(current)
                    current = (64*current+21) % modulus
                    count += 1
                need(current == observed[0] and set(observed) == set(range(row, modulus, 3)),
                     'one full fixed-row clock')
                for block in (0, 1, 7):
                    for digit in range(3):
                        low = ray_residue(target, row, block, depth+1)
                        high = ray_residue(target, row, block+digit*period, depth+1)
                        need((high-low-digit*3**depth) % 3**(depth+1) == 0, 'unit carry lift')
    need(count == 14760, 'complete finite target/address universe')
    print('all target residues at ternary depths1..4: address/affine-orbit checks', count)

    for target in (1, 5, 17, 89):
        for depth in range(1, 9):
            for residue in (0, 1, 2, 3**depth-1):
                row, block, period, first = address_family(target, residue, depth)
                for t in (first, first+1, first+2):
                    source, k = ray_source(target, row, block+period*t)
                    need(source > 1 and source % 3**depth == residue, 'root-offset guarded address family')
                    need(odd_step(source) == (target, k), 'address family literal inverse edge')
    for block in range(300):
        source, _ = ray_source(1, 0, block)
        need(valuation(source, 3) == 1+valuation(block+1, 3), 'root row-zero depth law')
    print('four certified-hub address families: depths1..8, four residues, three lifts each; 300 leaf-depth controls')

    # Fully expanded independent universe of row/block certificates.
    all_certs, frontier, census = [ROOT], [ROOT], []
    for depth in range(1, 6):
        next_frontier, leaves = [], []
        for parent in frontier:
            for row in range(3):
                for block in (0, 1):
                    if parent == ROOT and row == 1 and block == 0:
                        continue
                    child = extend(parent, row, block)
                    all_certs.append(child)
                    (leaves if row == 0 else next_frontier).append(child)
        census.append((depth, len(next_frontier), len(leaves)))
        frontier = next_frontier
    need(len(all_certs) == 1536, 'complete bounded certificate universe')
    values = [literal_certificate_check(cert) for cert in all_certs]
    need(len(set(values)) == len(values), 'first-hit codes are injective across depths')
    print('depth/internal/leaf census', census)
    print('expanded certificates', len(values), 'largest source bits', max(n.bit_length() for n in values))

    # Hostiles: root loop, terminal leaf, omitted hub precision, residue collisions.
    for parent, row, block in ((ROOT,1,0), (extend(ROOT,0,0),1,0)):
        failed = False
        try:
            extend(parent,row,block)
        except ValueError:
            failed = True
        need(failed, 'expected structural rejection')
    need(11 % 9 == 29 % 9, 'same one-generation hub address')
    paths = []
    for target in (11,29):
        path = [target]
        for _ in range(2):
            a0 = 2 if target % 3 == 1 else 1
            target = ((target << a0)-1)//3
            path.append(target)
        paths.append(path)
    need(paths == [[11,7,9],[29,19,25]], 'one more generation consumes one more ternary digit')
    for depth in range(1, 9):
        row, block, period, _ = address_family(5, 1, depth)
        first, _ = ray_source(5,row,block)
        second, _ = ray_source(5,row,block+period)
        need(first != second and (first-second) % 3**depth == 0, 'address does not retain block-height quotient')
    print('precision hostile', paths, '; root self-return and inverse continuation from21 rejected')
    for invalid in (Certificate(None,-1.0,False), Certificate(ROOT,1,0), Certificate(ROOT,2,-1),
                    Certificate(Certificate(ROOT,0,0),2,0), Certificate(ROOT,3,0)):
        failed = False
        try:
            audit_certificate(invalid)
        except ValueError:
            failed = True
        need(failed, 'tampered symbolic certificate rejected')
    for invalid_source in (True,1.0,0,2,-1):
        failed = False
        try:
            encode_source(invalid_source)
        except ValueError:
            failed = True
        need(failed, 'noninteger or nonpositive-odd source rejected at root boundary')

    # Completed sources outside the repeated-growth atlas, without a bank-novelty claim.
    for j in range(51):
        source = ((11 << (1+2*j))-1)//3
        cert = encode_source(source)
        need(ranks(cert) == (5,16+2*j), 'entire hub11 fibre has fixed odd first-hit rank')
        need(odd_step(source) == (11,1+2*j), 'hub11 fibre exact first edge')
        if j == 0:
            current, prefix = source, []
            while current >= source:
                current,k = odd_step(current)
                prefix.append(k)
            need(prefix == [1,1,2,3], 'source7 first-descent obstruction')
        else:
            need(1+2*j > 1 and source > 11, 'immediate descent outside every favorable repeated chart')
    bank = Path(__file__).resolve().parents[2]/'05-knowledge/results/reset_20260926_swaplift.out'
    inherited = []
    for line in bank.read_text().splitlines():
        if re.fullmatch(r'\d+(?: \d+){8}',line):
            row = list(map(int,line.split()))
            if row[-2:] == [7,7]:
                inherited.append(row[0])
    need(inherited == [25,89,153,217,281], 'source7 cylinder already belongs to the named bank')
    for t in range(100):
        source,current = 7+128*t,7+128*t
        word = []
        for i in range(4):
            current,k = odd_step(current)
            word.append(k)
            need((current < source) == (i == 3), 'inherited7mod128 exact first descent')
        need(word[:3] == [1,1,2] and word[-1] >= 3, 'nonperiodic exit prefix')
    print('hub11 fibre j0..50: ranks(5,16+2j); source7mod128 already in old-bank rows', inherited)

    # A huge source is never materialized. Every node still carries a proved edge to its parent.
    huge = encode_source(5)
    huge_blocks = []
    for i in range(20):
        block = 10**100 + i*10**30
        huge_blocks.append(block)
        huge = extend(huge, 1+i%2, block)
    terminal_leaf = extend(huge, 0, 10**100+21)
    need(ranks(huge)[0] == 21 and ranks(terminal_leaf)[0] == 22, 'symbolic first-hit odd rank')
    for cert in (huge, terminal_leaf):
        audit_certificate(cert)
        for node in chain(cert):
            k = exponent(node)
            for a in (1, 2, 5, 12):
                modulus = 3**(a+1)
                need((3*mod3(node,a+1)+1-pow(2,k,modulus)*mod3(node.parent,a+1)) % modulus == 0,
                     'symbolic ternary edge identity')
            for a in (1, 7, 64):
                modulus = 1 << a
                need((3*mod2(node,a)+1-pow(2,k,modulus)*mod2(node.parent,a)) % modulus == 0,
                     'symbolic binary edge identity')
        need(bit_bounds(cert)[0] > 10**102, 'source necessarily exceeds expansion scale')
    print('huge internal certificate odd/ordinary ranks',ranks(huge))
    print('huge internal source binary-length interval',bit_bounds(huge))
    print('huge internal residues mod2^64 / mod3^12',mod2(huge,64),mod3(huge,12))
    print('adding one row-zero leaf gives odd rank',ranks(terminal_leaf)[0],
          '; no expanded source was created')
    print('Universe: all14760 row/address pairs at depths1..4; 1536 expanded codes; two huge symbolic codes.')
    print('ALL CHECKS PASSED; global integer coverage remains OPEN')


if __name__ == '__main__':
    main()
