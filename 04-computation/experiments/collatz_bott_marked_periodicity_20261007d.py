"""Marked stabilization and an exact child-phase monodromy.

No astronomical source is expanded. Native phase packets are not ROOT proofs.
The Bott background is not a dependency of these elementary computations.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as Q
from math import comb
from itertools import product

import collatz_eight_bit_completion_20261007c as prior
import collatz_twoanchor_head_decoder_20261007b as words


BASE = prior.PHASE-8
STRIDE = 1 << 32
PREFIX = prior.RIGHT
PRECISION_SHIFT = 34
OFFSET = BASE+14
CONSTANT = 8589800907


def natural(n, minimum=0):
    if type(n) is not int or n < minimum:
        raise ValueError("exact integer in declared domain required")
    return n


def j_mod(parameter, bits):
    """Normalized endpoint (Y-5)/8, with full division precision retained."""
    natural(parameter)
    natural(bits)
    modulus = 1 << (bits+PRECISION_SHIFT)
    numerator = (pow(3, OFFSET+STRIDE*parameter, modulus)-CONSTANT) % modulus
    if numerator % (1 << PRECISION_SHIFT):
        raise ArithmeticError("marked family lost endpoint integrality")
    return numerator >> PRECISION_SHIFT


def endpoint_mod(parameter, bits):
    natural(parameter)
    natural(bits)
    return (8*j_mod(parameter, max(0,bits-3))+5) % (1 << bits)


def monodromy(value, steps, bits):
    """Exact parameter advance on j residues; no source-orbit iteration."""
    natural(value)
    natural(steps)
    natural(bits)
    modulus = 1 << bits
    power = pow(3, STRIDE*steps, modulus << PRECISION_SHIFT)
    quotient = (power-1) >> PRECISION_SHIFT
    return (power*value+CONSTANT*quotient) % modulus


def newton_degree(bits):
    natural(bits)
    return (bits+PRECISION_SHIFT-1)//PRECISION_SHIFT


def polynomial_mod(parameter, bits):
    """Minimal-degree integer-valued Newton polynomial at the given precision."""
    natural(parameter)
    natural(bits)
    modulus = 1 << bits
    w = (pow(3,STRIDE,modulus << PRECISION_SHIFT)-1) >> PRECISION_SHIFT
    scale = pow(3,OFFSET,modulus)
    value = j_mod(0,bits)
    for i in range(1,newton_degree(bits)+1):
        value += scale*(1 << (PRECISION_SHIFT*(i-1)))*pow(w,i,modulus)*comb(parameter,i)
    return value % modulus


def inverse_j(value, bits):
    """The unique parameter class modulo2^bits realizing a j residue."""
    natural(value)
    natural(bits)
    value %= 1 << bits
    t = 0
    for k in range(1,bits+1):
        choices = [s for s in (t,t+(1 << (k-1))) if j_mod(s,k)==value%(1 << k)]
        if len(choices) != 1:
            raise ArithmeticError("the proved binary isometry lost a unique lift")
        t = choices[0]
    return t


@dataclass(frozen=True)
class Pullback:
    word: tuple
    residue: int
    period: int
    strict_parameter_cut: int


def pullback(word):
    """Formal native continuation iff; strict first-hit cutoff is separate."""
    words.word_type(word)
    if not word:
        return Pullback((),0,1,0)
    if word[0] < 3:
        return None
    p,q,b = words.carrier(word)
    native = ((q-b)*pow(p,-1,2*q)) % (2*q)
    if native % 8 != 5:
        raise ArithmeticError("continuation did not respect the retained endpoint disk")
    bits = sum(word)-2
    residue = inverse_j((native-5)//8,bits)
    cut = max(0,(2*(32+sum(word))+1-BASE+STRIDE-1)//STRIDE)
    return Pullback(word,residue,1 << bits,cut)


def audit(packet):
    if type(packet) is not Pullback:
        raise ValueError("exact Pullback required")
    for value in (packet.residue,packet.period,packet.strict_parameter_cut):
        natural(value)
    if pullback(packet.word) != packet:
        raise ValueError("forged word/phase/cut packet")
    return packet


def contains(packet, parameter, strict=False):
    audit(packet)
    natural(parameter)
    if type(strict) is not bool:
        raise ValueError("exact first-hit requirement flag")
    return ((parameter-packet.residue) % packet.period == 0
            and (not strict or parameter >= packet.strict_parameter_cut))


def continuation_native(word, parameter):
    """Independent modular check of the source word's native endpoint guard."""
    words.word_type(word)
    natural(parameter)
    p,q,b = words.carrier(word)
    return (p*endpoint_mod(parameter,sum(word)+1)+b) % (2*q) == q


def literal_residue_control(parameter, bits):
    """Replay a modest residue representative, not the astronomical Mersenne."""
    natural(parameter)
    natural(bits)
    precision = 32+bits+3
    modulus = 1 << precision
    x = (2*pow(3,BASE-1+STRIDE*parameter,modulus)-1) % modulus
    x += modulus  # every intermediate stays above1 through the fixed cost32.
    for a in PREFIX:
        z = 3*x+1
        if (z & -z).bit_length()-1 != a:
            raise ArithmeticError("fixed actual prefix rejected the residue representative")
        x = z >> a
    return ((x-5)//8) % (1 << bits)


def matrix_type(matrix):
    if (type(matrix) is not tuple or not matrix or any(type(row) is not tuple
            or len(row)!=len(matrix) for row in matrix)
            or any(type(a) not in (int,Q) for row in matrix for a in row)):
        raise ValueError("square exact rational tuple matrix required")
    return matrix


def stabilize(matrix, size=16):
    matrix_type(matrix)
    natural(size,1)
    return tuple(tuple(matrix[i//size][j//size] if i%size==j%size else Q(0)
                       for j in range(len(matrix)*size)) for i in range(len(matrix)*size))


def marked_corner(matrix, size=16):
    matrix_type(matrix)
    natural(size,1)
    if len(matrix)%size:
        raise ValueError("marked block dimensions required")
    small = tuple(tuple(matrix[i*size][j*size] for j in range(len(matrix)//size))
                  for i in range(len(matrix)//size))
    if stabilize(small,size) != matrix:
        raise ValueError("operator lost the marked scalar-block structure")
    return small


def multiply(a,b):
    matrix_type(a);matrix_type(b)
    if len(a)!=len(b):raise ValueError("matching dimensions")
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(a)))
                       for j in range(len(a))) for i in range(len(a)))


def word_matrix(word):
    words.word_type(word)
    p,q,b=words.carrier(word)
    return ((Q(p,q),Q(b,q)),(Q(0),Q(1)))


def main():
    checks=0
    def check(ok):
        nonlocal checks
        checks+=1
        if not ok:raise ArithmeticError("exact control failed")
    def rejects(f):
        try:f()
        except (ValueError,TypeError):check(True)
        else:check(False)

    p,q,b=words.carrier(PREFIX)
    check((len(PREFIX),sum(PREFIX))==(15,32))
    check(CONSTANT==(5*q+p-b)//2 and OFFSET==BASE+len(PREFIX)-1)
    check((pow(3,STRIDE,1 << 35)-1)%(1 << 35)==1 << 34)
    # Independent actual finite residue readers, finite permutation and unique lifting.
    for bits in range(0,11):
        images=[j_mod(t,bits) for t in range(1 << bits)]
        check(sorted(images)==list(range(1 << bits)))
        for t,value in enumerate(images):
            check(inverse_j(value,bits)==t)
            check(value==literal_residue_control(t,bits))
            check(monodromy(value,1,bits)==j_mod(t+1,bits))
    for bits in (1,8,33,34,35,47,68,69,103):
        for t in (0,1,2,3,17,1 << 47,10**30):
            check(polynomial_mod(t,bits)==j_mod(t,bits))
            check(monodromy(j_mod(t,bits),17,bits)==j_mod(t+17,bits))
            check(j_mod(t+(1 << bits),bits)==j_mod(t,bits))
        check(j_mod(0,bits+1)!=j_mod(1 << bits,bits+1))
    difference=(j_mod(2,35)-2*j_mod(1,35)+j_mod(0,35))%(1 << 35)
    check(difference==1 << 34)
    check(newton_degree(34)==1 and newton_degree(47)==2 and newton_degree(69)==3)
    # Every first valuation>=3 exists, and later words pull back uniquely.
    tested=0
    for ell in range(1,4):
        for word in product(range(1,6),repeat=ell):
            record=pullback(word)
            check((record is None)==(word[0]<3))
            if record is None:continue
            audit(record);tested+=1
            for parameter in (record.residue,record.residue+record.period,record.residue+1):
                check(contains(record,parameter)==continuation_native(word,parameter))
            if ell>1:
                parent=pullback(word[:-1])
                check((record.residue-parent.residue)%parent.period==0)
    check(pullback(())==Pullback((),0,1,0))
    # New actual rule, independently owned. No execution of its main/discovery.
    import collatz_eight_child_routes_20261007d as child
    check(child.LEFT[:15]==PREFIX)
    suffix=child.LEFT[15:]
    refined=pullback(suffix)
    check((len(suffix),sum(suffix))==(24,49))
    check((refined.residue,refined.period,refined.strict_parameter_cut)==(0,1 << 47,0))
    for t in (0,1,2,1 << 46,1 << 47,3 << 47):
        check(contains(refined,t)==child.contains(BASE+STRIDE*t))
    check(j_mod(0,46)==j_mod(1 << 46,46))
    check(not contains(refined,1 << 46) and contains(refined,0))

    # Marked stabilization is lossless; homotopy/similarity discards guard/carry.
    m12,m21=word_matrix((1,2)),word_matrix((2,1))
    translate=((Q(1),Q(-2)),(Q(0),Q(1)))
    inverse=((Q(1),Q(2)),(Q(0),Q(1)))
    check(multiply(multiply(translate,m12),inverse)==m21)
    check(words.native_cell((1,2))==(11,16) and words.native_cell((2,1))==(9,16))
    for w in ((),(1,),(2,),(1,2),(2,1),PREFIX):
        m=word_matrix(w)
        for size in (1,2,16):
            check(marked_corner(stabilize(m,size),size)==m)
        for s in (Q(0),Q(1,2),Q(1)):
            check(1-s+s*m[0][0]>0)
    check(words.carrier((1,))[0]>words.carrier((1,))[1])
    check(words.carrier((2,))[0]<words.carrier((2,))[1])
    check((3*3+1)//2==5 and (3*9+1)//4==7)
    a=stabilize(m12);broken=[list(row) for row in a];broken[1][0]=1
    rejects(lambda: marked_corner(tuple(tuple(row) for row in broken)))
    rejects(lambda: stabilize(((True,0),(0,1))))
    rejects(lambda: j_mod(True,10))
    rejects(lambda: inverse_j(3,4.0))
    rejects(lambda: pullback((3,False)))
    rejects(lambda: audit(replace(refined,period=True)))
    rejects(lambda: audit(replace(refined,residue=1)))
    rejects(lambda: contains(refined,0,strict=1))

    print('PROVED marked endpoint j(t)=[3^(924745911+2^32*t)-8589800907]/2^34.')
    print('PROVED v2(j(t)-j(s))=v2(t-s); every2^h quotient is one cycle of length2^h.')
    print('Newton reader minimal degree ceil(h/34); at47bits exactlyquadratic; linear hostile at35bits,t2.')
    print('Continuation words begin>=3; costA gives unique tphase modulo2^(A-2), formal guards with strict cutoff retained.')
    print('FINITE-EXACT continuation words:',tested,'; full quotient permutations through10bits.')
    print('Actual39letter/cost81 child rule leaves24letter/cost49 suffix: t=0mod2^47.')
    print('Same low46bit coordinate t0/t2^46 disagrees on newguard; retain refinement bits.')
    print('Marked M tensor I16 recovers M; conjugate words12/21 have disjoint native cells11/9mod16.')
    print('Bott/Clifford stabilization is cited background, not an8-step Collatz return or positivity theorem.')
    print('No astronomical Mersenne expansion, ROOT discovery, or universal coverage assertion.')
    print('Exact checks',checks)


if __name__=='__main__':main()
