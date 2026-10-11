"""Exact controls for a real projective boundary with an arithmetic fibre.

No compactification theorem or universal convergence assertion is inferred
from a finite sample. The companion note proves the infinite statements.
"""
from fractions import Fraction as F
from itertools import product
import simplex_sphere_cube_geometry_20261009 as geometry
from collatz_debt_observer_integration_20261009 import replay_u

CHECKS = 0


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def integer(n, minimum=0):
    if type(n) is not int or n < minimum:
        raise ValueError('exact integer in stated domain required')
    return n


def terras(n):
    integer(n, 1)
    return (3*n+1)//2 if n % 2 else n//2


def source_chart(n):
    integer(n, 1)
    return (F(n,2*n+1),F(n,2*n+1),F(1,2*n+1))


def valuation(n):
    integer(n, 1)
    if n % 2 == 0:
        raise ValueError('odd input required')
    x=3*n+1
    return (x & -x).bit_length()-1


def valuation_progression(a):
    integer(a, 1)
    modulus=1 << (a+1)
    residue=((1 << a)-1)*pow(3,-1,modulus) % modulus
    return residue,modulus


def lift_residue(residue, depth, excess):
    integer(depth, 1); integer(excess, 1); integer(residue)
    modulus=1 << depth
    if residue >= modulus or residue % 2 == 0:
        raise ValueError('odd canonical residue required')
    # Grows without changing any of the depth stored arithmetic digits.
    return residue+modulus*(1 << excess)


def branch_matrix(parity):
    if type(parity) is not int or parity not in (0,1):
        raise ValueError('literal parity required')
    return ((1+2*parity,0,parity),(0,2,0),(0,0,2))


def recover_slope(matrix):
    return F(matrix[0][0],matrix[1][1])


def main():
    boundary=(F(1,2),F(1,2),F(0))
    for depth in range(1,10):
        for residue in range(1,1 << depth,2):
            for excess in (1,depth+2,2*depth+5):
                n=lift_residue(residue,depth,excess)
                point=source_chart(n)
                check(n % (1 << depth)==residue,'arithmetic address retained')
                check(max(abs(x-y) for x,y in zip(point,boundary))==F(1,2*n+1),
                      'exact real distance to the common boundary')
                check(point[0]/point[2]==n,'positive marker restores source')

    for a in range(1,17):
        residue,modulus=valuation_progression(a)
        previous=None
        for k in range(1,17):
            n=residue+modulus*(1 << k)
            check(valuation(n)==a,'every prescribed valuation approaches boundary')
            marker=source_chart(n)[2]
            if previous is not None:
                check(marker<previous,'marker tends down while valuation stays fixed')
            previous=marker

    for current,original,kind in product(range(1,45),range(1,24),('simplex','sphere','cube')):
        point=geometry.encode((current,original,1),kind)
        image=geometry.act(branch_matrix(current % 2),point)
        check(geometry.decode_monitor(image)==(terras(current),original),
              'arithmetic driver selects the actual homogeneous branch')
    for depth in range(1,11):
        for x in range(1,65):
            y=x+(1 << (depth+1))
            check((terras(y)-terras(x)) % (1 << depth)==0,
                  'one output digit of precision consumed')

    u,v=(1,3),(3,1)
    a,b=geometry.matrix_word(u),geometry.matrix_word(v)
    check(a[0][2]==5 and b[0][2]==11,'distinct ordered carries')
    check(recover_slope(a)==recover_slope(b)==F(9,16),'same boundary multiplier')
    for kind in ('simplex','sphere','cube'):
        for x,n in product(range(1,15),range(1,15)):
            infinity=geometry.encode((x,n,0),kind)
            check(geometry.act(a,infinity)==geometry.act(b,infinity),
                  'boundary forgets carry')
            finite=geometry.encode((x,n,1),kind)
            check(geometry.act(a,finite)!=geometry.act(b,finite),
                  'interior still separates carry')
    native=[]
    for word in (u,v):
        legal=[]
        for n in range(129,193,2):
            try:
                replay_u(n,word)
                legal.append(n)
            except ValueError:
                pass
        check(len(legal)==2,'two full native periods')
        native.append(tuple(n % 32 for n in legal))
    check(set(native[0]).isdisjoint(native[1]),'same slope, different native cells')

    for e in range(2,65):
        n=(1 << e)-1
        x=n
        for j in range(e):
            check(x==3**j*(1 << (e-j))-1 and x % 2==1,
                  'Mersenne shadow of the negative fixed point')
            if j>0:
                check(x>n,'no early Terras descent')
            x=terras(x)
        check(x==3**e-1 and x % 2==0,'exact end of the all-odd Terras prefix')

    slope=F(9,16)
    J=slope+1/slope
    check(slope*slope+1/(slope*slope)==J*J-2,
          'Chebyshev doubling describes the multiplier shadow')
    previous=None
    for k in range(1,65):
        original=4**k-1
        current=6**k-1
        check(current % (1 << k)==(1 << k)-1,
              'current arithmetic address approaches minus one')
        check(original % (1 << (2*k))==(1 << (2*k))-1,
              'original arithmetic address approaches minus one')
        check(F(current,original)>F(3,2)**k,
              'actual midpoint monitor approaches the current axis')
        # Both other normalized coordinates vanish as current/original grows.
        exact_error=F(original+1,current+original+1)
        if previous is not None:
            check(exact_error<previous,'approach to the fixed boundary is monotone')
        previous=exact_error
    for kind in ('simplex','sphere','cube'):
        fixed=geometry.encode((1,0,0),kind)
        check(geometry.act(branch_matrix(1),fixed)==fixed,
              'current-axis boundary is fixed under the odd branch')
    print('PROVED/FINITE-EXACT: real boundary collapse and retained arithmetic fibre')
    print('Odd residues through depth9; three unbounded lifts each; valuation1..16 progressions.')
    print('Marked branch transport: currents1..44, original1..23, all three charts.')
    print('Same slope9/16, carries5/11, native residues modulo32:',native)
    print('Mersenne exponents2..64: first E Terras parity bits are1; no earlier positive descent.')
    print('Mersenne midpoint monitors k1..64 approach the fixed boundary with two minus-one addresses.')
    print('The boundary kills carry; adding the 2-adic driver keeps branch dynamics continuous.')
    print('No source-specific ROOT proof or universal coverage theorem inferred.')
    print('PASS:',CHECKS,'exact checks')


if __name__=='__main__':
    main()
