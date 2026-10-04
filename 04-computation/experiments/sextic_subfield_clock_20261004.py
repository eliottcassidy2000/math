"""Exact F64 subfield norms, tensor ranks, and retained 63-clock memory.

Independent bit-polynomial arithmetic is cross-checked against the inherited
sixth-clock power-basis implementation. Every acceptance check survives -O.
"""
from collections import Counter, defaultdict
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import gcd
from pathlib import Path
import json


def need(ok,message):
    if not ok:
        raise ValueError(message)


MODULUS = (1<<6)|(1<<3)|1
T,B = 2,3


def raw_mul(x,y):
    out = 0
    while y:
        if y&1:
            out ^= x
        y >>= 1
        x <<= 1
        if x&64:
            x ^= MODULUS
    return out


MULTIPLICATION = tuple(tuple(raw_mul(x,y) for y in range(64)) for x in range(64))


def mul(x,y):
    return MULTIPLICATION[x][y]


def power(x,n):
    if n < 0:
        need(x != 0, "inverse requires nonzero element")
        return power(x,n%63)
    out = 1
    while n:
        if n&1:
            out = mul(out,x)
        x = mul(x,x)
        n >>= 1
    return out


PHI,GAMMA,RHO,DELTA = power(B,21),power(B,9),power(B,28),power(B,36)
LOG_PHI = {power(PHI,i):i for i in range(3)}
LOG_GAMMA = {power(GAMMA,i):i for i in range(7)}
TENSOR_BASIS = tuple(mul(a,power(GAMMA,j)) for a in (1,PHI) for j in range(3))


def from_tensor(mask):
    out = 0
    for i,basis in enumerate(TENSOR_BASIS):
        if (mask>>i)&1:
            out ^= basis
    return out


TENSOR_COORDINATES = {from_tensor(mask):mask for mask in range(64)}


def tensor_rank(x):
    mask = TENSOR_COORDINATES[x]
    first,second = mask&7,mask>>3
    if not (first or second):
        return 0
    return 1 if first == 0 or second == 0 or first == second else 2


def norms(x):
    return power(x,21),power(x,9)


def encode_clock(x):
    """Only small subfield logarithms are used; retain the missing ternary digit."""
    need(x != 0, "nonzero 63-clock state")
    a,c = norms(x)
    i,ell = LOG_PHI[a],LOG_GAMMA[c]
    reference = mul(power(RHO,i),power(DELTA,ell))
    residual = mul(x,power(reference,-1))
    need(residual in LOG_PHI, "joint norms leave exactly an F4 multiplicative fibre")
    return i,LOG_PHI[residual],ell


def decode_clock(code):
    i,j,ell = code
    return mul(power(RHO,i+3*j),power(DELTA,ell))


def multiply_codes(a,b):
    i,j,ell = a
    k,h,m = b
    return (i+k)%3,(j+h+(i+k)//3)%3,(ell+m)%7


def step_code(code):
    i,j,ell = code
    return (i+1)%3,(j+int(i == 2))%3,(ell+1)%7


def trace(x,degree=6):
    out,current = 0,x
    for _ in range(degree):
        out ^= current
        current = mul(current,current)
    return out


def binary_rank(rows):
    rows = [sum(bit<<j for j,bit in enumerate(row)) for row in rows]
    rank = 0
    for bit in range(5,-1,-1):
        pivot = next((i for i in range(rank,len(rows)) if (rows[i]>>bit)&1),None)
        if pivot is None:
            continue
        rows[rank],rows[pivot] = rows[pivot],rows[rank]
        for i in range(len(rows)):
            if i != rank and ((rows[i]>>bit)&1):
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def polynomial(x):
    terms = ["1" if i == 0 else "t" if i == 1 else f"t^{i}" for i in range(6) if (x>>i)&1]
    return "+".join(reversed(terms)) or "0"


def describe(x):
    return dict(bits=x,polynomial=polynomial(x))


def exact_order(x):
    need(x != 0, "zero has no multiplicative order")
    return next(k for k in (1,3,7,9,21,63) if power(x,k) == 1)


def inherited_multiplication_audit():
    path = Path(__file__).with_name("sixth_clock_branches_20261004.py")
    spec = spec_from_file_location("inherited_sixth_clock",path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    for x,y in product(range(64),repeat=2):
        a,b = tuple((x>>i)&1 for i in range(6)),tuple((y>>i)&1 for i in range(6))
        expected = sum(bit<<i for i,bit in enumerate(module.smul(a,b,2)))
        need(mul(x,y) == expected, "bit-polynomial arithmetic equals inherited power-basis multiplication")
    return 4096


def main():
    independent_products = inherited_multiplication_audit()
    orbit = tuple(power(B,k) for k in range(63))
    need(set(orbit) == set(range(1,64)), "b generates all 63 nonzero field states")
    need(exact_order(T) == 9 and power(T,3) == PHI^1 and power(T,6) == PHI,
         "the retained sextic coordinate is an order-nine sixth root of phi")
    need(PHI == 9 and GAMMA == 38 and mul(PHI,PHI) == PHI^1, "inherited subfield generators")
    f4 = {x for x in range(64) if power(x,4) == x}
    f8 = {x for x in range(64) if power(x,8) == x}
    need(len(f4) == 4 and len(f8) == 8 and f4&f8 == {0,1}, "subfields and prime-field intersection")
    need(f4 == {0,*LOG_PHI} and f8 == {0,*LOG_GAMMA}, "small-generator subfield descriptions")
    for subfield in (f4,f8):
        for x,y in product(subfield,repeat=2):
            need(x^y in subfield and mul(x,y) in subfield, "subfield operations close")

    need(len(TENSOR_COORDINATES) == 64, "six tensor products form a binary basis")
    simple_products = {mul(x,y) for x,y in product(f4,f8)}
    cubes = {power(x,3) for x in range(1,64)}
    need(simple_products == cubes|{0} and len(simple_products) == 22, "pure tensors give exactly zero plus 21 cubes")
    rank_counts = Counter(tensor_rank(x) for x in range(64))
    need(rank_counts == {0:1,1:21,2:42}, "binary two-by-three tensor-rank census")
    for x in range(1,64):
        need((tensor_rank(x) == 1) == (x in cubes) == (power(x,21) == 1), "rank one iff cube iff trivial F4 norm")
    cosets = [{mul(power(B,j),x) for x in cubes} for j in range(3)]
    need(set.union(*cosets) == set(range(1,64)) and sum(map(len,cosets)) == 63,
         "three disjoint linear images of the 21-point pure-tensor surface")
    need(exact_order(T) == 9 and tensor_rank(T) == 2, "rank-two does not imply primitive")
    order_counts = Counter(exact_order(x) for x in range(1,64))
    need(order_counts == {1:1,3:2,7:6,9:6,21:12,63:36}, "complete multiplicative order census")

    fibres = defaultdict(list)
    for x in range(1,64):
        pair = norms(x)
        fibres[pair].append(x)
        need(pair[0] in f4-{0} and pair[1] in f8-{0}, "norms land in the named subfields")
        need(power(x,3) == mul(pair[0],power(pair[1],-2)), "joint norms recover the cube exactly")
        need((exact_order(x) == 63) == (pair[0] != 1 and pair[1] != 1), "joint norm test retains primitivity")
    need(set(fibres) == set(product(f4-{0},f8-{0})) and {len(v) for v in fibres.values()} == {3},
         "joint norm is onto the 21-state product with three-element fibres")
    need(set(fibres[(1,1)]) == f4-{0}, "joint kernel is exactly F4 star")
    need(power(1^PHI,21) != power(1,21)^power(PHI,21), "norm is not an additive coordinate projection")
    order_three_target_lifts = fibres[(PHI,1)]
    need(all(exact_order(x) == 9 for x in order_three_target_lifts), "no multiplicative section of the norm map")
    sixth_fibres = defaultdict(list)
    for x in range(1,64):
        sixth_fibres[power(x,6)].append(x)
    need(set(sixth_fibres) == cubes and {len(v) for v in sixth_fibres.values()} == {3},
         "sixth-power extraction has the same threefold loss")
    need(set(sixth_fibres[PHI]) == {mul(T,power(PHI,j)) for j in range(3)}, "all three sixth roots of phi")

    need(RHO == power(T,5) and DELTA == mul(B,power(T,4)) and mul(RHO,DELTA) == B,
         "root-depth normalization agrees with the order-nine/order-seven factorization")
    need(exact_order(RHO) == 9 and exact_order(DELTA) == 7 and power(RHO,3) == PHI,
         "CRT phase generators")
    codes = set(product(range(3),range(3),range(7)))
    need({encode_clock(x) for x in range(1,64)} == codes, "all 63 memory codes occur")
    for k,x in enumerate(orbit):
        code = encode_clock(x)
        need(code == (k%3,(k%9)//3,k%7) and decode_clock(code) == x, "small-subfield-log decoder is exact")
        need(encode_clock(mul(B,x)) == step_code(code), "one clock step uses the retained ternary carry")
        after21 = encode_clock(mul(power(B,21),x))
        need(after21 == (code[0],(code[1]+1)%3,code[2]), "visible 21-clock return advances the hidden third sheet")
    for x,y in product(range(1,64),repeat=2):
        need(encode_clock(mul(x,y)) == multiply_codes(encode_clock(x),encode_clock(y)),
             "all 3969 products respect the ternary carry law")

    absolute_trace = tuple(trace(x) for x in range(64))
    need(set(absolute_trace) == {0,1}, "absolute trace has binary values")
    for x,y in product(range(64),repeat=2):
        need(absolute_trace[x^y] == absolute_trace[x]^absolute_trace[y], "trace is binary linear")
    walsh = tuple(tuple(1-2*absolute_trace[mul(y,x)] for x in range(64)) for y in range(64))
    for y,z in product(range(64),repeat=2):
        need(sum(a*b for a,b in zip(walsh[y],walsh[z])) == 64*int(y == z), "all trace-Walsh rows are orthogonal")
        need(walsh[y][mul(B,z)] == walsh[mul(B,y)][z], "clock action and dual-character permutation intertwine")
    basis4,basis8 = (1,PHI),(1,GAMMA,power(GAMMA,2))
    gram4 = [[trace(mul(x,y),2) for y in basis4] for x in basis4]
    gram8 = [[trace(mul(x,y),3) for y in basis8] for x in basis8]
    gram64 = [[absolute_trace[mul(x,y)] for y in TENSOR_BASIS] for x in TENSOR_BASIS]
    for i,j,k,ell in product(range(2),range(3),range(2),range(3)):
        need(gram64[3*i+j][3*k+ell] == gram4[i][k]*gram8[j][ell], "tensor trace pairing factors exactly")
    need(binary_rank(gram64) == 6, "tensor trace pairing is nondegenerate")
    trace_word = tuple(absolute_trace[x] for x in orbit)
    windows = {sum(trace_word[(k+i)%63]<<i for i in range(6)) for k in range(63)}
    need(windows == set(range(1,64)), "six consecutive trace bits recover every nonzero state")
    need(Counter(trace_word) == {0:31,1:32}, "one output bit is not a complete state")
    need(all(trace_word[(k+6)%63] == (trace_word[(k+4)%63]^trace_word[(k+3)%63]^trace_word[(k+1)%63]^trace_word[k])
             for k in range(63)), "inherited primitive LFSR recurrence for trace output")

    report = dict(status="PROVED elementary F64 maps; FINITE-EXACT full-field controls; no Collatz conjugacy",
        field_polynomial="t^6+t^3+1",independent_inherited_product_checks=independent_products,
        generators={name:describe(x) for name,x in [('t',T),('b',B),('phi',PHI),('gamma',GAMMA),('rho',RHO),('delta',DELTA)]},
        subfields=dict(F4=[describe(x) for x in sorted(f4)],F8=[describe(x) for x in sorted(f8)],intersection=[0,1]),
        tensor_basis=[describe(x) for x in TENSOR_BASIS],tensor_rank_counts=dict(sorted(rank_counts.items())),
        multiplicative_order_counts=dict(sorted(order_counts.items())),
        tensor_rank_by_order={str(order):dict(sorted(Counter(tensor_rank(x) for x in range(1,64) if exact_order(x)==order).items())) for order in sorted(order_counts)},
        norm_map=dict(output_count=len(fibres),fibre_size=3,kernel=[describe(x) for x in sorted(fibres[(1,1)])],
            cube_decoder="x^3=N_F4(x)*N_F8(x)^(-2)",sixth_roots_of_phi=[describe(x) for x in sorted(sixth_fibres[PHI])],
            no_section_hostile=dict(target=(PHI,1),lifts=[describe(x) for x in order_three_target_lifts],lift_orders=[9,9,9])),
        clock_code=dict(coordinates="(i,j,l): k mod9=i+3j and k mod7=l",step="(i+1 mod3,j+[i=2] mod3,l+1 mod7)",
            decoder="rho^(i+3j)*delta^l",codes=len(codes),product_checks=3969,
            first_states=[dict(k=k,state=describe(orbit[k]),code=encode_clock(orbit[k])) for k in range(12)]),
        trace_fourier=dict(character_count=64,nontrivial_character_clock_length=63,gram_F4=gram4,gram_F8=gram8,
            tensor_gram=gram64,trace_clock_word="".join(map(str,trace_word)),six_bit_windows=len(windows),
            one_bit_counts=dict(sorted(Counter(trace_word).items()))),
        boundaries=["F4 tensor F8 has dimension6; the additive Cartesian product has dimension5.",
                    "Pure tensor rank2 alone does not imply primitivity; t has order9.",
                    "The norm pair is nonlinear and loses one ternary phase digit.",
                    "Binary vector addition, multiplicative clock time, and Collatz iteration are different operations."])
    print(json.dumps(report,indent=2))
    print("PASS: all checks remain active under -O")


if __name__ == '__main__':
    main()
