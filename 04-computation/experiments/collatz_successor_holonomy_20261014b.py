"""Native Collatz block concatenation versus finite affine-group identity.

All production paths authenticate supplied sources.  Discovery is unnecessary.
The symbolic ROOT family is constructed explicitly, not asserted to cover inputs.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from math import gcd, lcm
import json


def integer(n, minimum=0):
    if type(n) is not int or n < minimum:
        raise ValueError("exact integer outside domain")
    return n


def odd(n):
    integer(n, 1)
    if n % 2 != 1:
        raise ValueError("positive odd integer required")
    return n


def word(w):
    if type(w) is not tuple:
        raise ValueError("exact word tuple required")
    for a in w:
        integer(a, 1)
    return w


def carrier(w):
    word(w)
    p, q, b = 1, 1, 0
    for a in w:
        p, q, b = 3*p, q*(1 << a), 3*b+q
    return p, q, b


def native_cell(w):
    p, q, b = carrier(w)
    return ((q-b)*pow(p, -1, 2*q)) % (2*q), 2*q


def replay(source, w):
    odd(source)
    word(w)
    n = source
    for a in w:
        if n == 1:
            raise ValueError("ROOT padding")
        z = 3*n+1
        actual = (z & -z).bit_length()-1
        if actual != a:
            raise ValueError("non-native source word")
        n = z >> a
    return n


@dataclass(frozen=True)
class Receipt:
    source: int
    target: int
    word: tuple


def receipt(source, w):
    return Receipt(source, replay(source, w), w)


def audit_receipt(r):
    if type(r) is not Receipt:
        raise ValueError("typed receipt required")
    odd(r.source)
    odd(r.target)
    word(r.word)
    if receipt(r.source, r.word) != r:
        raise ValueError("receipt does not authenticate")
    return r


def concatenate(left, right):
    audit_receipt(left)
    audit_receipt(right)
    if left.target != right.source:
        raise ValueError("integer seam mismatch")
    return receipt(left.source, left.word+right.word)


def discharge(r, supplied_root_word):
    audit_receipt(r)
    if replay(r.target, supplied_root_word) != 1:
        raise ValueError("supplied child word does not reach ROOT")
    return concatenate(r, receipt(r.target, supplied_root_word))


def unit_modulus(m):
    integer(m, 2)
    if gcd(m, 6) != 1:
        raise ValueError("affine group requires a modulus coprime to six")
    return m


def affine_mod(w, modulus):
    unit_modulus(modulus)
    p, q, b = carrier(w)
    inv = pow(q, -1, modulus)
    return p*inv % modulus, b*inv % modulus


def order(a, modulus):
    unit_modulus(modulus)
    integer(a)
    a %= modulus
    if gcd(a, modulus) != 1:
        raise ValueError("unit required")
    x = 1
    for k in range(1, modulus+1):
        x = x*a % modulus
        if x == 1:
            return k
    raise RuntimeError("finite unit group order bound failed")


def neutral_length(modulus):
    """Simultaneously neutralizes (3/2)^k and (9/16)^k."""
    unit_modulus(modulus)
    return lcm(order(3*pow(2,-1,modulus) % modulus, modulus),
               order(9*pow(16,-1,modulus) % modulus, modulus))


def paid_word(k):
    integer(k, 1)
    return (1,)*k + (2,)*(2*k)


def paid_data(k):
    integer(k, 1)
    p, q = 3**(3*k), 2**(5*k)
    b = q+p-2**(k+1)*3**(2*k)
    return p, q, b


def paid_source(k, index):
    integer(k, 1)
    integer(index)
    modulus = 2**(5*k+1)
    residue = (2**(k+1)*pow(3**k,-1,modulus)-1) % modulus
    return residue + modulus*index


def paid_receipt(source, k):
    odd(source)
    w = paid_word(k)
    r, modulus = native_cell(w)
    if source % modulus != r:
        raise ValueError("source outside the native paid cell")
    result = receipt(source, w)
    if result.target >= source:
        raise RuntimeError("proved paid-block inequality failed")
    return result


def recognize_paid(source):
    """Unique possible W_k address, determined from the unchanged source."""
    odd(source)
    k = ((source+1) & -(source+1)).bit_length()-2
    if k < 1:
        return None
    u = (source+1) >> (k+1)
    if (3**k*u-1) % 2**(4*k):
        return None
    return paid_receipt(source,k)


def principal_log4(target, depth):
    """Unique exponent modulo3^(depth-1), by three-candidate digit lifting."""
    integer(target)
    integer(depth, 1)
    modulus = 3**depth
    target %= modulus
    if target % 3 != 1:
        raise ValueError("principal 3-adic unit required")
    exponent, period = 0, 1
    for d in range(2, depth+1):
        mod = 3**d
        choices = [exponent+j*period for j in range(3)
                   if pow(4, exponent+j*period, mod) == target % mod]
        if len(choices) != 1:
            raise RuntimeError("principal-unit lift is not unique")
        exponent = choices[0]
        period *= 3
    return exponent, period


@dataclass(frozen=True)
class RootPlan:
    k: int
    index: int
    j: int
    period: int


def root_plan(k, index=0):
    integer(k, 1)
    integer(index)
    p, q, b = paid_data(k)
    denominator = 3*p
    target = (1+3*b*pow(q,-1,denominator)) % denominator
    phase, period = principal_log4(target, 3*k+1)
    # j>=2 suffices for positive source and excludes a padded final ROOT step.
    first = phase if phase >= 2 else phase+period
    return RootPlan(k, index, first+period*index, period)


def audit_plan(plan):
    if type(plan) is not RootPlan:
        raise ValueError("typed ROOT plan required")
    for n in (plan.k, plan.index, plan.j, plan.period):
        integer(n)
    if root_plan(plan.k, plan.index) != plan:
        raise ValueError("noncanonical ROOT phase")
    return plan


def source_residue(plan, modulus):
    audit_plan(plan)
    integer(modulus, 1)
    p, q, b = paid_data(plan.k)
    d = 3*p
    numerator = (q*(pow(4,plan.j,d*modulus)-1)-3*b) % (d*modulus)
    if numerator % d:
        raise RuntimeError("phase failed exact division")
    return numerator//d


def endpoint_residue(plan, modulus):
    audit_plan(plan)
    integer(modulus, 1)
    return ((pow(4,plan.j,3*modulus)-1)//3) % modulus


def materialize(plan, bit_cap=100000):
    audit_plan(plan)
    integer(bit_cap, 1)
    if 5*plan.k+2*plan.j > bit_cap:
        raise ValueError("explicit materialization cap; use modular readers")
    p, q, b = paid_data(plan.k)
    numerator = q*(4**plan.j-1)-3*b
    if numerator % (3*p):
        raise RuntimeError("ROOT phase failed")
    return numerator//(3*p)


def root_word_at_source(source, k):
    """Recognize this explicit family at the unchanged source; no orbit search."""
    odd(source)
    integer(k, 1)
    p, q, b = paid_data(k)
    value, rem = divmod(p*source+b, q)
    if rem or value <= 1:
        return None
    power = 3*value+1
    if power & (power-1):
        return None
    cost = power.bit_length()-1
    if cost < 4 or cost % 2:
        return None
    w = paid_word(k)+(cost,)
    try:
        return w if replay(source,w) == 1 else None
    except ValueError:
        return None


CHECKS = 0


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def rejects(call):
    try:
        call()
    except (ValueError, TypeError):
        check(True, "rejected")
        return
    check(False, "invalid input accepted")


def main():
    # Complete native residues for k1,k2, verified by independent literal steps.
    native_counts = []
    for k in (1,2):
        w = paid_word(k)
        residue, mod = native_cell(w)
        hits = []
        for n in range(1,mod,2):
            x, valid = n, True
            for a in w:
                if x == 1:
                    valid = False
                    break
                z = 3*x+1
                if z % 2**a or (z//2**a) % 2 != 1:
                    valid = False
                    break
                x = z//2**a
            if valid:
                hits.append(n)
            check(valid == (n == residue), "complete native cell")
        native_counts.append([k,mod,len(hits),residue])
    # A fixed family, k1..32, eight actual native members each.
    for k in range(1,33):
        p,q,b = paid_data(k)
        check(carrier(paid_word(k)) == (p,q,b), "closed carrier")
        fixed = Fraction(b,q-p)
        check(1 < fixed <= Fraction(23,5), "fixed point bound")
        for i in range(8):
            n = paid_source(k,i)
            r = paid_receipt(n,k)
            check(q*r.target == p*n+b, "affine endpoint")
            check(1 < r.target < n, "strict source payment")
            left = receipt(n,(1,)*k)
            right = receipt(left.target,(2,)*(2*k))
            check(concatenate(left,right) == r, "source-preserving block seam")
            check(recognize_paid(n) == r, "source determines its unique k")
    # Full affine identity, including translation, has both growth signs.
    quotient_rows = []
    for modulus in range(5,102):
        if gcd(modulus,6) != 1:
            continue
        k = neutral_length(modulus)
        growth = (1,)*k
        down = paid_word(k)
        check(affine_mod(growth,modulus) == (1,0), "full growing quotient identity")
        check(affine_mod(down,modulus) == (1,0), "full paid quotient identity")
        for t in range(1,4):
            n = 2**(k+1)*t-1
            h = replay(n,growth)
            check(h == 2*3**k*t-1 and h > n and (h-n)%modulus == 0,
                  "native growing identity family")
            source = paid_source(k,t-1)
            child = paid_receipt(source,k).target
            check(child < source and (child-source)%modulus == 0,
                  "native paid identity family")
        quotient_rows.append([modulus,k])
        for candidate_k in range(1,2*k+2):
            check((affine_mod(paid_word(candidate_k),modulus) == (1,0)) == (candidate_k%k == 0),
                  "complete neutral-index classification")
    # The same initial quotient state still permits opposite integer signs.
    growing_example = receipt(7,(1,1))
    paid_example = paid_receipt(8647,2)
    check((growing_example.target,paid_example.target) == (17,6157),
          "same-state opposite-sign endpoints")
    check(all(n%5 == 2 for n in (7,17,8647,6157)),
          "same initial and final quotient state")
    check(affine_mod(growing_example.word,5) == affine_mod(paid_example.word,5) == (1,0),
          "same full quotient map at the same quotient state")
    # Explicit infinite-complement law and a finite exact density control.
    for exponent in range(2,129):
        check(recognize_paid(2**exponent-1) is None, "Mersenne complement")
        k = exponent-1
        z = 3**k-1
        v = (z & -z).bit_length()-1
        check(v == (1 if k%2 else 2+((k & -k).bit_length()-1)) and v < 4*k,
              "LTE obstruction")
    bound = 8191
    observed = sum(recognize_paid(n) is not None for n in range(1,bound+1,2))
    predicted = 0
    for k in range(1,(bound+1).bit_length()):
        residue, modulus = native_cell(paid_word(k))
        if residue <= bound:
            predicted += 1+(bound-residue)//modulus
    check(observed == predicted, "disjoint native-union count")
    # Principal log validated against exhaustive exponent lists at small depths.
    for depth in range(1,8):
        modulus = 3**depth
        image = {pow(4,j,modulus):j for j in range(3**(depth-1))}
        check(len(image) == 3**(depth-1), "full principal group")
        for target, j in image.items():
            check(principal_log4(target,depth) == (j,3**(depth-1)), "independent discrete log")
    grounded_rows = []
    for k in (1,2):
        for index in range(3):
            plan = root_plan(k,index)
            n = materialize(plan)
            expected = paid_word(k)+(2*plan.j,)
            check(root_word_at_source(n,k) == expected, "unchanged-source ROOT membership")
            check(replay(n,expected) == 1, "literal first-hit ROOT")
            for modulus in (1,2,3,5,9,16,27,64,81,125,729,1001):
                check(source_residue(plan,modulus) == n%modulus, "full-denominator source reader")
                check(endpoint_residue(plan,modulus) == ((4**plan.j-1)//3)%modulus,
                      "endpoint reader")
            grounded_rows.append({"k":k,"index":index,"j":plan.j,
                                  "source_bits":n.bit_length(),"rank":3*k+1,
                                  "cost":5*k+2*plan.j})
    for k in range(1,17):
        for index in (0,1,10**20+7):
            plan = root_plan(k,index)
            p,q,b = paid_data(k)
            for modulus in (2,9,27,32,729,1001):
                n = source_residue(plan,modulus)
                h = endpoint_residue(plan,modulus)
                check((p*n+b-q*h)%modulus == 0, "symbolic carrier identity")
            residue, modulus = native_cell(paid_word(k))
            check(source_residue(plan,modulus) == residue, "symbolic native phase")
    check(root_plan(1) == RootPlan(1,0,10,27), "closed k1 phase")
    check(materialize(root_plan(1)) == 414251, "least displayed grounded source")
    check(root_word_at_source(43,1) is None, "paid does not mean this ROOT terminal")
    check(root_word_at_source(414253,1) is None, "adjacent source cannot borrow phase")
    check(discharge(paid_receipt(43,1),(4,1,1,2,3,4)).target == 1, "supplied independent terminal")
    check(carrier((1,2)) == (9,8,5) and carrier((2,1)) == (9,8,7), "order keeps carry")
    rejects(lambda: replay(11,(2,1)))
    rejects(lambda: concatenate(receipt(43,(1,)),receipt(49,(2,))))
    rejects(lambda: replay(1,(2,)))
    for n in (True,43.0,0,2,Fraction(43)):
        rejects(lambda n=n: paid_receipt(n,1))
    rejects(lambda: paid_receipt(45,1))
    rejects(lambda: audit_receipt(Receipt(43,37.0,(1,2,2))))
    rejects(lambda: audit_receipt(Receipt(43,35,(1,2,2))))
    rejects(lambda: audit_plan(replace(root_plan(1),k=True)))
    rejects(lambda: audit_plan(replace(root_plan(1),j=10.0)))
    rejects(lambda: audit_plan(replace(root_plan(1),j=11)))
    rejects(lambda: source_residue(root_plan(1),0))
    rejects(lambda: materialize(root_plan(10),100))
    rejects(lambda: neutral_length(6))
    print(json.dumps({
        "status":"PROVED scoped families; FINITE-EXACT controls; universal coverage OPEN",
        "checks":CHECKS,
        "native_cell_controls":native_counts,
        "paid_family_controls":{"k_inclusive":[1,32],"members_per_k":8},
        "native_union":{"relative_odd_density":"1/31","with_immediate_descent":"33/62",
                        "finite_odd_bound":bound,"finite_count":observed,
                        "mersenne_exponents_checked":[2,128]},
        "full_affine_identity_moduli_and_lengths":quotient_rows,
        "same_initial_state_control":{"modulus":5,"residue":2,"k":2,
                                      "growing":[7,17],"paid":[8647,6157]},
        "grounded_literal_controls":grounded_rows,
        "symbolic_controls":{"k_inclusive":[1,16],"indices":[0,1,10**20+7]},
        "closed_grounded_family":{"source":"(32*4^(10+27s)-101)/81",
                                  "word":"(1,2,2,20+54s)","s_minimum":0},
        "scope":"Construction varies the source. Group identity, payment and ROOT grounding remain separate predicates."
    },indent=2))


if __name__ == "__main__":
    main()
