"""Independent ordinary-Collatz realization filter for golden phase cycles.

The phase cycles and their periodic golden lifts are supplied by the sibling
module. This script does not import its arithmetic decoder: ordinary branch
composition, rational fixed points, and full cyclic parity replay are rebuilt.
Run: python -X utf8 -B 04-computation/experiments/denominator_arithmetic_filter_20261004.py
"""

from collections import Counter
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path
import json


MODULI = tuple(range(1, 41)) + (64, 76, 81, 105)


def need(condition, message):
    if not condition:
        raise ValueError(message)


def load_phase_module():
    path = Path(__file__).with_name("denominator_fibre_recursion_20261004.py")
    spec = spec_from_file_location("golden_phase_source", path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def pair_sign(a, b):
    """Exact sign of a+b*phi; uses only integer comparisons."""
    c = 2*a+b
    if b == 0:
        return (c > 0)-(c < 0)
    if c == 0 or (c > 0) == (b > 0):
        return (b > 0)-(b < 0)
    squared_difference = c*c-5*b*b
    return ((c > 0)-(c < 0))*((squared_difference > 0)-(squared_difference < 0))


def golden_word(q, initial, expected_phases):
    a, b = initial
    states, bits = [], []
    for phase in expected_phases:
        need((a % q, b % q) == phase, "golden lift follows supplied phase orbit")
        need(gcd(a,b,q) == 1, "exact scalar golden denominator")
        states.append((a,b))
        digit = int(pair_sign(b-q,a+b) > 0)  # Upper boundary convention.
        bits.append(digit)
        a,b = b-digit*q,a+b
    need((a,b) == initial, "golden cycle closes at exact phase period")
    return tuple(bits), tuple(states)


def arithmetic_fixed_point(bits):
    """Compose T_1(x)=3x+1 and T_0(x)=x/2, then solve T_word(x)=x."""
    power3, zeros, carry = 1, 0, 0
    for bit in bits:
        need(bit in (0,1), "binary input")
        if bit:
            power3 *= 3
            carry = 3*carry+2**zeros
        else:
            zeros += 1
    denominator = 2**zeros-power3
    need(denominator != 0, "nonempty word has nonzero return denominator")
    return Fraction(carry,denominator), dict(ones=sum(bits),zeros=zeros,
                                           carry=carry,raw_denominator=denominator)


def rotate_to_smallest_absolute(orbit, bits):
    anchor = min(range(len(orbit)),key=lambda i:(abs(orbit[i]),orbit[i]))
    return orbit[anchor:]+orbit[:anchor], bits[anchor:]+bits[:anchor]


def decode_and_check(bits):
    need(bool(bits), "nonempty cyclic word")
    need(all(not(bits[i] and bits[(i+1)%len(bits)]) for i in range(len(bits))),
         "cyclic golden word has no adjacent ones")
    root, transfer = arithmetic_fixed_point(bits)
    state = root
    orbit = []
    for bit in bits:
        need(state.denominator % 2 == 1, "rational state has odd denominator")
        need(state.numerator % 2 == bit, "actual rational 2-adic parity agrees with prescribed bit")
        orbit.append(state)
        state = 3*state+1 if bit else state/2
    need(state == root, "literal rational cycle closes")
    need(len(set(orbit)) == len(bits), "arithmetic cycle has exact golden word period")
    denominators = {x.denominator for x in orbit}
    need(len(denominators) == 1 and gcd(next(iter(denominators)),6) == 1,
         "rational cycle denominator is constant and coprime to six")
    integral = all(x.denominator == 1 for x in orbit)
    need(integral == (root.denominator == 1), "root and every-state integrality agree")
    orbit, canonical_bits = rotate_to_smallest_absolute(tuple(orbit),tuple(bits))
    return dict(root=orbit[0],orbit=orbit,bits=canonical_bits,denominator=next(iter(denominators)),
                integral=integral,transfer=transfer)


def odd_compression(decoded):
    """Keep odd states and the exact number of zero branches after each one."""
    bits, orbit = decoded["bits"], decoded["orbit"]
    positions = [i for i,bit in enumerate(bits) if bit]
    need(bool(positions), "zero cycle has no odd compression")
    word = tuple((positions[(i+1)%len(positions)]-position-1) % len(bits)
                 for i,position in enumerate(positions))
    odd_states = tuple(orbit[i] for i in positions)
    for i,(state,a) in enumerate(zip(odd_states,word)):
        need(a >= 1, "positive odd valuation")
        need((3*state+1)/2**a == odd_states[(i+1)%len(word)],
             "compressed rational odd cycle closes with exact valuations")
    return word,odd_states


def favorable_rotation(word):
    """Rotate an expanding cycle after its unique minimum prefix product.

    Fraction comparison is exact: no logarithm or floating-point threshold
    is used to select the minimum or certify any prefix.
    """
    need(bool(word) and all(a >= 1 for a in word), "nonempty valuation word")
    need(3**len(word) > 2**sum(word), "negative anchor has expanding full multiplier")
    ratios, A = [], 0
    for i,a in enumerate(word):
        ratios.append(Fraction(3**i,2**A))
        A += a
    need(len(set(ratios)) == len(ratios), "partial products cannot tie")
    cut = min(range(len(word)),key=ratios.__getitem__)
    rotated = word[cut:]+word[:cut]
    A = 0
    for i,a in enumerate(rotated,1):
        A += a
        need(3**i > 2**A, "every rotated nonempty prefix expands")
    return cut,rotated


def ordinary_odd_step(n):
    need(n > 0 and n % 2 == 1, "positive odd integer input")
    out,a = 3*n+1,0
    while out % 2 == 0:
        out //= 2
        a += 1
    return out,a


def rational_anchor_control():
    """Finite replay of the inherited -19/11 compiler's first cylinder."""
    least_orbit,least_word = None,None
    for k in range(64):
        source = 999+1024*k
        state,orbit,word = source,[source],[]
        for _ in range(6):
            state,a = ordinary_odd_step(state)
            orbit.append(state)
            word.append(a)
        need(tuple(word[:5]) == (1,1,2,1,1) and word[-1] >= 4,
             "exact rational-shadow prefixes and guarded final exit")
        need(all(x > source for x in orbit[1:-1]) and 0 < orbit[-1] < source,
             "first odd descent is at step six")
        if k == 0:
            least_orbit,least_word = orbit,word
    need(least_orbit == [999,1499,2249,1687,2531,3797,89], "least cylinder source control")
    need(least_word == [1,1,2,1,1,7], "least source exact valuation word")
    weak_state,weak_orbit = 487,[487]
    for _ in range(6):
        weak_state,_ = ordinary_odd_step(weak_state)
        weak_orbit.append(weak_state)
    need(weak_orbit == [487,731,1097,823,1235,1853,695], "weak-budget hostile")
    need(all(x > 487 for x in weak_orbit[1:]), "legal shadow alone does not pay for descent")
    return dict(cylinder="999 mod 1024",finite_lifts_checked=64,
                least_source_odd_orbit=least_orbit,least_source_valuation_word=least_word,
                weak_budget_hostile=weak_orbit)


def compile_rational_anchor(word,anchor,q):
    """Instantiate the inherited compiler, then independently replay three lifts."""
    need(anchor < 0 and anchor.denominator > 1, "negative rational-only anchor")
    h,d = -anchor.numerator,anchor.denominator
    P,Q,B = 1,1,0
    for a in word:
        P,B,Q = 3*P,3*B+Q,Q*2**a
        need(P > Q, "compiler requires every nonempty prefix to expand")
    need(Fraction(-B,P-Q) == anchor, "chart fixed point agrees with rational cycle")
    g = gcd(B,P-Q)
    need((h,d) == (B//g,(P-Q)//g) and gcd(h*d,6) == 1,
         "reduced positive anchor numerator and denominator are units at two and three")
    m = 1
    while Q**m <= h:
        m += 1
    Pm,Qm = P**m,Q**m
    t = 1
    while 2**t*(Qm-h) <= Pm-h:
        t += 1
    A = sum(word)
    modulus = 2**(A*m+t)
    beta = (h*pow(Pm,-1,2**t)) % 2**t
    residue = ((beta*Qm-h)*pow(d,-1,modulus)) % modulus
    need(residue > 0 and residue % 2 == 1, "least positive cylinder source is odd")
    endpoints,extra_halvings = [],[]
    nominal = word*m
    for lift in range(3):
        source = residue+lift*modulus
        need((d*source+h) % Qm == 0, "source decodes an integral displacement coefficient")
        b = (d*source+h)//Qm
        need(b > 0 and b % 2**t == beta, "source satisfies complete compiler congruence")
        z = b*Pm-h
        reduced_z,tau = z,0
        while reduced_z % 2 == 0:
            reduced_z //= 2
            tau += 1
        need(tau >= t and reduced_z % d == 0, "guarded endpoint is integral")
        endpoint = reduced_z//d
        state = source
        actual_word,orbit = [],[source]
        for _ in nominal:
            state,a = ordinary_odd_step(state)
            actual_word.append(a)
            orbit.append(state)
        need(tuple(actual_word) == nominal[:-1]+(nominal[-1]+tau,),
             "literal iteration agrees with the complete compiled valuation word")
        need(state == endpoint and 0 < endpoint < source,
             "affine endpoint and literal replay pay the original-source budget")
        need(all(x > source for x in orbit[1:-1]), "compiled exit is the exact first odd descent")
        endpoints.append(endpoint)
        extra_halvings.append(tau)
    return dict(golden_denominator=q,anchor=str(anchor),word=list(word),h=h,d=d,
                m=m,t=t,residue=residue,modulus=modulus,modulus_exponent=A*m+t,
                first_descent_odd_step=len(word)*m,endpoint=endpoints[0],
                three_lift_endpoints=endpoints,three_lift_extra_halvings=extra_halvings)


def main():
    phase = load_phase_module()
    rows, integers, smallest_hostiles, all_denominators = [], [], [], Counter()
    cycle_count = phase_points = rational_states = 0
    signs,negative_denominators,negative_odd_periods = Counter(),Counter(),Counter()
    negative_examples,retained_anchor = [],None
    rotated_phase_changes = 0
    arbitrary_negative_phases_checked = arbitrary_phase_changes = 0
    compiler_instances = []
    seen_words = set()
    for q in MODULI:
        candidates = []
        if q == 1:
            # Three zero-phase lifts produce two cycles; canonical lift alone keeps only zero.
            candidates = [((0,),((0,0),),"zero"), ((1,0),((1,0),(-1,1)),"upper boundary")]
        else:
            for phases in phase.phase_cycles(q):
                x = phase.periodic_lift(q,phases[0])
                initial = tuple(int(value*q) for value in x)
                need(all(Fraction(initial[i],q) == x[i] for i in (0,1)), "integral numerator lift")
                bits,states = golden_word(q,initial,phases)
                candidates.append((bits,states,"primitive phase"))
                phase_points += len(phases)
        arithmetic_denominators = Counter()
        q_signs = Counter()
        integer_roots = []
        for bits,states,tag in candidates:
            word_key = min(bits[i:]+bits[:i] for i in range(len(bits)))
            need(word_key not in seen_words, "distinct golden cycles have distinct cyclic parity words")
            seen_words.add(word_key)
            decoded = decode_and_check(bits)
            cycle_count += 1
            rational_states += len(bits)
            denominator = decoded["denominator"]
            arithmetic_denominators[denominator] += 1
            all_denominators[denominator] += 1
            root = decoded["root"]
            sign = "zero" if root == 0 else "negative" if root < 0 else "positive"
            kind = sign+("_integer" if decoded["integral"] else "_rational_only")
            signs[kind] += 1
            q_signs[kind] += 1
            need(all((x > 0)-(x < 0) == (root > 0)-(root < 0) for x in decoded["orbit"]),
                 "every cycle has one arithmetic sign")
            if root < 0:
                word,odd_states = odd_compression(decoded)
                cut,rotated = favorable_rotation(word)
                rotated_phase_changes += int(cut != 0)
                need(decoded["bits"][0] == 1 and cut == 0,
                     "smallest-absolute negative state is already an expanding-prefix odd anchor")
                for shift in range(len(word)):
                    shifted = word[shift:]+word[:shift]
                    shifted_cut,_ = favorable_rotation(shifted)
                    arbitrary_negative_phases_checked += 1
                    arbitrary_phase_changes += int(shifted_cut != 0)
                anchor = odd_states[cut]
                formal = anchor
                for a in rotated:
                    formal = (3*formal+1)/2**a
                need(formal == anchor, "favorable phase retains the same rational cycle")
                if not decoded["integral"]:
                    compiler_instances.append(compile_rational_anchor(rotated,anchor,q))
                    negative_denominators[denominator] += 1
                    negative_odd_periods[len(word)] += 1
                    if len(negative_examples) < 8:
                        negative_examples.append(dict(golden_denominator=q,root=str(root),
                            arithmetic_denominator=denominator,odd_word=word,rotation_index=cut,
                            favorable_anchor=str(anchor),favorable_word=rotated))
                if root == Fraction(-19,11):
                    expected_bits = tuple(map(int,"1010100"))
                    need(q == 29 and decoded["bits"] == expected_bits and word == (1,1,2),
                         "retained rational anchor has the claimed exact golden fibre and word")
                    offset = next(i for i in range(len(bits)) if bits[i:]+bits[:i] == expected_bits)
                    need(states[offset] == (-4,20), "Theta(-19/11)=(-4+20*phi)/29")
                    retained_anchor = dict(golden_denominator=q,root=str(root),
                        arithmetic_denominator=denominator,parity_word="1010100",
                        golden_numerator_pair=list(states[offset]),odd_word=word,
                        odd_cycle=[str(x) for x in odd_states],favorable_rotation_index=cut)
            if decoded["integral"]:
                root = int(decoded["root"])
                integer_roots.append(root)
                integers.append(dict(golden_denominator=q,tag=tag,root=root,
                                     ordinary_period=len(bits),orbit=[int(x) for x in decoded["orbit"]],
                                     parity_word="".join(map(str,decoded["bits"])),arithmetic_denominator=1))
            elif len(smallest_hostiles) < 8:
                smallest_hostiles.append(dict(golden_denominator=q,
                    root=str(decoded["root"]),arithmetic_denominator=denominator,
                    ordinary_period=len(bits),orbit=[str(x) for x in decoded["orbit"]],
                    parity_word="".join(map(str,decoded["bits"]))))
        rows.append(dict(q=q,cycles=len(candidates),points=sum(len(c[0]) for c in candidates),
                         periods=dict(sorted(Counter(len(c[0]) for c in candidates).items())),
                         integer_roots=sorted(integer_roots),
                         arithmetic_sign_counts=dict(sorted(q_signs.items())),
                         arithmetic_denominators={str(d):c for d,c in sorted(arithmetic_denominators.items())}))

    actual = {(row["golden_denominator"],row["root"]) for row in integers}
    need(actual == {(1,0),(1,-1),(2,1),(11,-5),(76,-17)}, "integer cycles in the declared finite denominator universe")
    need(sum(row["cycles"] for row in rows) == cycle_count, "cycle census total")
    need(sum(row["points"] for row in rows) == rational_states == phase_points+3, "boundary-adjusted point total")
    need(retained_anchor is not None, "q=29 retained rational anchor appears in complete census")
    need(len(compiler_instances) == signs["negative_rational_only"],
         "every negative rational-only cycle has a replayed compiler instance")
    need(sum(signs.values()) == cycle_count, "arithmetic sign census")
    positive_word,_ = odd_compression(decode_and_check((1,0,0)))
    need(positive_word == (2,) and 3 < 2**sum(positive_word),
         "positive cycle 1 is a hostile to universal expanding-prefix rotation")
    report = dict(status="FINITE-EXACT arithmetic-realization filter; no all-denominator or global Collatz coverage claim",
                  golden_denominator_universe=list(MODULI),
                  counts=dict(golden_cycles=cycle_count,golden_periodic_points=rational_states,
                              primitive_nonzero_phase_points=phase_points,
                              integer_cycles=len(integers),rational_only_cycles=cycle_count-len(integers),
                              distinct_arithmetic_denominators=len(all_denominators)),
                  arithmetic_sign_counts=dict(sorted(signs.items())),
                  negative_anchor_gate=dict(
                      status="PROVED inherited favorable-rotation lemma; FINITE-EXACT application to this census",
                      cycles_checked=signs["negative_integer"]+signs["negative_rational_only"],
                      rational_only_cycles_checked=signs["negative_rational_only"],
                      canonical_odd_phases_changed=rotated_phase_changes,
                      arbitrary_cyclic_odd_phases_checked=arbitrary_negative_phases_checked,
                      arbitrary_odd_phases_changed=arbitrary_phase_changes,
                      rational_only_odd_period_histogram=dict(sorted(negative_odd_periods.items())),
                      rational_only_denominator_histogram={str(d):c for d,c in sorted(negative_denominators.items())},
                      first_rational_examples=negative_examples,
                      retained_q29_anchor=retained_anchor,
                      compiler_instances=compiler_instances,
                      compiler_integer_source_replays=3*len(compiler_instances),
                      inherited_compiler_control=rational_anchor_control(),
                      preserved="The cyclic word and rational cycle, with a phase change retaining every exact valuation.",
                      additional_guard="A positive-integer source, congruence guards, and original-source exit budget are still required."),
                  all_integer_cycles=integers,
                  first_rational_hostiles=smallest_hostiles,
                  exact_census_by_golden_denominator=rows,
                  missing_gate="A periodic golden cycle carries a rational parity realization; integrality is an additional tested predicate.")
    print(json.dumps(report,indent=2))
    print("PASS: all acceptance checks remain active under -O")


if __name__ == "__main__":
    main()
