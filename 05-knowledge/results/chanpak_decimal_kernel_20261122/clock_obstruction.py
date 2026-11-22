"""Check the Farey clock obstructions, their finite replacements, and the E22 seam."""

from math import gcd, lcm, prod

from e15_pair import D, EDGE, flipped, join, literal_count, ordinal, p, star_rows, strong, tagged_data
from fermat_receiver import address_height, literal_receiver, single_runs
from local_word_inverse import steer


SUPERSINGULAR = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71)
E61_LOWER = (1 << 60) * ((1 << 61) - 1) - 1
# (left rows, left tags, right replacement rows, replacement tags)
FIXTURES = (
    ((2, 4, 33, 7, 15, 27), 4, (4, 9, 18, 21, 35, 15), 1),
    ((62, 60, 56, 48, 32, 64, 31), 7, (4, 5, 96, 7, 15, 27, 59), 18),
)


def order_two(odd_modulus):
    exponent, a = 1, 2 % odd_modulus
    while a != 1:
        a = a * 2 % odd_modulus
        exponent += 1
    return exponent


def check_short_simultaneous_word(chart, x, y):
    odd_modulus = prod(SUPERSINGULAR)//2
    h = chart[0]
    common = gcd(h, odd_modulus)
    base = address_height(chart, x, y)-h*p(y)
    assert (E61_LOWER-base) % common == 0
    reduced = odd_modulus//common
    target_register = (E61_LOWER-base)//common * pow(2*h//common, -1, reduced) % reduced
    word = "1" + steer(odd_modulus, (1,0), (y,target_register))
    value = address_height(chart, x, int(word, 2))
    assert len(word) == 7562 and (value-E61_LOWER) % (2*odd_modulus) == 0
    print("short simultaneous word", len(word), "bits at height coefficient", h)


def check_clocks():
    full_modulus = prod(SUPERSINGULAR)
    odd_modulus = full_modulus // 2
    repeat_width = lcm(*(order_two(q) for q in SUPERSINGULAR[1:]))
    assert pow(2, repeat_width, odd_modulus) == 1
    for left, tags, right, rtags in FIXTURES:
        a, v = tagged_data(left, tags), tagged_data(right, rtags)
        den, num = a[0], a[0]-a[4]
        right_height, right_numerator = v[0], v[5]
        assert den*right_numerator-num*right_height == -1
        period = order_two(den)
        assert all((num * (pow(2, n, den)+1)) % den not in (1, den-1)
                   for n in range(1, period+1))
        ad = join(a, tagged_data(D, 35))
        flipped_ad = join(a, tagged_data(flipped(D, *EDGE), 35))
        local = ((3,2,0,1), (5,3,3,0)) if den == 15 else ((3,2,0,1), (11,1,0,2))
        for prime, n, x, y in local:
            rv = tagged_data(*literal_receiver((False,)*n))
            lower_local = address_height(join(flipped_ad, rv), x, y)
            upper_local = address_height(join(ad, rv), x, y)
            assert (lower_local-E61_LOWER) % prime == (upper_local-E61_LOWER-2) % prime == 0
        lower = join(flipped_ad, v)
        upper = join(ad, v)
        assert lower[:-1] == upper[:-1] and lower[-1]+2 == upper[-1]
        assert gcd(*lower[:4]) == 1
        common = gcd(odd_modulus, lower[0])
        x, y = next((x, y) for x in range(0, 2*common, 2) for y in range(0, 2*common, 2)
                    if (E61_LOWER-address_height(lower, x, y)) % common == 0)
        offset = address_height(lower, x, y)
        copies = ((E61_LOWER-offset)//common * pow(lower[0]//common, -1, odd_modulus//common)
                  % (odd_modulus//common))
        assert (offset+lower[0]*copies-E61_LOWER) % odd_modulus == 0
        assert (offset-E61_LOWER) % 2 == 0
        assert not single_runs(lower, E61_LOWER)
        check_short_simultaneous_word(lower, x, y)
        print("denominator", den, "A data", a, "A+D", ad, "replacement", v,
              "lower chart", lower, "reset (x,y,width,copies)", (x,y,repeat_width,copies))


def check_e22():
    left = (2, 4, 1)
    middle = (4, 17, 34, 7, 13, 27)
    right = (310, 364, 56, 465, 482, 456, 389, 7, 132)
    mark = (2,3)
    tag = 1 | (22 << len(left)) | (420 << (len(left)+len(middle)))
    x, y = 0, 185591
    chart_inputs = tagged_data(left,1), tagged_data(right,420)
    center = (1 << 21)*((1 << 22)-1)
    hosts, profiles = [], []
    for side, mid in enumerate((middle, flipped(middle, *mark))):
        lower_left, upper_right = chart_inputs
        chart = join(join(lower_left, tagged_data(mid,22)), upper_right)
        host = star_rows(ordinal(ordinal(left, mid), right), tag, x, y)
        value = address_height(chart, x, y)
        # Split the transitive suffix into singleton blocks to retain literal contact
        # counts without a subset DP of order its entire width.
        lengths = (len(left), len(middle), len(right)) + (1,)*y.bit_length()
        assert literal_count(host, lengths) == value == center-1+2*side
        assert len(host) == 37 and strong(host)
        hosts.append(host)
        profiles.append(chart)
    assert flipped(hosts[0], 5, 6) == hosts[1]
    assert profiles[0][:-1] == profiles[1][:-1] and profiles[0][-1]+2 == profiles[1][-1]
    assert chart_inputs[0][0]*chart_inputs[1][5]-(chart_inputs[0][0]-chart_inputs[0][4])*chart_inputs[1][0] == 1
    assert gcd(*profiles[0][:4]) == 3 and profiles[0][0]%9 != 0
    assert tuple(v%3 for v in profiles[0]) == (0,0,0,0,0,1,2)
    base = address_height(profiles[0], 1, 36)
    assert base == 6557140035
    copies = 3332926730265235
    assert (base+profiles[0][0]*copies-E61_LOWER) % (prod(SUPERSINGULAR)//2) == 0
    assert not single_runs(profiles[0], E61_LOWER)
    check_short_simultaneous_word(profiles[0], 1, 36)
    print("E22", center, "base profiles", profiles, "physical arc (5,6)")


if __name__ == "__main__":
    check_e22()
    check_clocks()
