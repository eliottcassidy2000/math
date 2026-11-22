"""Test all single-run insertions in the three exponent-31 certificates."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

source = Path(__file__).resolve().parents[1] / 'chanpak_e31_physical_20261122' / 'check.py'
spec = spec_from_file_location('exponent31', source)
e31 = module_from_spec(spec)
spec.loader.exec_module(e31)

start, target = (e31.perfect_shape(p) for p in (31, 61))
hits = []
for index, (core, tag, edge, reader, reader_tag, y, order) in enumerate(e31.WITNESSES):
    deletion = e31.ordinal(core, reader)
    mark = tag | (reader_tag << len(core))
    data = e31.profile(deletion, mark)
    twin = e31.profile(e31.flipped(deletion, *edge), mark)
    h, c = data[0], sum(data[2])
    assert h == twin[0] and c == sum(twin[2])
    for s in range(1, y.bit_length() + 1):
        b = 1 << s
        u, v = divmod(y, b)
        reverse = int(format(v, f'0{s}b')[::-1], 2)
        for digit in (0, 1):
            if digit:
                slope = (u+1)*(c*b + h*(b-1-reverse))
                extra = h * (u % 2 == 0 and 2*v < b)
            else:
                slope = u*(c*b + h*(b-reverse))
                extra = h * (u % 2 == 1 and 2*v >= b)
            for width in (1, 2, 6):
                address = ((((u+digit) << width)-digit) << s) + v
                center = start + extra + slope*((1 << width)-1)
                assert e31.count_from_profile(data, 0, address) == center-1
                assert e31.count_from_profile(twin, 0, address) == center+1
            for m0 in range(12):
                if (start + extra + slope*((1 << m0)-1) - target) % 4095 == 0:
                    hits.append((index, s, digit, slope, extra, m0))
slope = 2304679764206100516
assert hits == [(2, s, 1, slope, 0, 6) for s in (4, 5, 6)]
for prime, valuation in ((2, 2), (3, 4), (5, 1), (7, 1), (13, 1)):
    difference = target - (start + 63*slope)
    assert difference % (prime**valuation) == 0 and difference % (prime**(valuation+1))
print('Three equivalent cuts give the unique ray, at width 6 modulo 12.')
print('First physical center:', start + 63*slope, 'order:', 66)
