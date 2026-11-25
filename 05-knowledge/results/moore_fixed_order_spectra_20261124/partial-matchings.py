from math import prod
for preset, exponent in [(0,16256),(1,16241)]:
    m = 1250
    b = m + 1 - 4*59 - 3
    e = 52929
    match_orders = prod((l*(l-375)+49)//50
                        for l in range(376,1626-preset))
    num, den = match_orders*b**e, (b+1)**e
    assert den << exponent < num < den << (exponent+1)
    print('forced backward pairs:',preset,'exponent:',exponent)
