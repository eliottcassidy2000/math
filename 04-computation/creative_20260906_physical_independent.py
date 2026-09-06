#!/usr/bin/env python3
"""Independent safe-wall cell sweep, no producer imports."""
from fractions import Fraction as F

H = (1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13)


def require(ok, why):
    if not ok:
        raise RuntimeError(why)


def norm(x):
    x %= 1
    return min(x, 1-x)


def safe(y):
    return all(norm(h*y) >= F(1, 14) for h in H)


def main():
    walls = {F(0), F(1)}
    for h in H:
        walls.update(F(14*k+s, 14*h) for k in range(h+1) for s in (-1, 1)
                     if 0 <= 14*k+s <= 14*h)
    walls = sorted(walls)
    intervals = [(a, b) for a, b in zip(walls, walls[1:]) if safe((a+b)/2)]
    points = [(y, y) for y in walls if safe(y)]
    components = []
    for a, b in sorted(intervals+points):
        if components and a <= components[-1][1]:
            components[-1] = (components[-1][0], max(b, components[-1][1]))
        else:
            components.append((a, b))
    expected = [(F(1,14),F(1,14)),(F(15,182),F(13,154)),
                (F(3,14),F(3,14)),(F(5,14),F(5,14)),
                (F(29,70),F(41,98)),(F(85,182),F(27,56)),
                (F(29,56),F(97,182)),(F(57,98),F(41,70)),
                (F(9,14),F(9,14)),(F(11,14),F(11,14)),
                (F(141,154),F(167,182)),(F(13,14),F(13,14))]
    require(components == expected, 'entire safe set including isolated points')
    mass = sum((b-a for a, b in components), F(0))
    width = max(b-a for a, b in components)
    require(mass == F(5939,140140) and width == F(11,728), 'mass and width')
    frontier = ((F(2,49),F(1,49)),(F(138,3325),F(37,3325)),
                (F(12,287),F(2,287)),(F(78,1855),F(2,371)),
                (F(20,469),F(2,469)))
    require(all(mass >= m or width >= b for m, b in frontier), 'joint entry')
    require(mass < F(20,469) and width < F(1,49), 'both scalar gates fail')
    row = (1, 67)+tuple(2*h for h in H)
    clearance = min(norm(v*F(9,34)) for v in row)
    require(clearance == F(2,17), 'literal full-row phase')
    print('PHYSICAL_BODY_INDEPENDENT_WALL_SWEEP')
    print(f'wall_points={len(walls)} intervening_cells={len(walls)-1}')
    print(f'components={len(components)} isolated={sum(a==b for a,b in components)}')
    print(f'mass={mass} width={width}')
    print(f'full_row_phase=9/34 clearance={clearance}')
    print('joint_certificate=PASS both_scalar_certificates=FAIL')
    print('PASS')


if __name__ == '__main__':
    main()
