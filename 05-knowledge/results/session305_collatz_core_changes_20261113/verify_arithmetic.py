"""Replay S305 arithmetic. No external packages; checks also under -O."""
def require(cond):
    if not cond:
        raise ValueError("Arithmetic certificate failed")


HOLES = {7, 21}
EXCEPTION = {
    9: ("-+", (9, 13, 5)),
    19: ("+-+++--+--", (19, 29, 43, 65, 49, 37, 55, 41, 31, 23, 17)),
    27: ("++-", (27, 41, 31, 23)),
    29: ("-+++--+-", (29, 43, 65, 49, 37, 55, 41, 31, 23)),
    37: ("--+", (37, 55, 41, 31)),
    75: ("++---", (75, 113, 85, 127, 95, 71)),
    117: ("---", (117, 175, 131, 49)),
}


def step(n, sign):
    z = 3*n + sign
    a = (z & -z).bit_length() - 1
    m = z >> a
    return m, a


def safe(m, a, sign):
    return m not in HOLES and all((m << i)-sign not in HOLES for i in range(1, a+1))


def desc(n):
    require(n > 1 and n % 2 and n not in HOLES | {5})
    s = 1 if n % 4 == 1 else -1
    m, a = step(n, s)
    if a >= 6:
        signs = [-s]*3
    elif safe(m, a, s):
        signs = [s]
    else:
        signs = [1 if t == "+" else -1 for t in EXCEPTION[n][0]]
    values = [n]
    vals = []
    for sign in signs:
        m, a = step(values[-1], sign)
        require(a <= 5 and safe(m, a, sign))
        values.append(m)
        vals.append(a)
    require(len(signs) <= 10 and values[-1] < n)
    require(19 * 3 * max(values[:-1]) <= 195 * n)
    return tuple(values), tuple(vals)


def main():
    bad_ladders = set()
    for h in HOLES:
        for s in [-1, 1]:
            d = h+s
            i = (d & -d).bit_length()-1
            m = d >> i
            if m % 3:
                for a in range(i, 6):
                    z = (m << a)-s
                    if z % 3 == 0 and z//3 not in HOLES:
                        bad_ladders.add(z//3)
    require(bad_ladders == {5, 27, 29, 117})
    for n, (_, vals) in EXCEPTION.items():
        require(desc(n)[0] == vals)
    require(step(5, -1)[0] == 7 and not safe(*step(5, 1), 1))
    for n in range(3, 100001, 2):
        if n not in HOLES | {5}:
            vals, vs = desc(n)
            if max(step(n, 1)[1], step(n, -1)[1]) >= 6:
                require(vs == (1, 2, 2))
    layer = {19}
    for expected in [29, 43, 65, 49, 37, 55, 41, 31, 23, 17]:
        layer = {m for n in layer for s in (-1, 1)
                 for m, a in [step(n, s)] if a <= 5 and safe(m, a, s)}
        require(min(layer) == expected)
    for length in range(1, 101):
        mod = 1 << (length+3)
        n = (-pow(5, -1, mod)) % mod + mod
        for i in range(length):
            s = 1 if i % 2 == 0 else -1
            y, a = step(n, s)
            require(a == 1 and y > n and step(n, -s)[1] == 3)
            n = y
    print("PASS: exact exceptions, 100000 control, bounded-lookahead obstructions")


if __name__ == "__main__":
    main()
