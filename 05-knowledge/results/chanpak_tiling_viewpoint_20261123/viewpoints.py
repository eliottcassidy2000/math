"""Finite certificates for the path-permutation group and the one-backward-arc family."""
from itertools import product


def path_words(rows):
    full = (1 << len(rows))-1

    def extend(path, available):
        if not available:
            yield path
            return
        choices = available if not path else available & rows[path[-1]]
        while choices:
            bit = choices & -choices
            choices -= bit
            yield from extend(path + (bit.bit_length()-1,), available ^ bit)

    yield from extend((), full)


def odd(word):
    return sum(a > b for i, a in enumerate(word) for b in word[i+1:]) % 2


def board(rows, word):
    return tuple((rows[word[i]] >> word[j]) & 1 for i in range(len(word))
                 for j in range(i+2, len(word)))


def window_certificate():
    # Labels 0,...,5 follow the rim. Apex 6 beats 1,3,5, loses to 0,2,4.
    # The baseline (0,6,1,2,3,4,5) is odd. An even path beginning in 0 or
    # ending in 5 can use the unchanged remainder of the rim.
    chords = [(i,j) for i in range(6) for j in range(i+2,6)]
    survivors = []
    for mask in range(1 << len(chords)):
        rows = [2,4,8,16,32,0,0]
        for i in range(6):
            if i % 2:
                rows[6] |= 1 << i
            else:
                rows[i] |= 1 << 6
        for bit, (i,j) in enumerate(chords):
            if (mask >> bit) & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
        if not any((p[0] == 0 or p[-1] == 5) and not odd(p)
                   for p in path_words(rows)):
            survivors.append(mask)
    assert survivors == [855, 863, 967, 975]
    assert all(mask & 1 and mask & (1 << 9) for mask in survivors)
    return survivors


def small_strong(n):
    # Enumerate the whole staircase, then quotient by its actual path presentations.
    chords = [(i,j) for i in range(n) for j in range(i+2,n)]
    classes = {}
    for flags in product((0,1), repeat=len(chords)):
        rows = [1 << (i+1) for i in range(n-1)] + [0]
        for forward,(i,j) in zip(flags,chords):
            a,b = (i,j) if forward else (j,i)
            rows[a] |= 1 << b
        paths = list(path_words(rows))
        closed = [p for p in paths if (rows[p[-1]] >> p[0]) & 1]
        if not closed:
            continue
        views = {board(rows,p) for p in paths}
        # Reflection of chord (i,j) is (n-1-j,n-1-i).
        reflected = [chords.index((n-1-j,n-1-i)) for i,j in chords]
        blue = sum(tuple(t[k] for k in reflected) == t for t in views)
        key = min(views)
        value = (len(paths), len(views), blue, len({odd(p) for p in paths}))
        assert classes.get(key, value) == value
        classes[key] = value
    return sorted(classes.values())


def backward_arc(n):
    rows = [sum(1 << j for j in range(i+1,n)) for i in range(n)]
    rows[0] ^= 1 << (n-1)
    rows[n-1] ^= 1
    paths = list(path_words(rows))
    closed = [p for p in paths if rows[p[-1]] >> p[0] & 1]
    signed = lambda family: sum(1-2*odd(p) for p in family)
    assert (len(paths), signed(paths)) == (1+2**(n-2), 1+(-1)**(n-1)*2**((n-1)//2))
    assert (len(closed), signed(closed)) == (n, n if n%2 else 0)
    blue = [p for p in paths if all(p[i] == n-1-p[-1-i] for i in range(n))]
    assert len(blue) == 1 + (2**((n-2)//2) if n%2 == 0 else 0)
    assert signed(blue) == 1


if __name__ == "__main__":
    print("Six-vertex window masks", window_certificate())
    for n in (3,4,5):
        classes = small_strong(n)
        print("Strong order", n, "(paths, tilings, blue, number of signs)", classes)
        assert len(classes) == {3:1,4:1,5:6}[n]
        assert [r for r in classes if r[-1]==1] == {3:[(3,1,1,1)],4:[],5:[(15,5,1,1)]}[n]
    for n in range(4,13):
        backward_arc(n)
    print("One-backward-arc formulas verified through order 12")
