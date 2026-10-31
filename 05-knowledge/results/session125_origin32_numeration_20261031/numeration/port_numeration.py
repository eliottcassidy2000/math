"""Session 125: exact Chung--Graham arithmetic on 1,3,8,21,... .

Public words are most-significant-digit first; coefficient arrays used by the
local rewrite engine are least-significant first. No floating point or assert
statements are used. These are integer encodings, not Hamiltonicity witnesses.
"""

from collections import deque


def _nat(n, name="n"):
    if not isinstance(n, int) or isinstance(n, bool):
        raise TypeError(f"{name} must be an integer")
    if n < 0:
        raise ValueError(f"{name} must be nonnegative")


def fibonacci(n):
    _nat(n)
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def weights(length):
    """Return A_0,...,A_(length-1), where A_i=F_(2i+2)."""
    _nat(length, "length")
    out = []
    a, b = 1, 3
    for _ in range(length):
        out.append(a)
        a, b = b, 3 * b - a
    return out


def valid(word, locked=False):
    """Leading zeros allowed; empty word represents zero."""
    if not isinstance(word, str):
        return False
    for char in word:
        if char not in "012" or (locked and char == "2"):
            return False
        if char == "0":
            locked = False
        elif char == "2":
            locked = True
    return True


def value(word):
    """Evaluate any 0/1/2 word, even if it is noncanonical."""
    if not isinstance(word, str) or any(c not in "012" for c in word):
        raise ValueError("word must be a string over 0,1,2")
    n, shifted = 0, 0
    for char in word:
        d = int(char)
        n, shifted = shifted + d, 3 * shifted - n + 3 * d
    return n


def encode(n):
    """Greedy Chung--Graham representation, with encode(0)='0'."""
    _nat(n)
    if n == 0:
        return "0"
    ws = [1]
    previous, current = 0, 1
    while 3 * current - previous <= n:
        previous, current = current, 3 * current - previous
        ws.append(current)
    digits = []
    for weight in reversed(ws):
        digit, n = divmod(n, weight)
        digits.append(str(digit))
    return "".join(digits)


def successor(word):
    """Lexicographic next legal word, independent of integer evaluation."""
    if not valid(word):
        raise ValueError("successor needs a legal Chung--Graham word")
    digits = list("0" + (word.lstrip("0") or "0"))
    maxima = []
    locked = False
    for digit in digits:
        maxima.append(1 if locked else 2)
        if digit == "0":
            locked = False
        elif digit == "2":
            locked = True
    for index in range(len(digits) - 1, -1, -1):
        if int(digits[index]) < maxima[index]:
            digits[index] = str(int(digits[index]) + 1)
            digits[index + 1:] = ["0"] * (len(digits) - index - 1)
            return "".join(digits).lstrip("0") or "0"
    raise RuntimeError("prepended zero should always permit a successor")


def normalize(coefficients, with_trace=False):
    """Normalize finite nonnegative LSD-first coefficients by local carries.

    The potential sum(coefficients) strictly decreases on every rewrite.
    Input coefficients need not be at most two. The optional trace records
    pre/post potentials and the applicable identity, not giant intermediate
    arrays. It provides evidence of the proved termination invariant.
    """
    ds = list(coefficients)
    for index, coefficient in enumerate(ds):
        _nat(coefficient, f"coefficient[{index}]")
    if not ds:
        ds = [0]
    trace = []

    def add(index, amount=1):
        while len(ds) <= index:
            ds.append(0)
        ds[index] += amount

    while True:
        before = sum(ds)
        high_digit = next((i for i, d in enumerate(ds) if d >= 3), None)
        if high_digit is not None:
            i = high_digit
            ds[i] -= 3
            add(i + 1)
            if i:
                add(i - 1)
            trace.append({"rule": "three", "i": i,
                          "before": before, "after": sum(ds)})
            continue
        lower_two = None
        offending_pair = None
        for i, d in enumerate(ds):
            if d == 0:
                lower_two = None
            elif d == 2:
                if lower_two is not None:
                    offending_pair = lower_two, i
                    break
                lower_two = i
        if offending_pair is None:
            break
        low, high = offending_pair
        ds[low] -= 2
        ds[high] -= 2
        for i in range(low + 1, high):
            ds[i] -= 1
        add(high + 1)
        if low:
            add(low - 1)
        trace.append({"rule": "two_ones_two", "low": low, "high": high,
                      "before": before, "after": sum(ds)})
    while len(ds) > 1 and ds[-1] == 0:
        ds.pop()
    answer = "".join(map(str, reversed(ds)))
    return (answer, trace) if with_trace else answer


def zeckendorf(n):
    """Canonical binary word on weights F_2,F_3,..., MSD first."""
    _nat(n)
    if n == 0:
        return "0"
    fs = [1, 2]
    while fs[-1] <= n:
        fs.append(fs[-1] + fs[-2])
    while fs[-1] > n:
        fs.pop()
    out = []
    for f in reversed(fs):
        digit = int(n >= f)
        out.append(str(digit))
        n -= digit * f
    return "".join(out)


def zeckendorf_value(word):
    if not isinstance(word, str) or any(c not in "01" for c in word):
        raise ValueError("Zeckendorf word must be binary")
    a, b, total = 1, 2, 0
    for digit in reversed(word):
        total += int(digit) * a
        a, b = b, a + b
    return total


def split_into_zeckendorf(word):
    """Burns' splitting bridge: d=min(d,1)+max(d-1,0) on even indices.

    This does not itself return the Zeckendorf word for the sum: ordinary
    Zeckendorf addition/normalization is still necessary.
    """
    if not valid(word):
        raise ValueError("split needs a legal Chung--Graham word")
    lower, upper = [], []
    for d in reversed(word):
        lower.extend((str(int(d != "0")), "0"))
        upper.extend((str(int(d == "2")), "0"))
    def trim(ds):
        return "".join(reversed(ds)).lstrip("0") or "0"
    return trim(lower), trim(upper)


def parity(word):
    """Return n mod 2 via the MSB finite-state arithmetic tracker."""
    if not valid(word):
        raise ValueError("parity needs a legal Chung--Graham word")
    n, shifted = 0, 0
    for char in word:
        d = int(char)
        n, shifted = (shifted + d) % 2, (n + shifted + d) % 2
    return n


def parity_states():
    """Reachable legality x parity states, omitting the rejection sink."""
    start = (False, 0, 0)
    witnesses = {start: ""}
    queue = deque([start])
    while queue:
        locked, n, shifted = queue.popleft()
        for d in range(2 if locked else 3):
            new_locked = False if d == 0 else locked or d == 2
            state = (new_locked, (shifted + d) % 2,
                     (n + shifted + d) % 2)
            if state not in witnesses:
                witnesses[state] = witnesses[(locked, n, shifted)] + str(d)
                queue.append(state)
    return witnesses
