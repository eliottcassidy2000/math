"""Exact controls for reset_composition.md; no third-party dependencies."""
from fractions import Fraction as F
from itertools import product
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def data(word):
    A, Q, B = 1, 1, 0
    a, p, d = F(1), F(0), F(0)
    for letter in word:
        require(type(letter) is int and letter >= 1, "invalid letter")
        slope = F(3, 2**letter)
        A, Q, B = 3*A, 2**letter*Q, 3*B+Q
        a, p, d = (slope*a, slope*p+(letter == 1),
                   slope*d+(1-F(4, 2**letter) if letter >= 3 else 0))
    require(a == F(A, Q), "multiplier metadata")
    require(p-d == F(B, Q)+a-1, "shifted constant")
    residue = ((Q-B)*pow(A, -1, 2*Q)) % (2*Q)
    return A, Q, B, a, p, d, residue


def suffix_data(word):
    K = sum(word)
    running = 0
    p, d = F(0), F(0)
    for i, letter in enumerate(word):
        running += letter
        suffix = F(3**(len(word)-i-1), 2**(K-running))
        if letter == 1:
            p += suffix
        elif letter >= 3:
            d += (1-F(4, 2**letter))*suffix
    return F(3**len(word), 2**K), p, d


def replay(n, word):
    require(type(n) is int and n > 0 and n % 2 == 1, "invalid source")
    for letter in word:
        t = 3*n+1
        actual = (t & -t).bit_length()-1
        require(actual == letter, "source outside exact cylinder")
        n = t >> letter
    return n


def first_descent(n):
    word, h = [], n
    for _ in range(10000):
        t = 3*h+1
        letter = (t & -t).bit_length()-1
        word.append(letter)
        h = t >> letter
        if h < n:
            return tuple(word), h
    raise RuntimeError("first-descent control cap reached")


def main():
    words = [w for m in range(6) for w in product(range(1, 5), repeat=m)]
    cuts, cylinder_replays, forward, reciprocal = 0, 0, 0, 0
    for w in words:
        A, Q, B, a, p, d, r = data(w)
        require((a, p, d) == suffix_data(w), "independent suffix formula")
        require((Q*d).denominator == 1 and int(Q*d) % 4 == 0,
                "debt numerator divisibility")
        if w:
            require(d != a-1 and d != 1-a, "nonempty equality boundary")
        forward += d <= a-1
        reciprocal += d <= 1-a
        for lift in range(3):
            n = r+2*Q*lift
            endpoint = replay(n, w)
            require(endpoint == F(A*n+B, Q), "native endpoint")
            require(endpoint-1 == a*(n-1)+p-d, "native shifted endpoint")
            cylinder_replays += 1
        for cut in range(len(w)+1):
            u, v = w[:cut], w[cut:]
            Au, Qu, Bu, au, pu, du, ru = data(u)
            Av, Qv, Bv, av, pv, dv, rv = data(v)
            require((a, p, d) == (av*au, av*pu+pv, av*du+dv), "cocycle")
            require(a-1-d == av*(au-1-du)+(av-1-dv), "forward slack")
            require(1-a-d == av*(1-au-du)+(1-av-dv), "inverse slack")
            require(d/a == du/au+(1/au)*(dv/av), "normalized inverse law")
            if du <= au-1 and dv <= av-1:
                require(d <= a-1, "forward cone closure")
            if du <= 1-au and dv <= 1-av:
                require(d <= 1-a, "reciprocal cone closure")
            su = replay(ru, u)
            require((rv-su) % 2 == 0, "odd endpoint guard")
            t = ((rv-su)//2*pow(Au, -1, Qv)) % Qv if Qv != 1 else 0
            require((ru+2*Qu*t) % (2*Q) == r, "domain pullback")
            cuts += 1
    hostile_words = [(1, 1, 1, 3), (1, 1, 1, 1, 3),
                     (3, 1, 1, 1, 1), (1, 3), (3, 1)]
    examples = []
    for w in hostile_words:
        A, Q, B, a, p, d, r = data(w)
        examples.append(dict(word=list(w), source=r, endpoint=replay(r, w),
                             slope=str(a), debt=str(d),
                             forward_cone=d <= a-1, inverse_cone=d <= 1-a))
    require(examples[0]["source"] == 79 and examples[0]["endpoint"] == 101,
            "expansion hostile")
    require(examples[1]["forward_cone"] and not examples[2]["forward_cone"],
            "order hostile")
    require(not examples[3]["inverse_cone"], "mixed-cone hostile")
    require(data((1,))[5] <= data((1,))[3]-1, "forward factor")
    require(data((3,))[5] <= 1-data((3,))[3], "inverse factor")
    require(data((1, 3))[5] < 19-11, "inverse bound is not necessary")
    pairs, root_exclusions = 0, 0
    for n in range(3, 1002, 2):
        u, h = first_descent(n)
        if h == 1:
            root_exclusions += 1
            continue
        v, k = first_descent(h)
        au, pu, du = suffix_data(u)
        av, pv, dv = suffix_data(v)
        require(du < n-h and dv < h-k, "control local payment")
        require(av < 1, "descending slope")
        composed = av*du+dv
        require(composed < n-k, "composed source payment")
        require(composed/F(n-k) <= max(du/F(n-h), dv/F(h-k)), "ratio bound")
        require(replay(n, u+v) == k, "source-preserving concatenation")
        pairs += 1
    print(json.dumps(dict(status="FINITE-EXACT", words=len(words), cuts=cuts,
                          cylinder_replays=cylinder_replays,
                          forward_cone_words=forward,
                          reciprocal_cone_words=reciprocal,
                          paid_descending_pairs=pairs,
                          first_endpoint_root_exclusions=root_exclusions,
                          examples=examples), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
