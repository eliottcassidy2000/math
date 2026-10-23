"""Independent sparse signed-path and integer-clock controls."""
from collections import Counter
from pathlib import Path
import json


def check(value, message):
    if not value:
        raise RuntimeError(message)


def signed_refinement_rounds(n):
    colors = [0]*n
    rounds = 0
    while len(set(colors)) < n:
        signatures = []
        for v in range(n):
            counts = Counter()
            if v+1 < n:
                counts[colors[v+1]] += 1
            if v:
                counts[colors[v-1]] -= 1
            signatures.append((colors[v], tuple(sorted((c, x) for c, x in counts.items() if x))))
        names = {s: i for i, s in enumerate(sorted(set(signatures)))}
        next_colors = [names[s] for s in signatures]
        check(len(set(next_colors)) > len(set(colors)), "signed refinement progresses")
        colors = next_colors
        rounds += 1
    return rounds


def run():
    rows = []
    for n in range(2, 257):
        rounds = signed_refinement_rounds(n)
        check(rounds == n//2, "signed-path observer exact round count")
        period = 1 << (n-1).bit_length()
        state, mask = 1, (1 << n)-1
        first_return = None
        for t in range(1, period+1):
            state = (state ^ (state << 1)) & mask
            if state == 1:
                first_return = t
                break
        check(first_return == period, "Pascal snapshot exact period")
        rows.append({"n": n, "signed_rounds": rounds, "pascal_period": period})
    higher_modulus_rows = []
    for n in range(2, 33):
        for s in range(1, 5):
            modulus = 1 << s
            period = 1 << (s-1+(n-1).bit_length())
            initial = [1]+[0]*(n-1)
            state = initial[:]
            first_return = None
            for t in range(1, period+1):
                state = [state[0]]+[(state[j]+state[j-1]) % modulus
                                    for j in range(1, n)]
                if state == initial:
                    first_return = t
                    break
            check(first_return == period, "exact power-of-two modulus clock period")
            higher_modulus_rows.append({"n": n, "s": s, "modulus": modulus,
                                        "period": period, "first_return": first_return})
    return {"status": "FINITE-EXACT", "scope": "signed-path observer alone; orders2..256",
            "not_claimed": "sharpness of actual two-tournament refinement at every order",
            "rows": rows, "higher_modulus_scope": "n2..32,s1..4; chronological integer residues",
            "higher_modulus_rows": higher_modulus_rows}


if __name__ == "__main__":
    result = run()
    Path(__file__).with_suffix(".json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print("PASS: 255 signed-path round controls, 255 binary Pascal periods, 124 higher-modulus periods")
