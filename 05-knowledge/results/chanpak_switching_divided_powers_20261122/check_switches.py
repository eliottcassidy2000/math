"""Compute Boolean differences of directed path and cycle counts in one switching orbit."""

from collections import Counter
from itertools import combinations

ORDER = 9
CODE = 0x7103D59DC


def predecessor_masks(code, cut):
    incoming = [0] * ORDER
    for i, (a, b) in enumerate(combinations(range(ORDER), 2)):
        forward = ((code >> i) ^ (cut >> a) ^ (cut >> b)) & 1
        u, v = (a, b) if forward else (b, a)
        incoming[v] |= 1 << u
    return incoming


def count_paths_cycles(incoming):
    size = 1 << ORDER
    paths = [[0] * ORDER for _ in range(size)]
    anchored = [[0] * ORDER for _ in range(size)]
    for v in range(ORDER):
        paths[1 << v][v] = 1
    anchored[1][0] = 1
    for subset in range(1, size):
        if subset.bit_count() == 1:
            continue
        vertices = subset
        while vertices:
            bit = vertices & -vertices
            vertices ^= bit
            end = bit.bit_length() - 1
            previous = subset ^ bit
            choices = previous & incoming[end]
            path_count = cycle_prefixes = 0
            while choices:
                last = choices & -choices
                choices ^= last
                v = last.bit_length() - 1
                path_count += paths[previous][v]
                cycle_prefixes += anchored[previous][v]
            paths[subset][end] = path_count
            anchored[subset][end] = cycle_prefixes
    # Every directed spanning cycle has one ordering starting at vertex 0.
    return sum(paths[-1]), sum(anchored[-1][v] for v in range(ORDER)
                               if incoming[0] & (1 << v))


def boolean_transform(values, sign):
    result = values.copy()
    for i in range(ORDER - 1):
        for s in range(len(result)):
            if s & (1 << i):
                result[s] += sign * result[s ^ (1 << i)]
    return result


def show_checks(name, values, shift, modulus, radius):
    differences = boolean_transform(values, -1)
    for s, delta in enumerate(differences):
        exponent = max(0, (s.bit_count() + 1) // 2 - shift)
        assert delta % (1 << exponent) == 0

    # Values outside the sample ball are placeholders.
    samples = [value if s.bit_count() <= radius else 0 for s, value in enumerate(values)]
    coefficients = boolean_transform(samples, -1)
    truncated = [delta if s.bit_count() <= radius else 0 for s, delta in enumerate(coefficients)]
    reconstructed = boolean_transform(truncated, 1)
    assert all((a-b) % modulus == 0 for a, b in zip(reconstructed, values))

    nonzero = Counter(s.bit_count() for s, delta in enumerate(differences)
                      if delta % modulus)
    print(f"{name}: base={values[0]}, mod {modulus} nonzero differences by degree "
          f"{[nonzero[j] for j in range(ORDER)]}")
    for j in range(ORDER - 1, 0, -1):
        witness = next(((s, delta) for s, delta in enumerate(differences)
                        if s.bit_count() == j and delta % modulus), None)
        if witness:
            s, delta = witness
            labels = [i+1 for i in range(ORDER-1) if s & (1 << i)]
            recovered = reconstructed[s]
            print(f"  highest visible degree {j}, switched labels {labels}: "
                  f"difference {delta}, value {values[s]}, reconstructed {recovered % modulus}")
            break


def main():
    code = CODE
    print(f"order {ORDER}, lexicographic orientation code {code:#x}, vertex 0 fixed")
    # A cut and its complement produce the same tournament, so fix spin 0.
    counts = [count_paths_cycles(predecessor_masks(code, mask << 1))
              for mask in range(1 << (ORDER - 1))]
    h, hc = map(list, zip(*counts))
    show_checks("H", h, 0, 24, 4)
    show_checks("H", h, 0, 12, 2)
    show_checks("HC", hc, 1, 4, 4)
    assert {value % 24 for value in h} == set(range(1, 24, 2))
    for residue in (7, 21):
        print(f"least count in residue {residue} mod 24: "
              f"{min(value for value in h if value % 24 == residue)}")
    joint = boolean_transform([value - 6 * cycles for value, cycles in counts], -1)
    assert all(delta == 0 for s, delta in enumerate(joint) if s.bit_count() > 4)
    print("all 256 cuts recovered mod 24 from 163 cuts, mod 12 from 37 cuts")
    print("the degree-five and higher differences of H - 6 HC vanish over the integers")


if __name__ == "__main__":
    main()
