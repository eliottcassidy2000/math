"""All seven-sheet K7 covers admitting a branch five-cycle, in a fixed-anchor gauge."""
from itertools import permutations
from construction import compose, disjoint

IDENTITY = tuple(range(7))
ROTATE_BRANCH = (0, 2, 3, 4, 5, 1, 6)
ROTATE_SHEET = (1, 2, 3, 4, 0, 5, 6)
SEEDS = ((6, 1), (1, 2), (1, 3))
PAIR_TYPES = ((0, 6), (0, 1), (6, 1), (1, 2), (1, 3))


def inverse(p):
    return tuple(p.index(x) for x in IDENTITY)


def conjugate(h, p):
    return compose(h, compose(p, inverse(h)))


def orbit_assignment(i, j, p):
    edges = {}
    for _ in range(5):
        edges[i, j] = p
        edges[j, i] = inverse(p)
        i, j = ROTATE_BRANCH[i], ROTATE_BRANCH[j]
        p = conjugate(ROTATE_SHEET, p)
    return edges


def compatible(edges):
    # Simultaneous rotation and path reversal identify all other branch pairs.
    for i, j in PAIR_TYPES:
        routes = []
        if (i, j) in edges:
            routes.append(edges[i, j])
        for k in IDENTITY:
            if (i, k) in edges and (k, j) in edges:
                routes.append(compose(edges[k, j], edges[i, k]))
        if any(len(set(images)) != len(routes) for images in zip(*routes)):
            return False
    return True


def enumerate_covers():
    all_perms = list(permutations(IDENTITY))
    centralizer = [h for h in all_perms if compose(h, ROTATE_SHEET) == compose(ROTATE_SHEET, h)]
    anchor = {(0, i): IDENTITY for i in range(1, 7)}
    anchor.update({(i, 0): IDENTITY for i in range(1, 7)})
    domains = []
    for idx, (i, j) in enumerate(SEEDS):
        options, seen = [], set()
        for p in all_perms:
            if not disjoint(p, IDENTITY) or p in seen:
                continue
            if idx == 0:
                seen.update(conjugate(h, p) for h in centralizer)
            orbit = orbit_assignment(i, j, p)
            if compatible(anchor | orbit):
                options.append(orbit)
        domains.append(options)
    assert [len(d) for d in domains] == [40, 720, 720]
    surviving_pairs = final_tests = full = 0
    for fixed in domains[0]:
        for adjacent in domains[1]:
            partial = anchor | fixed | adjacent
            if compatible(partial):
                surviving_pairs += 1
                for diagonal in domains[2]:
                    final_tests += 1
                    if compatible(partial | diagonal):
                        full += 1
    assert (surviving_pairs, final_tests, full) == (114, 82080, 0)
    print("domains 40,720,720; partial pairs 114; full tuples 82080; solutions 0")


if __name__ == "__main__":
    enumerate_covers()
