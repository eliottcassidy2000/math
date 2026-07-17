#!/usr/bin/env python3
"""Exact referee for THM-933, the sharp local-density block-gluing theorem.

The theorem separates two jobs.

* Inside a block B, retain the safe density delta(B) and the exact centered
  primitive discrepancy q(B).
* Across blocks, pay q(B_r) once for each component made by earlier teeth.

All interval arithmetic below uses Fraction.  No sampling is used.

Tournament Analysis uses certified blocks as vertices.  For two blocks A,B,
the pairwise observable is the second-block boundary debt

    debt(A then B) = q(B) M(A).

The gauge orients A -> B when that debt is no larger than debt(B then A),
with increasing physical scale as the tie Hamiltonian path.  This quotient
preserves the two-block density product and boundary debt.  It destroys the
within-block endpoint order, which is retained exactly in the q sidecar.
"""

from collections import Counter
from bisect import bisect_right
from fractions import Fraction
from itertools import combinations, permutations


F = Fraction


def measure(intervals):
    return sum((b - a for a, b in intervals), F(0))


def merge_intervals(intervals):
    merged = []
    for a, b in sorted(intervals):
        if not a < b:
            continue
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return merged


def complement_intervals(intervals):
    out = []
    cursor = F(0)
    for a, b in merge_intervals(intervals):
        if cursor < a:
            out.append((cursor, a))
        cursor = max(cursor, b)
    if cursor < 1:
        out.append((cursor, F(1)))
    return out


def danger_intervals(speed, denominator):
    """The split lift to [0,1] of ||speed*t|| < 1/denominator."""
    assert speed >= 1 and denominator >= 3
    half_width = F(1, denominator * speed)
    intervals = []
    for residue in range(speed):
        left = F(residue, speed) - half_width
        right = F(residue, speed) + half_width
        if left < 0:
            intervals.extend([(F(0), right), (left + 1, F(1))])
        elif right > 1:
            intervals.extend([(left, F(1)), (F(0), right - 1)])
        else:
            intervals.append((left, right))
    return merge_intervals(intervals)


def safe_intervals(speeds, denominator):
    danger = []
    for speed in speeds:
        danger.extend(danger_intervals(speed, denominator))
    return complement_intervals(danger)


def primitive_data(safe):
    """Return delta, q=osc primitive, and all breakpoint primitive values."""
    density = measure(safe)
    points = sorted({F(0), F(1), *(endpoint for interval in safe for endpoint in interval)})
    primitive = F(0)
    values = [(F(0), primitive)]
    safe_index = 0
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        while safe_index < len(safe) and safe[safe_index][1] <= midpoint:
            safe_index += 1
        inside = (
            safe_index < len(safe)
            and safe[safe_index][0] < midpoint < safe[safe_index][1]
        )
        slope = 1 - density if inside else -density
        primitive += slope * (right - left)
        values.append((right, primitive))
    assert primitive == 0
    primitive_values = [value for _, value in values]
    discrepancy = max(primitive_values) - min(primitive_values)
    return density, discrepancy, values


class SafePrefix:
    """Exact cumulative safe measure with logarithmic arc queries."""

    def __init__(self, safe):
        self.starts = [left for left, _ in safe]
        self.ends = [right for _, right in safe]
        self.cumulative = [F(0)]
        for left, right in safe:
            self.cumulative.append(self.cumulative[-1] + right - left)
        self.total = self.cumulative[-1]

    def cdf(self, point):
        index = bisect_right(self.starts, point) - 1
        if index < 0:
            return F(0)
        clipped = min(max(point - self.starts[index], F(0)),
                      self.ends[index] - self.starts[index])
        return self.cumulative[index] + clipped

    def arc(self, start, length):
        assert 0 <= start < 1 and 0 <= length <= 1
        end = start + length
        if end <= 1:
            return self.cdf(end) - self.cdf(start)
        return self.total - self.cdf(start) + self.cdf(end - 1)


def eta_at_length(safe, length):
    """Exact fixed-scale local density floor eta(length)."""
    assert 0 < length <= 1
    candidates = {F(0)}
    for interval in safe:
        for endpoint in interval:
            candidates.add(endpoint % 1)
            candidates.add((endpoint - length) % 1)
    prefix = SafePrefix(safe)
    retained = min(prefix.arc(start, length) for start in candidates)
    return retained / length


def primitive_deficit_extremizer(certificate):
    """An arc attaining q = length * (delta - eta(length))."""
    values = certificate["values"]
    maximum = max(value for _, value in values)
    minimum = min(value for _, value in values)
    candidates = []
    for start, value_start in values:
        if value_start != maximum:
            continue
        for end, value_end in values:
            if value_end != minimum:
                continue
            length = (end - start) % 1
            if length == 0:
                length = F(1)
            candidates.append((length, start, end))
    length, start, end = min(candidates)
    eta = eta_at_length(certificate["safe"], length)
    assert certificate["q"] == length * (certificate["density"] - eta)
    return length, eta, start, end


def block_certificate(speeds, denominator, label):
    safe = safe_intervals(speeds, denominator)
    density, discrepancy, values = primitive_data(safe)
    return {
        "label": label,
        "speeds": tuple(speeds),
        "denominator": denominator,
        "safe": safe,
        "density": density,
        "q": discrepancy,
        "M": sum(speeds),
        "values": values,
    }


def verify_endpoint_extremizers(certificate):
    """Exhaust the finite endpoint-pair certificate for a modest block."""
    values = [value for _, value in certificate["values"]]
    minimum = min(v - u for u in values for v in values)
    maximum = max(v - u for u in values for v in values)
    assert minimum == -certificate["q"]
    assert maximum == certificate["q"]
    return len(values) ** 2


def scaled_certificate(template, scale, denominator, label, materialize=True):
    base = block_certificate(template, denominator, f"{label}-template")
    result = {
        "label": label,
        "speeds": tuple(scale * speed for speed in template),
        "denominator": denominator,
        "density": base["density"],
        "q": base["q"] / scale,
        "M": scale * base["M"],
        "template": base,
    }
    if materialize:
        direct = block_certificate(result["speeds"], denominator, label)
        assert direct["density"] == result["density"]
        assert direct["q"] == result["q"]
        assert direct["M"] == result["M"]
        result.update({"safe": direct["safe"], "values": direct["values"]})
    return result


def product(values):
    answer = F(1)
    for value in values:
        answer *= value
    return answer


def gluing_bound(certificates):
    density_product = product(certificate["density"] for certificate in certificates)
    debt = F(0)
    earlier_teeth = 0
    terms = []
    for index, certificate in enumerate(certificates):
        if index:
            suffix = product(
                later["density"] for later in certificates[index + 1 :]
            )
            term = certificate["q"] * earlier_teeth * suffix
            debt += term
            terms.append((certificate["label"], earlier_teeth, suffix, term))
        earlier_teeth += certificate["M"]
    return density_product - debt, density_product, debt, terms


def gluing_bound_with_component_caps(certificates, component_caps):
    """Closed ledger using certified prefix-component caps instead of tooth sums."""
    assert len(component_caps) == len(certificates) - 1
    density_product = product(certificate["density"] for certificate in certificates)
    debt = F(0)
    terms = []
    for index, certificate in enumerate(certificates[1:], start=1):
        suffix = product(later["density"] for later in certificates[index + 1 :])
        term = certificate["q"] * component_caps[index - 1] * suffix
        debt += term
        terms.append((certificate["label"], component_caps[index - 1], suffix, term))
    return density_product - debt, density_product, debt, terms


def verify_block_recurrence(blocks, certificates, denominator):
    prior_speeds = []
    previous_measure = F(1)
    rows = []
    for index, (block, certificate) in enumerate(zip(blocks, certificates)):
        if index == 0:
            component_count = 1
            step_floor = certificate["density"]
        else:
            previous_safe = safe_intervals(prior_speeds, denominator)
            component_count = len(previous_safe)
            assert component_count <= sum(prior_speeds)
            step_floor = (
                certificate["density"] * previous_measure
                - certificate["q"] * component_count
            )
        prior_speeds.extend(block)
        current_safe = safe_intervals(prior_speeds, denominator)
        current_measure = measure(current_safe)
        assert current_measure >= step_floor
        rows.append((index + 1, component_count, sum(prior_speeds), current_measure, step_floor))
        previous_measure = current_measure
    return rows, previous_measure


def tournament_fingerprint(certificates):
    labels = [certificate["label"] for certificate in certificates]
    scale_key = {
        certificate["label"]: min(certificate["speeds"])
        for certificate in certificates
    }

    def debt(first, second):
        return second["q"] * first["M"]

    def edge(first, second):
        forward = debt(first, second)
        backward = debt(second, first)
        if forward < backward:
            return first["label"], second["label"]
        if backward < forward:
            return second["label"], first["label"]
        ordered = sorted(
            [first["label"], second["label"]], key=lambda x: (scale_key[x], x)
        )
        return ordered[0], ordered[1]

    directed = {}
    pair_rows = []
    for first, second in combinations(certificates, 2):
        source, target = edge(first, second)
        directed[(source, target)] = True
        pair_rows.append(
            (
                first["label"],
                second["label"],
                debt(first, second),
                debt(second, first),
                source,
            )
        )

    scores = {label: 0 for label in labels}
    for source, _ in directed:
        scores[source] += 1
    score_histogram = Counter(scores.values())

    cycles = 0
    for triple in combinations(labels, 3):
        out = {label: 0 for label in triple}
        for first, second in combinations(triple, 2):
            if (first, second) in directed:
                out[first] += 1
            else:
                out[second] += 1
        cycles += sorted(out.values()) == [1, 1, 1]

    reach = {a: {a} for a in labels}
    for source, target in directed:
        reach[source].add(target)
    changed = True
    while changed:
        changed = False
        for label in labels:
            expanded = set().union(*(reach[item] for item in tuple(reach[label])))
            if not expanded <= reach[label]:
                reach[label] |= expanded
                changed = True
    unused = set(labels)
    scc_sizes = []
    while unused:
        seed = next(iter(unused))
        component = {label for label in unused if seed in reach[label] and label in reach[seed]}
        scc_sizes.append(len(component))
        unused -= component

    paths = []
    for ordering in permutations(labels):
        if all((ordering[i], ordering[i + 1]) in directed for i in range(len(ordering) - 1)):
            paths.append(ordering)

    density_edges = set()
    for first, second in combinations(certificates, 2):
        if first["density"] > second["density"]:
            density_edges.add((first["label"], second["label"]))
        elif second["density"] > first["density"]:
            density_edges.add((second["label"], first["label"]))
        else:
            density_edges.add(tuple(sorted((first["label"], second["label"]))))
    edge_flips = sum(
        ((a, b) in directed) != ((a, b) in density_edges)
        for a, b in combinations(labels, 2)
    )

    all_bounds = []
    lookup = {certificate["label"]: certificate for certificate in certificates}
    for ordering in permutations(labels):
        bound, _, _, _ = gluing_bound([lookup[label] for label in ordering])
        all_bounds.append((bound, ordering))
    all_bounds.sort(reverse=True)

    return {
        "pairs": pair_rows,
        "scores": scores,
        "score_histogram": dict(sorted(score_histogram.items())),
        "cycles": cycles,
        "scc_sizes": sorted(scc_sizes, reverse=True),
        "paths": paths,
        "edge_flips": edge_flips,
        "best_order": all_bounds[0],
        "all_bounds": all_bounds,
    }


print("=" * 78)
print("THM-933 EXACT REFEREE: SHARP LOCAL-DENSITY BLOCK GLUING")
print("=" * 78)

print("\n[1] Primitive discrepancy and scale covariance")
templates = (
    ("A", [1, 2, 3, 4, 5], 1),
    ("B", [1, 2, 4, 5], 30),
    ("C", [1, 3, 5, 7], 900),
)
certificates = []
blocks = []
for label, template, scale in templates:
    certificate = scaled_certificate(template, scale, 14, label)
    certificates.append(certificate)
    blocks.append(list(certificate["speeds"]))
    endpoint_checks = verify_endpoint_extremizers(certificate["template"])
    print(
        f"  {label}: template={template}, scale={scale}, "
        f"delta={certificate['density']} ({float(certificate['density']):.9f}), "
        f"q={certificate['q']} ({float(certificate['q']):.9f}), "
        f"M={certificate['M']}; endpoint pairs={endpoint_checks}; "
        "direct scale check=PASS"
    )

print("\n[2] Exact non-lacunary 13-speed block composition")
all_speeds = [speed for block in blocks for speed in block]
assert len(all_speeds) == 13 and len(set(all_speeds)) == 13
rows, actual = verify_block_recurrence(blocks, certificates, 14)
bound, density_product, debt, terms = gluing_bound(certificates)
component_caps = [row[1] for row in rows[1:]]
component_bound, _, component_debt, component_terms = (
    gluing_bound_with_component_caps(certificates, component_caps)
)
assert actual >= bound > 0
assert actual >= component_bound > bound
for stage, components, prefix_sum, current, step_floor in rows:
    print(
        f"  stage {stage}: previous components={components}, "
        f"prefix tooth sum={prefix_sum}, actual={float(current):.9f}, "
        f"one-step floor={float(step_floor):.9f}"
    )
print(f"  speeds = {all_speeds}")
print(f"  product density = {density_product} ({float(density_product):.9f})")
for label, earlier_teeth, suffix, term in terms:
    print(
        f"  debt at {label}: q*M_previous*suffix = {term} "
        f"({float(term):.9f}); M_previous={earlier_teeth}, suffix={suffix}"
    )
print(f"  THM-933 lower bound = {bound} ({float(bound):.9f})")
print(
    f"  exact-component strengthening = {component_bound} "
    f"({float(component_bound):.9f}); caps={component_caps}, debt={component_debt}"
)
print(f"  exact common-safe measure = {actual} ({float(actual):.9f})")
print("  verdict: POSITIVE; not covered by a runner-by-runner 7-lacunary hypothesis")

print("\n[3] Singleton corollaries (exact constants)")
for speed in (1, 2, 7, 19):
    singleton = block_certificate([speed], 14, f"x={speed}")
    assert singleton["density"] == F(6, 7)
    assert singleton["q"] == F(6, 49 * speed)
print("  delta({x})=6/7 and q({x})=6/(49x): PASS at x=1,2,7,19")
a = F(6, 7)
r7_bound = a**13 - a * (1 - a**12) / 6
assert r7_bound == F(199455593, 13841287201) > 0
assert 6**12 > 7**11
print(
    f"  R>=7: mu > {r7_bound} ({float(r7_bound):.9f}); "
    f"integer check 6^12 - 7^11 = {6**12 - 7**11} > 0"
)
uniform_checks = []
for runner_count in range(3, 501):
    safe_fraction = F(runner_count - 2, runner_count)
    uniform_checks.append(8 * safe_fraction ** (runner_count - 2) > 1)
assert all(uniform_checks)
print("  uniform R>=8 criterion: exact PASS for n=3..500; analytic proof uses e^2<8")

print("\n[4] THM-928(C) certificate receives its missing local sidecar")
packet_928 = [
    300, 406, 511, 652, 862, 963, 1074,
    1357, 1459, 1571, 1776, 1991, 2208,
]
packet_certificate = block_certificate(packet_928, 13, "THM-928(C)")
bonf5 = F(624500321285438209432944191, 15959187221015235325636692600)
assert 0 < bonf5 <= packet_certificate["density"]
print(f"  stronger radius = 1/13 (hence also an LRC(14) certificate)")
print(f"  exact BONF5 lower density d = {bonf5} ({float(bonf5):.9f})")
print(
    f"  exact safe density delta = {packet_certificate['density']} "
    f"({float(packet_certificate['density']):.9f})"
)
print(
    f"  exact primitive q sidecar = {packet_certificate['q']} "
    f"({float(packet_certificate['q']):.9f})"
)
dual_length, dual_eta, dual_start, dual_end = primitive_deficit_extremizer(
    packet_certificate
)
probe_length = F(4, 300)
probe_eta = eta_at_length(packet_certificate["safe"], probe_length)
assert packet_certificate["q"] >= probe_length * (
    packet_certificate["density"] - probe_eta
)
print(
    f"  eta/q dual extremizer: ell={dual_length} ({float(dual_length):.9f}), "
    f"eta(ell)={dual_eta} ({float(dual_eta):.9f}), start={dual_start}, end={dual_end}"
)
print(
    f"  exact duality q=ell*(delta-eta(ell)): PASS; "
    f"Opus probe eta(4/300)={probe_eta} ({float(probe_eta):.9f})"
)
print(
    f"  safe components={len(packet_certificate['safe'])}, "
    f"tooth cap M={packet_certificate['M']}; d<=delta and scale law make it glueable"
)

print("\n[5] Pulled Opus-S333 fixed-scale 7+6 composition")
opus_first_speeds = [300, 406, 511, 652, 862, 963, 1074]
opus_second_template = [862, 963, 1074, 1357, 1459, 1571]
opus_second_scale = 300
opus_first = block_certificate(opus_first_speeds, 13, "Opus-B7")
opus_second = block_certificate(opus_second_template, 13, "Opus-B6-template")
opus_template_length = F(2, 300)
opus_actual_length = opus_template_length / opus_second_scale
opus_eta = eta_at_length(opus_second["safe"], opus_template_length)
opus_components = len(opus_first["safe"])
opus_sharp_bound = opus_eta * (
    opus_first["density"] - opus_components * opus_actual_length
)
opus_weaker_bound = (
    opus_first["density"] * opus_eta - opus_components * opus_actual_length
)
assert 0 <= opus_eta <= 1
assert opus_sharp_bound > opus_weaker_bound > 0
print(
    f"  B1 density={opus_first['density']} ({float(opus_first['density']):.9f}), "
    f"components={opus_components}"
)
print(
    f"  eta_B2(1/45000)=eta_template(2/300)={opus_eta} "
    f"({float(opus_eta):.9f})"
)
print(
    f"  sharp G1 bound eta*(mu-kappa*ell)={opus_sharp_bound} "
    f"({float(opus_sharp_bound):.9f})"
)
print(
    f"  weaker Opus ledger mu*eta-kappa*ell={opus_weaker_bound} "
    f"({float(opus_weaker_bound):.9f}); both certify 13 non-lacunary speeds"
)

print("\n[6] Tournament Analysis on block vertices")
fingerprint = tournament_fingerprint(certificates)
print("  observable: debt(A then B)=q(B)M(A); gauge chooses smaller debt")
for first, second, forward, backward, source in fingerprint["pairs"]:
    print(
        f"  {first}/{second}: forward={float(forward):.9f}, "
        f"reverse={float(backward):.9f}; edge source={source}"
    )
print(f"  scores={fingerprint['scores']}; histogram={fingerprint['score_histogram']}")
print(
    f"  directed 3-cycles={fingerprint['cycles']}; SCC sizes={fingerprint['scc_sizes']}; "
    f"Hamiltonian paths={len(fingerprint['paths'])}: {fingerprint['paths']}"
)
print(f"  edge flips versus density-only gauge={fingerprint['edge_flips']}")
best_bound, best_order = fingerprint["best_order"]
print(f"  best full-ledger order={best_order}, bound={float(best_bound):.9f}")
assert best_order == ("A", "B", "C") and best_bound == bound

print("\n[7] Referee verdict")
print("  primitive interval lemma: EXACT (finite endpoint extrema)")
print("  scale covariance: EXACT (direct recomputation on all three blocks)")
print("  component cap and every recurrence step: PASS")
print("  coarse and exact-component closed ledgers: PASS")
print("  primitive and fixed-scale positive 13-speed witnesses: PASS")
print("  lacunary improvement R>=7 and uniform R>=8: PASS")
print("  assumption challenge: blocks preserve (delta,q,M), not endpoint order")
print("  SCOPE: proves the composition theorem and new families, not arbitrary LRC(14) alone")
