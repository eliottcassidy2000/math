"""Patch the frozen THM-4231 residual by all 45 proved literal boundary pairs."""

import contextlib
import io
import runpy
from collections import defaultdict

CANONICAL = "04-computation/lrc14_all_fixed_outsider_ray_symmetric_postprocess_thm4231.py"

BOUNDARY_LITERAL = {
    (744, 824), (744, 822), (744, 821), (744, 820), (744, 818),
    (744, 817), (744, 815), (744, 814), (744, 813), (744, 812),
    (744, 811), (744, 810), (744, 809), (744, 805), (744, 803),
    (744, 800), (766, 800), (744, 798), (744, 794), (744, 793),
    (744, 791), (744, 790), (744, 789), (744, 787), (744, 780),
    (765, 780), (766, 780), (768, 780), (616, 777), (616, 776),
    (616, 775), (744, 775), (616, 774), (616, 773), (744, 773),
    (616, 772), (616, 771), (721, 771), (616, 770), (721, 770),
    (744, 770), (750, 770), (765, 770), (766, 770), (768, 770),
}

INHERITED_LITERAL = {
    (1, 542),
    (49, 50),
    (50, 51),
}


def fnv_add(value, word):
    for shift in range(0, 64, 8):
        value ^= (word >> shift) & 0xFF
        value = (value * 0x100000001B3) & ((1 << 64) - 1)
    return value


def edge_fnv(edges):
    value = 0xCBF29CE484222325
    for a, b in edges:
        value = fnv_add(value, a)
        value = fnv_add(value, b)
    return value


sink = io.StringIO()
with contextlib.redirect_stdout(sink):
    namespace = runpy.run_path(CANONICAL)

qs = namespace["qs"]
kappa = namespace["kappa"]
original = sorted(
    (a, b)
    for a in qs
    for b in qs
    if a < b and b < kappa[a] and a < kappa[b]
)
assert original == namespace["residual_edges"]
assert len(original) == 181_242
assert edge_fnv(original) == 0x8A4E1370FB023907

original_set = set(original)
assert len(BOUNDARY_LITERAL) == 45
assert BOUNDARY_LITERAL == {edge for edge in original if edge[1] >= 770}
assert INHERITED_LITERAL <= original_set
assert BOUNDARY_LITERAL.isdisjoint(INHERITED_LITERAL)

combined_closures = BOUNDARY_LITERAL | INHERITED_LITERAL
remainder = [edge for edge in original if edge not in combined_closures]
layers = defaultdict(list)
for edge in remainder:
    layers[edge[1]].append(edge)
max_endpoint = max(layers)
next_layer = layers[max_endpoint]

print("LRC14_LITERAL_BOUNDARY_45_CUMULATIVE_PATCHED_RESIDUAL")
print(f"ORIGINAL COUNT {len(original)} FNV {edge_fnv(original):016x}")
print(
    f"BOUNDARY_LITERAL COUNT {len(BOUNDARY_LITERAL)} "
    f"INHERITED_LITERAL COUNT {len(INHERITED_LITERAL)} "
    f"COMBINED_CLOSURES COUNT {len(combined_closures)}"
)
print(
    f"REMAINDER COUNT {len(remainder)} FNV {edge_fnv(remainder):016x} "
    f"MAX_ENDPOINT {max_endpoint} CUTOFF {max_endpoint + 1}"
)
print(
    f"NEXT_LAYER ENDPOINT {max_endpoint} COUNT {len(next_layer)} "
    + "EDGES " + " ".join(f"{a},{b}" for a, b in next_layer)
)
print("BOUNDARY_EQUALS_ALL_ORIGINAL_EDGES_WITH_ENDPOINT_GE_770 YES")
print("INHERITED 1,542 49,50 50,51")
print("VERDICT EXACT_PATCHED_GRAPH LRC14_OPEN")
