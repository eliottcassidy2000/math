"""Exact companion for THM-3227's residue contact trie and first carry."""

from itertools import product
import ast
import hashlib
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = ROOT / (
    "01-canon/theorems/"
    "THM-3221-selected-root-osculating-separation-and-minimal-jet-prime-carry.md"
)
DEPENDENCY_SHA256 = "6ae707482cf73ebbd995722fe056d093df37d3aea0af6f58d6989e18d0b06f84"


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(payload).hexdigest()


require(lf_sha256(DEPENDENCY) == DEPENDENCY_SHA256, "THM-3221 dependency hash")
syntax_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax_tree))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax_tree)
)
require(assert_nodes == 0, "assert statements are optimization-sensitive")
require(float_literals == 0, "floating literals are forbidden")


def mul(left, right, degree):
    out = [0] * (degree + 1)
    for i, a in enumerate(left[: degree + 1]):
        for j, b in enumerate(right[: degree + 1 - i]):
            out[i + j] += a * b
    return out


def compose(outer, inner, degree):
    out = [0] * (degree + 1)
    power = [0] * (degree + 1)
    power[0] = 1
    for order, coefficient in enumerate(outer[: degree + 1]):
        for j in range(degree + 1):
            out[j] += coefficient * power[j]
        if order < degree:
            power = mul(power, inner, degree)
    return out


def compose_mod(outer, inner, degree, p):
    return [value % p for value in compose(outer, inner, degree)]


def inverse_tangent_mod(poly, degree, p):
    inverse = [0] * (degree + 1)
    inverse[1] = 1
    for order in range(2, degree + 1):
        error = compose_mod(poly, inverse, order, p)[order]
        inverse[order] = (-error) % p
    identity = [0] * (degree + 1)
    identity[1] = 1
    require(compose_mod(poly, inverse, degree, p) == identity, "left inverse")
    require(compose_mod(inverse, poly, degree, p) == identity, "right inverse")
    return inverse


def contact(word_a, word_b):
    for index, (a, b) in enumerate(zip(word_a, word_b), start=2):
        if a != b:
            return index
    return 10**9


# Full coefficient-word tries are sharp and ultrametric.
ultrametric_checks = 0
full_trie_checks = 0
for alphabet_size, length in ((2, 4), (3, 3), (4, 2), (5, 2), (9, 1)):
    words = list(product(range(alphabet_size), repeat=length))
    require(len(words) == alphabet_size**length, "full residue word bank")
    # Every prefix node has at most q children, and the full bank attains q.
    for prefix_length in range(length):
        prefixes = {word[:prefix_length] for word in words}
        for prefix in prefixes:
            children = {
                word[prefix_length]
                for word in words
                if word[:prefix_length] == prefix
            }
            require(len(children) == alphabet_size, "sharp q-ary branching")
            full_trie_checks += 1
    # Sum(children-1) over internal nodes equals leaves-1.
    label_count = sum(
        alphabet_size**prefix_length * (alphabet_size - 1)
        for prefix_length in range(length)
    )
    require(label_count == len(words) - 1, "minimal affine-label census")
    for a in words:
        for b in words:
            nu_ab = contact(a, b)
            for c in words:
                nu_bc = contact(b, c)
                nu_ac = contact(a, c)
                require(nu_ac >= min(nu_ab, nu_bc), "contact ultrametric")
                if nu_ab != nu_bc:
                    require(nu_ac == min(nu_ab, nu_bc), "strict ultrametric")
                ultrametric_checks += 1


# Extension-field branching is q, while every nonzero additive digit has order p.
extension_additive_order_checks = 0
for p, dimension in ((2, 2), (2, 3), (3, 2)):
    zero = (0,) * dimension
    for value in product(range(p), repeat=dimension):
        if value == zero:
            continue
        running = zero
        for step in range(1, p + 1):
            running = tuple((a + b) % p for a, b in zip(running, value))
            if step < p:
                require(running != zero, "premature extension-field reset")
        require(running == zero, "extension-field additive order p")
        extension_additive_order_checks += 1


# Nonlinear integral source charts preserve reduced contact and its weight.
coordinate_checks = 0
for p in (3, 5, 7):
    degree = 6
    for rho in range(1, p):
        inverse_rho = pow(rho, -1, p)
        inverse_chart = [0] * (degree + 1)
        inverse_chart[1] = inverse_rho
        for order in range(2, degree + 1):
            inverse_chart[order] = (rho + 2 * order) % p
        for m in range(2, degree + 1):
            base = [0] * (degree + 1)
            base[1] = 1
            for order in range(2, m):
                base[order] = (order * rho + 1) % p
            left = list(base)
            right = list(base)
            left[m] = 0
            right[m] = 1
            left_new = [
                rho * value % p
                for value in compose_mod(left, inverse_chart, degree, p)
            ]
            right_new = [
                rho * value % p
                for value in compose_mod(right, inverse_chart, degree, p)
            ]
            require(
                contact(left[2:], right[2:]) == m,
                "original contact order",
            )
            require(
                contact(left_new[2:], right_new[2:]) == m,
                "new contact order",
            )
            require(
                (right_new[m] - left_new[m]) % p == pow(rho, 1 - m, p),
                "weighted leading difference",
            )
            coordinate_checks += 1


# Earlier p-divisible jets do not contaminate the first residue-visible carry.
primitive_carry_checks = 0
for p in (2, 3, 5, 7):
    for m in range(2, 9):
        for seed in range(3):
            for unit in range(1, p):
                germ = [0] * (m + 1)
                germ[1] = 1
                for order in range(2, m):
                    germ[order] = p * (seed + order)
                germ[m] = unit + p * (seed + 1)
                iterate = [0] * (m + 1)
                iterate[1] = 1
                for _ in range(p):
                    iterate = compose(germ, iterate, m)
                difference = list(iterate)
                difference[1] -= 1
                require(
                    all(value % p == 0 for value in difference),
                    "residue reset",
                )
                require(
                    (difference[m] // p) % p == germ[m] % p,
                    "primitive divided carry",
                )
                primitive_carry_checks += 1


# Ramified smallness is insufficient: a Nottingham tail can precede /p.
ramified_boundary_checks = 0
for lower in range(1, 6):
    for visible in range(3):
        germ = [0, 1, lower, 0, 0, visible]
        iterate = [0, 1, 0, 0, 0, 0]
        for _ in range(3):
            iterate = compose(germ, iterate, 5)
        require(
            iterate[5] == 3 * visible + 10 * lower**4,
            "ramified Nottingham coefficient",
        )
        ramified_boundary_checks += 1


# Reduction is not a graph contraction: it can postpone and recreate splits.
delayed = (
    (0, 0, 0),
    (3, 0, 1),
    (6, 1, 0),
)
require(all(contact(delayed[i], delayed[j]) == 2 for i in range(3) for j in range(i)), "generic depth-two split")
delayed_mod3 = tuple(tuple(value % 3 for value in word) for word in delayed)
require(contact(delayed_mod3[0], delayed_mod3[2]) == 3, "delayed depth-three split")
require(contact(delayed_mod3[0], delayed_mod3[1]) == 4, "delayed depth-four split")


print("dependency_thm3221_sha256=%s" % DEPENDENCY_SHA256)
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("ultrametric_triple_checks=%d" % ultrametric_checks)
print("sharp_full_trie_node_checks=%d" % full_trie_checks)
print("affine_label_census=leaves_minus_one")
print("extension_additive_order_checks=%d" % extension_additive_order_checks)
print("nonlinear_coordinate_checks=%d" % coordinate_checks)
print("primitive_carry_checks=%d" % primitive_carry_checks)
print("ramified_boundary_checks=%d" % ramified_boundary_checks)
print("branch_bound=residue_field_size")
print("carry_order=residue_characteristic")
print("delayed_resplitting_hostile=generic_depth2_to_residue_depth3_and4")
print("ramified_hostile=p3_e5:u+pi*u2+u5_has_u5_return_3+10*pi4")
print("scope=selected_simple_root_residue_trie_not_global_branch_selector")
print("all_exact_checks=PASS")
