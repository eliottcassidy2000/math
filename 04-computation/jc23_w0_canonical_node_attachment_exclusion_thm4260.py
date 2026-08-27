#!/usr/bin/env python3
"""Exact compiler for the final THM-4249/4253 W=0 incidence frontier.

This packet is intentionally outside the repository.  It loads three frozen
repository computations directly from a named git revision, checks their raw
SHA-256 hashes, reconstructs the 176/104 map-class quotient and all 1,512
class--torsion-ratio incidences, and evaluates the hidden involution projector
in exact rational-function arithmetic over two independently chosen good finite
fields.

The proof lane is the q=397 embedding.  For a hidden projection of degree
12*K, the characteristic-zero pole ledger bounds its reduced odd denominator
by 4*K-1.  The specialized denominator attains that bound, so specialization
has lost no denominator degree.  Its reciprocal gcd is exactly t^2-1, the two
iota-fixed points on the excluded Z=0 wall.  Therefore no additional
characteristic-zero common factor exists.  The q=577 embedding is an
independent hostile control.

The native finite carrier is a bipartite graph:

    map class  --  admissible CM-torsion unit orbit (marked ratio).

There is no intrinsic oriented pairwise observable, so no tournament is used.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import subprocess
import sys
import types


DEFAULT_REVISION = "d526261faf8fb33a3b210b19c7d0039b8ef18ca8"

DEPENDENCIES = {
    "projector": (
        "04-computation/jc23_w0_cyclic_projector_squeeze_thm4249.py",
        "64cefef1ab610cdeab05eeaaeff25ae03bb2c69095f734e86d69b92fdccfea10",
    ),
    "torsion": (
        "04-computation/jc23_w0_cyclic_projector_squeeze_thm4249_independent_audit.py",
        "1c1ae0d47f5218af5978cb840c0f6f9c564a6df338a7b650700cbca774e5e3c4",
    ),
    "attachment": (
        "04-computation/jc23_w0_hidden_degree12_attachment_audit_thm4247.py",
        "9dc1f8614db388f463acb93951285fb4b8245cf20c976d1d2129cb883ccf9c28",
    ),
}

PROOF_EMBEDDING = (397, 157, 161, 27)
HOSTILE_EMBEDDINGS = ((577, 57, 224, 25),)
EMBEDDINGS = (PROOF_EMBEDDING,) + HOSTILE_EMBEDDINGS


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def git_bytes(repo: Path, revision: str, path: str) -> bytes:
    result = subprocess.run(
        ["git", "-C", str(repo), "show", f"{revision}:{path}"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return result.stdout


def load_frozen_module(
    repo: Path, revision: str, name: str, path: str, expected_hash: str
) -> types.ModuleType:
    source = git_bytes(repo, revision, path)
    actual_hash = sha256(source).hexdigest()
    require(actual_hash == expected_hash, f"dependency hash changed: {path}")
    module = types.ModuleType(name)
    module.__file__ = f"{repo}@{revision}:{path}"
    sys.modules[name] = module
    exec(compile(source, module.__file__, "exec"), module.__dict__)
    return module


def resolved_revision(repo: Path, revision: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", revision],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.stdout.strip()


def fraction_token(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def torsion_token(orbit: set[tuple[Fraction, Fraction]]) -> tuple[tuple[Fraction, Fraction], ...]:
    return tuple(sorted(orbit))


def torsion_token_text(token: tuple[tuple[Fraction, Fraction], ...]) -> str:
    return ";".join(
        f"{fraction_token(left)},{fraction_token(right)}" for left, right in token
    )


def e_text(value: tuple[int, int]) -> str:
    return f"{value[0]},{value[1]}"


def vector_text(vector: tuple[tuple[int, int], ...]) -> str:
    return ";".join(e_text(value) for value in vector)


def counter_text(counter: Counter) -> str:
    return ",".join(f"{key}:{counter[key]}" for key in sorted(counter))


def hidden_profile(projector, vector):
    _, b_value, c_value, d_value = vector
    hidden_a = projector.e_add(
        projector.e_scale(2, b_value), projector.e_mul(projector.OMEGA2, d_value)
    )
    hidden_b = projector.e_add(projector.e_scale(2, c_value), d_value)
    hidden_degree = projector.hidden_degree(hidden_a, hidden_b)
    require(hidden_degree % 12 == 0, "hidden degree lost divisibility by twelve")
    return hidden_a, hidden_b, projector.e_norm(d_value), hidden_degree // 12


def compile_classes(projector):
    _, _, residual, _, _, _, _ = projector.enumerate_full_residual()
    class_representatives = {}
    profile_counts = {}
    raw_counts = {degree: len(residual[degree]) for degree in (34, 42)}

    for degree in (34, 42):
        vectors = residual[degree]
        if degree == 42:
            norm_three = {
                vector
                for vector in vectors
                if hidden_profile(projector, vector)[2:] == (3, 13)
            }
            require(len(norm_three) == 672, "THM-4253 vector deletion changed")
            deleted_representatives, deleted_sizes = projector.symmetry_orbits(norm_three)
            require(
                len(deleted_representatives) == 28
                and deleted_sizes == Counter({24: 28}),
                "THM-4253 orbit deletion changed",
            )
            vectors = vectors - norm_three

        representatives, sizes = projector.symmetry_orbits(vectors)
        expected_classes = 176 if degree == 34 else 104
        require(
            len(representatives) == expected_classes
            and sizes == Counter({24: expected_classes}),
            f"degree-{degree} refined class quotient changed",
        )
        class_representatives[degree] = representatives
        profile_counts[degree] = Counter(
            hidden_profile(projector, vector)[2:] for vector in representatives
        )

        # On a_u=0, precomposition by rho:(x,y)->(-x,-y) fixes f,g and
        # negates v.  In the glued basis this is
        #   (0,b,c,d) -> (0,b+omega^2*d,c+d,-d).
        # It equals postcomposition by -omega^2 after T^2, so rho fixes each
        # source-target symmetry class individually and exchanges the two
        # radical C12 node orbits without moving to a different map class.
        for vector in representatives:
            a_value, b_value, c_value, d_value = vector
            require(a_value == projector.ZERO, "refined residual regained a u-coordinate")
            rho_vector = (
                projector.ZERO,
                projector.e_add(b_value, projector.e_mul(projector.OMEGA2, d_value)),
                projector.e_add(c_value, d_value),
                projector.e_neg(d_value),
            )
            tau_squared = projector.tau_vector(projector.tau_vector(vector))
            classwise_rho = projector.unit_scale_vector(
                projector.e_neg(projector.OMEGA2), tau_squared
            )
            require(rho_vector == classwise_rho, "rho=(-omega^2)T^2 failed classwise")

    require(raw_counts == {34: 4224, 42: 3168}, "THM-4249 residual changed")
    require(
        profile_counts[34]
        == Counter({
            (4, 10): 36,
            (7, 9): 52,
            (13, 7): 32,
            (16, 6): 24,
            (19, 5): 24,
            (25, 3): 8,
        }),
        "degree-34 class profile changed",
    )
    require(
        profile_counts[42]
        == Counter({(9, 11): 24, (12, 10): 36, (21, 7): 32, (27, 5): 12}),
        "degree-42 refined class profile changed",
    )
    return class_representatives, profile_counts


def compile_incidence_graph(projector, torsion, class_representatives):
    pi = projector.e_sub(projector.OMEGA2, projector.ONE)
    kernel_pi = torsion.torsion_kernel(pi)
    kernel_two = torsion.torsion_kernel((2, 0))
    kernel_three = torsion.torsion_kernel((3, 0))
    excluded_small = kernel_pi | kernel_two
    ratio_one_third_orbit = kernel_three - kernel_pi
    require(
        len(ratio_one_third_orbit) == 6
        and torsion.torsion_orbits(ratio_one_third_orbit)
        == [ratio_one_third_orbit],
        "ratio-one-third torsion orbit changed",
    )

    raw_edges = []
    class_records = []
    tokens_by_degree = {}
    for degree in (34, 42):
        degree_edges = []
        for index, vector in enumerate(class_representatives[degree]):
            hidden_a, hidden_b, norm_d, hidden_k = hidden_profile(projector, vector)
            d_value = vector[3]
            annihilator = projector.e_mul(pi, d_value)
            orbits = [
                orbit
                for orbit in torsion.torsion_orbits(torsion.torsion_kernel(annihilator))
                if orbit.isdisjoint(excluded_small)
            ]
            if degree == 42:
                orbits = [orbit for orbit in orbits if orbit != ratio_one_third_orbit]
            class_id = f"C{degree}_{index:03d}"
            tokens = sorted(torsion_token(orbit) for orbit in orbits)
            class_records.append(
                {
                    "degree": degree,
                    "class_id": class_id,
                    "vector": vector,
                    "hidden_a": hidden_a,
                    "hidden_b": hidden_b,
                    "norm_d": norm_d,
                    "hidden_k": hidden_k,
                    "tokens": tokens,
                }
            )
            degree_edges.extend((class_id, token) for token in tokens)
        raw_edges.extend((degree, class_id, token) for class_id, token in degree_edges)
        tokens_by_degree[degree] = {token for _, token in degree_edges}

    all_tokens = sorted({token for _, _, token in raw_edges})
    ratio_ids = {token: f"R{index:03d}" for index, token in enumerate(all_tokens)}
    edges = sorted(
        (degree, class_id, ratio_ids[token]) for degree, class_id, token in raw_edges
    )
    ratio_rows = sorted((ratio_ids[token], torsion_token_text(token)) for token in all_tokens)

    edge_counts = Counter(degree for degree, _, _ in edges)
    require(edge_counts == Counter({34: 864, 42: 648}), "1,512 edge split changed")
    require(len(edges) == 1512, "incidence total changed")
    require(
        {degree: len(tokens_by_degree[degree]) for degree in (34, 42)}
        == {34: 55, 42: 34},
        "ratio envelope sizes changed",
    )
    require(
        len(tokens_by_degree[34] & tokens_by_degree[42]) == 7
        and len(all_tokens) == 82,
        "cross-degree ratio overlap changed",
    )

    classes_per_ratio = {
        degree: Counter(ratio_id for edge_degree, _, ratio_id in edges if edge_degree == degree)
        for degree in (34, 42)
    }
    require(
        Counter(classes_per_ratio[34].values())
        == Counter({8: 12, 12: 18, 16: 12, 24: 6, 26: 6, 60: 1}),
        "degree-34 ratio incidence multiplicities changed",
    )
    require(
        Counter(classes_per_ratio[42].values()) == Counter({12: 9, 16: 18, 36: 7}),
        "degree-42 ratio incidence multiplicities changed",
    )

    class_rows = []
    for record in sorted(class_records, key=lambda row: row["class_id"]):
        ratio_list = ",".join(ratio_ids[token] for token in record["tokens"])
        class_rows.append(
            "\t".join(
                (
                    record["class_id"],
                    str(record["degree"]),
                    vector_text(record["vector"]),
                    e_text(record["hidden_a"]),
                    e_text(record["hidden_b"]),
                    str(record["norm_d"]),
                    str(record["hidden_k"]),
                    str(len(record["tokens"])),
                    ratio_list,
                )
            )
        )
    edge_rows = [f"{degree}\t{class_id}\t{ratio_id}" for degree, class_id, ratio_id in edges]
    ratio_text_rows = [f"{ratio_id}\t{token}" for ratio_id, token in ratio_rows]

    return {
        "class_records": class_records,
        "class_rows": class_rows,
        "edge_rows": edge_rows,
        "ratio_rows": ratio_text_rows,
        "tokens_by_degree": tokens_by_degree,
        "classes_per_ratio": classes_per_ratio,
        "class_sha256": sha256("\n".join(class_rows).encode("ascii")).hexdigest(),
        "edge_sha256": sha256("\n".join(edge_rows).encode("ascii")).hexdigest(),
        "ratio_sha256": sha256("\n".join(ratio_text_rows).encode("ascii")).hexdigest(),
    }


def audit_reciprocal_denominators(projector, attachment, class_records, embeddings):
    groups = defaultdict(list)
    for record in class_records:
        groups[(record["degree"], record["hidden_k"])].append(
            (record["hidden_a"], record["hidden_b"])
        )

    expected_group_counts = {
        (34, 3): 8,
        (34, 5): 24,
        (34, 6): 24,
        (34, 7): 32,
        (34, 9): 52,
        (34, 10): 36,
        (42, 5): 12,
        (42, 7): 32,
        (42, 10): 36,
        (42, 11): 24,
    }
    require(
        {key: len(value) for key, value in groups.items()} == expected_group_counts,
        "hidden projector group counts changed",
    )

    rows = []
    for embedding in embeddings:
        prime, zeta12, rho, scale = embedding
        for degree, hidden_k in sorted(groups):
            representatives = groups[(degree, hidden_k)]
            expected_degree = 4 * hidden_k - 1
            denominator_degrees, gcd_degrees, gcd_polynomials, digest = (
                attachment.finite_field_attachment_audit(
                    prime,
                    zeta12,
                    rho,
                    scale,
                    representatives,
                    expected_degree,
                )
            )
            require(
                Counter(denominator_degrees)
                == Counter({expected_degree: len(representatives)}),
                "specialized denominator did not attain the pole bound",
            )
            require(
                Counter(gcd_degrees) == Counter({2: len(representatives)}),
                "reciprocal gcd gained a non-wall factor",
            )
            expected_gcd = (1, 0, prime - 1)
            require(
                Counter(gcd_polynomials) == Counter({expected_gcd: len(representatives)}),
                "reciprocal gcd is not the monic wall t^2-1",
            )
            rows.append(
                "\t".join(
                    (
                        str(prime),
                        str(degree),
                        str(hidden_k),
                        str(len(representatives)),
                        str(expected_degree),
                        "2",
                        "1,0,-1",
                        digest,
                    )
                )
            )
    return rows, sha256("\n".join(rows).encode("ascii")).hexdigest()


def print_report(
    revision: str,
    profile_counts,
    graph,
    audit_rows,
    audit_sha256: str,
    proof_only: bool,
    full_ledger: bool,
) -> None:
    degree34_class_edges = Counter(
        len(record["tokens"])
        for record in graph["class_records"]
        if record["degree"] == 34
    )
    degree42_class_edges = Counter(
        len(record["tokens"])
        for record in graph["class_records"]
        if record["degree"] == 42
    )
    degree34_ratio_degrees = Counter(graph["classes_per_ratio"][34].values())
    degree42_ratio_degrees = Counter(graph["classes_per_ratio"][42].values())

    print("JC W=0 FINAL 1512 CANONICAL-NODE COMPILER")
    print("status=FINITE_EXACT_PROOF_CANDIDATE_RELATIVE_TO_THM4247_4249_4253")
    print(f"revision={revision}")
    for name in sorted(DEPENDENCIES):
        path, digest = DEPENDENCIES[name]
        print(f"dependency_{name}={path} sha256={digest}")
    print("carrier=bipartite_class_ratio_incidence_graph tournament=NOT_USED")
    print("orientation_gauge=NONE ties=NOT_APPLICABLE")
    print("node_gauge=one_C12_orbit_of_24_radical_choices rho_shift_(3,2)_exchanges_other_orbit")
    print("node_class_firewall=rho=(-omega^2)T^2_on_a_u0 fixes_each_of_280_classes PASS")
    print("projector=ell=m-m_after_iota coefficients=(2b+omega^2*d,2c+d)")
    print("projector_eigenline=integral_hidden_P_L(m)=ell collapse_implies_ell(Q_j)=O")
    print("denominator_form=d_ell(t)=t*D(t^2) char0_degree_bound=4K-1")
    print("wall_factor=t^2-1 iota_fixed_points=t:+-1 marked_wall=Z/U:0")
    print(f"class_profiles_d34={dict(sorted(profile_counts[34].items()))}")
    print(f"class_profiles_d42={dict(sorted(profile_counts[42].items()))}")
    print("classes=d34:176 d42:104 total:280 all_source_target_orbits_size24")
    print("incidences=d34:864 d42:648 total:1512")
    print("ratio_vertices=d34:55 d42:34 cross_degree_overlap:7 union:82")
    print(f"class_incidence_degrees_d34={counter_text(degree34_class_edges)}")
    print(f"class_incidence_degrees_d42={counter_text(degree42_class_edges)}")
    print(f"ratio_class_degrees_d34={counter_text(degree34_ratio_degrees)}")
    print(f"ratio_class_degrees_d42={counter_text(degree42_ratio_degrees)}")
    print(f"class_ledger_sha256={graph['class_sha256']}")
    print(f"edge_ledger_sha256={graph['edge_sha256']}")
    print(f"ratio_token_ledger_sha256={graph['ratio_sha256']}")
    print(f"embedding_mode={'proof_only_q397' if proof_only else 'proof_plus_hostile'}")
    for row in audit_rows:
        print(f"denominator_audit={row}")
    print(f"denominator_audit_sha256={audit_sha256}")
    print("proof_embedding=q397 every_denominator_attains_4K-1 reciprocal_gcd=t^2-1")
    print("hostile_embeddings=q577 same_exact_result" if not proof_only else "hostile_embeddings=SKIPPED_BY_FLAG")
    print("consequence=all_280_hidden_projectors_fail_off_wall hence_all_1512_incidences_excluded")
    print("candidate_theorem=S34_empty_and_S42_empty_on_W0_gate")
    print("scope=relative_W0_gate_only other_hidden_Hom_parameters_M12_entry_JC2_DC2_OPEN")
    print("checks=PASS exact_integer_lattice exact_torsion_module exact_GF(t)_group_law no_floating_point")

    if full_ledger:
        print("BEGIN_RATIO_TOKENS")
        print("ratio_id\tcanonical_mu6_orbit_in_O_tensor_Q/O")
        print("\n".join(graph["ratio_rows"]))
        print("END_RATIO_TOKENS")
        print("BEGIN_CLASS_LEDGER")
        print("class_id\tdegree\tvector_u_f_g_h\thidden_A\thidden_B\tN(d)\tK\tincidence_count\tratio_ids")
        print("\n".join(graph["class_rows"]))
        print("END_CLASS_LEDGER")
        print("BEGIN_EDGE_LEDGER")
        print("degree\tclass_id\tratio_id")
        print("\n".join(graph["edge_rows"]))
        print("END_EDGE_LEDGER")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--revision", default=DEFAULT_REVISION)
    parser.add_argument(
        "--proof-only",
        action="store_true",
        help="run only the logically sufficient q=397 embedding",
    )
    parser.add_argument(
        "--full-ledger",
        action="store_true",
        help="print all 82 ratio tokens, 280 classes, and 1,512 edges",
    )
    arguments = parser.parse_args()
    repo = arguments.repo.resolve()
    revision = resolved_revision(repo, arguments.revision)

    loaded = {}
    for name, (path, expected_hash) in DEPENDENCIES.items():
        loaded[name] = load_frozen_module(
            repo, revision, f"jc1512_{name}", path, expected_hash
        )

    class_representatives, profile_counts = compile_classes(loaded["projector"])
    graph = compile_incidence_graph(
        loaded["projector"], loaded["torsion"], class_representatives
    )
    embeddings = (PROOF_EMBEDDING,) if arguments.proof_only else EMBEDDINGS
    audit_rows, audit_hash = audit_reciprocal_denominators(
        loaded["projector"], loaded["attachment"], graph["class_records"], embeddings
    )
    print_report(
        revision,
        profile_counts,
        graph,
        audit_rows,
        audit_hash,
        arguments.proof_only,
        arguments.full_ledger,
    )


if __name__ == "__main__":
    main()
