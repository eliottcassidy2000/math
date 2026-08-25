#!/usr/bin/env python3
"""Primary exact certificate and optional full replay for THM-4137.

The default mode verifies the frozen 16-shard certificate, recomputes the
small exact controls, and emits the canonical transcript in under a second.
Pass ``--full`` to regenerate both nauty streams, compile the pinned exact
contracted-block kernel, and replay all 9,355,949 strong order-ten classes.
"""

from argparse import ArgumentParser
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from hashlib import sha256
import importlib.util
from json import dumps
import os
from pathlib import Path
import re
import shutil
import subprocess
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
BASE_PY = ROOT / "04-computation/tournament_strong_centrality_through_order_eight_thm4131.py"
BASE_CPP = ROOT / "04-computation/tournament_strong_centrality_order_nine_thm4135_independent_audit.cpp"
FAMILY_PY = ROOT / "04-computation/tournament_strong_cyclic_substitution_centrality_counterexample_thm4133.py"

EXPECTED_BASE_PY_SHA256 = "6b195b0379d1ae3e5d215aa1c495f7180daeecae189df86269d07ef855867881"
EXPECTED_BASE_CPP_SHA256 = "5ed81bd2ab3054d1a05471c42fe724c14eb0df278dbc3ee3b7d6d38253fff530"
EXPECTED_FAMILY_PY_SHA256 = "7bd4c518464d4c48baf9cb9c1c8c2012a9f79f029c8e07141c0e51c338ffeeb2"
EXPECTED_TRANSFORMED_CPP_SHA256 = "e860bdb634ddf93dc01120be8218a139d31e70a1d8658e32f94047197d0bf4bd"
EXPECTED_GENTOURNG_SHA256 = "89df605922cc574b28688248b7c256d24342cc615f887e89b2d096038970c110"
EXPECTED_ALL_RAW_SHA256 = "af2873154068897522bc15477d989b0577877d2bbebc08aea3082353e5378b67"
EXPECTED_STRONG_RAW_SHA256 = "47bcaa3ef6272261dee3092735b47b3d2154d882aae6a8420c964cd3ef7289b7"
EXPECTED_CERTIFICATE_SHA256 = "c57ff382a892fd7136714605ed73940ef7bb0dba39af4cf5e98ff816ee5ffa6c"
EXPECTED_SEMANTIC_SHA256 = "39a9adba9b2c4d3917e07ef24df62fc02d9203c53e32b8d0de7f02b61e66c114"


# Each row is
# (res, classes, raw_sha256, exact_profile_sha256, rational_fail,
#  coset_fail, central_minus_outer_margin, worst_rho_num, worst_rho_den,
#  worst_multiplicity, theta_zero, coset_reorders, actual_noncentral,
#  actual_optimizer_histogram_from_the_kernel).
SHARDS = (
    (0, 462493, "bcb1efb3f578e33592de01d07576fe1202065e6835c3c7f6fd766caa36db0078", "132c566d500ee12adb181f4ff7bb022e337e9f760c51decd73294ccec44a594f", 0, 0, 24, 36301, 39843, 1, 344, 0, 153615, "(-2):62039,(-2,-4):3,(-4):923,(0):308142,(0,-2):346,(2):89183,(2,0):390,(4):1454,(4,2):13"),
    (1, 511673, "f4e0761c7521c38e6d67b516f4f18134288b15992d83a17141c7a63fd75af770", "40138cecc422ac8134f9525cab3a72fa630e9ff39088a2a71e8951ed767a49a2", 0, 0, 34, 100564247, 109818111, 1, 381, 0, 188467, "(-2):106176,(-2,-4):38,(-4):1899,(0):321957,(0,-2):730,(2):78807,(2,0):519,(4):1536,(4,2):11"),
    (2, 378196, "e48a31182a9942158936026db5c9ab9f005e2735e3c4294a936040ef740525db", "f43fa8e6891b43a9c926006ed660e4152de0aa61116b285b1068c37611e89c71", 0, 0, 38, 69596401, 78319475, 1, 401, 0, 130354, "(-2):65140,(-2,-4):15,(-4):1097,(0):246461,(0,-2):433,(2):63148,(2,0):948,(4):945,(4,2):9"),
    (3, 504700, "082339aa1b2f7f6c0958922cd12d7868cdd7adf60200ddfe7bcae106c488a3d9", "57399b129334d5e6a584f05b8317ca7adbb07298d07b8f6db9da096f1498cf64", 0, 0, 26, 1528856, 1716745, 1, 520, 0, 168611, "(-2):94098,(-2,-4):9,(-4):1604,(0):334609,(0,-2):1023,(2):71183,(2,0):457,(4):1688,(4,2):29"),
    (4, 633698, "cff3a7c69b54c779063f4cbc3620c2f30a4c9832bebb0211eb77e2c1e3f51c91", "eb1823826aa9d25c67e710eae527e359aa73f9742f46910864fb17a1912c39c9", 0, 0, 46, 41340474, 47480165, 1, 606, 0, 206005, "(-2):114394,(-2,-4):32,(-4):1782,(0):425777,(0,-2):1309,(2):87865,(2,0):607,(4):1863,(4,2):69"),
    (5, 591812, "5d19fd75887f87e7d668021fb86658cdd4ab9c7b5eff1e0fc43227fde5e71ccf", "0f777df9bbb7daf6735ec290e750eb5cd355238f35953da71c3b0b8cb46a8900", 0, 0, 36, 40775000, 45270237, 1, 550, 0, 192493, "(-2):104254,(-2,-4):12,(-4):1548,(0):397991,(0,-2):558,(2):85982,(2,0):770,(4):671,(4,2):26"),
    (6, 599797, "997d90279c5854d58388486029e8cc1c96357f8d8c2f0d1fb11dae99c3bda66d", "0217fd52fc654f9324f95a2e609a712d98bd7fb5f878e9723d3ec8f46ec6893e", 0, 0, 46, 64258222, 74708353, 1, 512, 0, 193525, "(-2):73919,(-2,-4):8,(-4):727,(0):405287,(0,-2):354,(2):117729,(2,0):631,(4):1132,(4,2):10"),
    (7, 700320, "1e56fbb2e118f7d5477db1b4c44c08114dd90093a635ee2d0e0ef4efd7ebebed", "e02b512b12d257d629cd5527c5ff69448a74d83abce20d069c0fed13835f3efc", 0, 0, 58, 40735345, 48518483, 1, 590, 0, 220586, "(-2):100408,(-2,-4):9,(-4):1153,(0):477829,(0,-2):774,(2):117915,(2,0):1131,(4):1090,(4,2):11"),
    (8, 566706, "eef8c0c2098bdbd35f63e07f59b1384e823655f9558d273a154b7ec26b0c4d74", "0ba5b3746cadbfe5568ebc72fbe4dd2bf5bb826aca97a708760fa946d763c92a", 0, 0, 24, 36301, 39843, 1, 522, 0, 190670, "(-2):68339,(-2,-4):7,(-4):695,(0):374844,(0,-2):538,(2):119807,(2,0):654,(4):1810,(4,2):12"),
    (9, 521118, "7f618c88f567ecf97b278365ce466c34143ecf6ec61803aad96b53648ab9ab3e", "828190c3ab64a1374d1bb590981fa961439ca4e4969ec3bcddb146206a5b6085", 0, 0, 34, 100575517, 109819011, 1, 651, 0, 197338, "(-2):78937,(-2,-4):44,(-4):1104,(0):322997,(0,-2):248,(2):114272,(2,0):535,(4):2960,(4,2):21"),
    (10, 555062, "e7568fba7e50812dfd794ebbd6dbd42d43f442b86ea033cf00ba81c0a4a5ef40", "2418427be3b638d10d9dfd90b9e3900ec0fd222a43ec62b74eb4659c24c8483c", 0, 0, 30, 34499248, 37325237, 1, 493, 0, 182844, "(-2):83589,(-2,-4):5,(-4):1165,(0):370378,(0,-2):1268,(2):96343,(2,0):572,(4):1720,(4,2):22"),
    (11, 525107, "97af7bb95425dd74f99965bd933a0539203f6939134577e15fc2857ca307531c", "a7f9b394ff8143454e0b1962e78a1aa8772a82a328dbd9d4f43f21ce90aade55", 0, 0, 50, 77823956, 90270449, 1, 509, 0, 169609, "(-2):92582,(-2,-4):22,(-4):1364,(0):353961,(0,-2):688,(2):74384,(2,0):849,(4):1247,(4,2):10"),
    (12, 787448, "549fd1c4c1829a89f71ffdceda1672cf72a261c18fb3a358a5e2283e71e01d90", "de8d387c7bd06c3a5dc3edabaea09e4c7b3739670c30108c52a1da27c6e8f10b", 0, 0, 36, 97613103, 107670631, 1, 695, 0, 267442, "(-2):158774,(-2,-4):49,(-4):2306,(0):518397,(0,-2):785,(2):105087,(2,0):824,(4):1195,(4,2):31"),
    (13, 702125, "b8e6027f86726702b47ff4f9f2419684d61cae66fd5ec30f5faae8296533c176", "fb92c8de94fe073c306efbb9d95277929de66f71660c5030be774974b759a77a", 0, 0, 52, 120431808, 139918237, 1, 651, 0, 231542, "(-2):127721,(-2,-4):9,(-4):1711,(0):468961,(0,-2):609,(2):100860,(2,0):1013,(4):1222,(4,2):19"),
    (14, 660406, "115107cd34c149d38de4b6b62bcfa14453ef0d460a630316b967c4c9b76542dd", "0f0b5c7749e9387fd56f59bd60dd3458625d6e60f18716457fc3ace5d24f1ecd", 0, 0, 30, 34499248, 37325237, 1, 514, 0, 225271, "(-2):109085,(-2,-4):49,(-4):1571,(0):433610,(0,-2):797,(2):113514,(2,0):728,(4):1042,(4,2):10"),
    (15, 655288, "7c2689f37ff83d7bdb4678d9a5ef3b2d6596bbf34ae014cf8491d2c64fed0598", "4e728e8d86c57e0d7c1d42ca02e641cd5d19cc2f25ac91fba0136685ca10affc", 0, 0, 60, 17819725, 21271501, 1, 660, 0, 228600, "(-2):111357,(-2,-4):10,(-4):1704,(0):425432,(0,-2):712,(2):114733,(2,0):544,(4):778,(4,2):18"),
)


WORST_PACKETS = (
    {
        "label": "111111101110011111110111110111111111111111101",
        "H_W_D4_Chd": (1431, 19557, 111975711, -29570784),
        "theta_rho": ((-34499248, 37325237), (34499248, 37325237)),
        "optimizers": ((0,), (0,), (-4,)),
        "layers": ((1, 8, (13234965, 2173), 6093, 6933, 6), (2, 6, (40554497, 4346), 9333, 11805, 2), (3, 4, (536384935, 45633), 11755, 16047, 2), (4, 2, (406414943, 30422), 13361, 19839, 2), (5, 0, (215178219, 15211), 14147, 21923, 2), (6, -2, (1288243325, 91266), 14117, 23377, 2), (7, -4, (201794477, 15211), 13267, 24201, 2), (8, -6, (50411425, 4346), 11601, 23175, 2), (9, -8, (59418751, 6519), 9115, 13041, 2)),
    },
    {
        "label": "101111110111111111111111110011111011101111111",
        "H_W_D4_Chd": (1431, 19557, 111975711, 29570784),
        "theta_rho": ((34499248, 37325237), (34499248, 37325237)),
        "optimizers": ((0,), (0,), (4,)),
        "layers": ((1, 8, (59418751, 6519), 9115, 13041, 2), (2, 6, (50411425, 4346), 11601, 23175, 2), (3, 4, (201794477, 15211), 13267, 24201, 2), (4, 2, (1288243325, 91266), 14117, 23377, 2), (5, 0, (215178219, 15211), 14147, 21923, 2), (6, -2, (406414943, 30422), 13361, 19839, 2), (7, -4, (536384935, 45633), 11755, 16047, 2), (8, -6, (40554497, 4346), 9333, 11805, 2), (9, -8, (13234965, 2173), 6093, 6933, 6)),
    },
)


MIN_MARGIN_PACKETS = (
    ("111111101110011111101111111111111111111111101", (765, 13095, 50202180, -13068360), (-36301, 39843), (0,), (0,), (-2,), 24),
    ("101111110111111111111111111111110011101111111", (765, 13095, 50202180, 13068360), (36301, 39843), (0,), (0,), (2,), 24),
)


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def digest(value):
    return sha256(dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def file_digest(path):
    answer = sha256()
    with open(path, "rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            answer.update(block)
    return answer.hexdigest()


def parse_histogram(text):
    answer = Counter()
    for key, count in re.findall(r"(\([^)]*\)):(\d+)", text):
        entries = tuple(sorted(int(value) for value in re.findall(r"-?\d+", key)))
        answer[entries] += int(count)
    return answer


def shard_dict(row):
    (residue, classes, raw_hash, profile_hash, rational_fail, coset_fail,
     margin, rho_num, rho_den, worst_count, theta_zero, reorder,
     actual_noncentral, actual_histogram) = row
    return {
        "residue_mod_16": residue,
        "classes": classes,
        "raw_sha256": raw_hash,
        "exact_profile_sha256": profile_hash,
        "failures": (rational_fail, coset_fail),
        "minimum_coset_margin": margin,
        "worst_rho": (rho_num, rho_den),
        "worst_multiplicity": worst_count,
        "theta_zero": theta_zero,
        "coset_reorders": reorder,
        "actual_noncentral": actual_noncentral,
        "actual_histogram": tuple(sorted(parse_histogram(actual_histogram).items())),
    }


def load_module(name, path, expected_hash):
    require(file_digest(path) == expected_hash, f"pinned dependency hash: {path.name}")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exact_controls():
    base = load_module("thm4131_base_for_4137", BASE_PY, EXPECTED_BASE_PY_SHA256)
    controls = base.named_controls()
    require(controls["code2_n6"]["rational_t"] == (2,), "rational hostile")
    require(controls["code140_n6"]["coset_t"] == (4,), "coset hostile")
    require(controls["code20_n6"]["rational_t"] == (0,)
            and controls["code20_n6"]["coset_t"] == (0,)
            and controls["code20_n6"]["actual_t"] == (2,), "actual-max hostile")
    family = load_module("thm4133_family_for_4137", FAMILY_PY, EXPECTED_FAMILY_PY_SHA256)
    near = family.compact_family_row(7)
    require(near["order"] == 10 and near["rho"] == (97613103, 107670631),
            "structured order-ten near-boundary")
    require(near["rational_t"] == near["coset_t"] == (0,)
            and near["actual_t"] == (-4,)
            and near["coset_central_minus_outer"] == 36,
            "structured control optimizer row")
    return {
        "code2_n6": (controls["code2_n6"]["rational_t"], controls["code2_n6"]["coset_t"]),
        "code140_n6": (controls["code140_n6"]["rational_t"], controls["code140_n6"]["coset_t"]),
        "code20_n6": (controls["code20_n6"]["rational_t"], controls["code20_n6"]["coset_t"], controls["code20_n6"]["actual_t"]),
        "thm4133_order10": (near["packet"], near["rho"], near["rational_t"], near["coset_t"], near["actual_t"], near["coset_central_minus_outer"]),
    }


def certificate_summary(rows):
    require(tuple(row[0] for row in rows) == tuple(range(16)), "all shard residues")
    parsed = tuple(shard_dict(row) for row in rows)
    require(sum(row["classes"] for row in parsed) == 9355949, "strong universe count")
    require(all(row["failures"] == (0, 0) for row in parsed), "central support-floor failures")
    require(all(row["coset_reorders"] == 0 for row in parsed), "rational/coset optimizer equality")
    require(min(row["minimum_coset_margin"] for row in parsed) == 24, "strict coset margin")
    require(sum(row["theta_zero"] for row in parsed) == 8599, "zero-tilt count")
    require(sum(row["actual_noncentral"] for row in parsed) == 3146972,
            "actual-max scope boundary")
    actual = Counter()
    for row in parsed:
        actual.update(dict(row["actual_histogram"]))
    expected_actual = {
        (-4,): 22353, (-4, -2): 321, (-2,): 1550812,
        (-2, 0): 11172, (0,): 6186633, (0, 2): 11172,
        (2,): 1550812, (2, 4): 321, (4,): 22353,
    }
    require(dict(actual) == expected_actual, "actual optimizer histogram")
    worst = max(Fraction(*row["worst_rho"]) for row in parsed)
    multiplicity = sum(
        row["worst_multiplicity"] for row in parsed
        if Fraction(*row["worst_rho"]) == worst
    )
    require(worst == Fraction(34499248, 37325237) and multiplicity == 2,
            "sharp reversal-paired tilt")
    return {
        "shards": parsed,
        "class_counts": (9733056, 9355949),
        "failures": (0, 0),
        "rational_histogram": (((0,), 9355949),),
        "coset_histogram": (((0,), 9355949),),
        "coset_reorders": 0,
        "minimum_strict_coset_margin": 24,
        "theta_zero": 8599,
        "worst_rho": (34499248, 37325237),
        "worst_multiplicity": 2,
        "actual_noncentral_only": 3146972,
        "actual_histogram": tuple(sorted(expected_actual.items())),
    }


def transformed_kernel():
    require(file_digest(BASE_CPP) == EXPECTED_BASE_CPP_SHA256, "pinned C++ kernel hash")
    source = BASE_CPP.read_text()
    replacements = (
        ("constexpr int N = 9;", "constexpr int N = 10;"),
        ("rat((N - 3) * std::llabs(p.Chdx4), 4 * p.D4x4)",
         "rat((N - 3) * std::llabs(p.Chdx4), 2 * p.D4x4)"),
        ("const uint16_t central_mask = uint16_t((1u << 4) | (1u << 5)); // t=+1,-1",
         "const uint16_t central_mask = uint16_t(1u << 5); // t=0"),
        ("if (x.t == -1 || x.t == 1)", "if (x.t == 0)"),
    )
    for old, new in replacements:
        require(source.count(old) == 1, f"unique kernel transform: {old}")
        source = source.replace(old, new)
    require(sha256(source.encode()).hexdigest() == EXPECTED_TRANSFORMED_CPP_SHA256,
            "transformed order-ten kernel hash")
    return source


def find_program(variable, candidates):
    if os.environ.get(variable):
        return os.environ[variable]
    for candidate in candidates:
        found = shutil.which(candidate)
        if found:
            return found
    raise RuntimeError(f"missing executable for {variable}: {candidates}")


def compile_kernel(source, directory):
    compiler = find_program("THM4137_CXX", ("clang++", "g++"))
    source_path = Path(directory) / "thm4137_n10.cpp"
    binary_path = Path(directory) / "thm4137_n10"
    source_path.write_text(source)
    command = [
        compiler, "-O3", "-std=c++17", "-Xpreprocessor", "-fopenmp",
        "-I/opt/homebrew/opt/libomp/include", "-I/opt/homebrew/opt/openssl/include",
        str(source_path), "-L/opt/homebrew/opt/libomp/lib", "-lomp",
        "-L/opt/homebrew/opt/openssl/lib", "-lcrypto", "-o", str(binary_path),
    ]
    completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    require(completed.returncode == 0,
            "kernel compilation failed: " + completed.stderr.decode(errors="replace"))
    return str(binary_path)


def stream_count_digest(executable, strong):
    command = [executable, "-q"] + (["-c"] if strong else []) + ["10"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    hasher = sha256()
    count = 0
    while True:
        block = process.stdout.read(1 << 20)
        if not block:
            break
        hasher.update(block)
        count += block.count(b"\n")
    stderr = process.stderr.read().decode(errors="replace")
    require(process.wait() == 0, "gentourng stream failed: " + stderr)
    return count, hasher.hexdigest()


def parse_kernel_output(residue, output):
    fields = {}
    for line in output.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            fields.setdefault(key, []).append(value)
    needed = (
        "classes", "strong_classes", "raw_stream_sha256",
        "ordered_full_profile_sha256", "rational_central_failures",
        "coset_central_failures", "minimum_strict_coset_margin", "worst_rho",
        "worst_multiplicity", "theta_zero", "coset_reorders_rational",
        "actual_noncentral_only", "actual_histogram", "status",
    )
    require(all(key in fields for key in needed), f"complete shard output {residue}")
    rho_num, rho_den = map(int, fields["worst_rho"][0].split("/"))
    return (
        residue, int(fields["classes"][0]), fields["raw_stream_sha256"][0],
        fields["ordered_full_profile_sha256"][0],
        int(fields["rational_central_failures"][0]),
        int(fields["coset_central_failures"][0]),
        int(fields["minimum_strict_coset_margin"][0]), rho_num, rho_den,
        int(fields["worst_multiplicity"][0]), int(fields["theta_zero"][0]),
        int(fields["coset_reorders_rational"][0]),
        int(fields["actual_noncentral_only"][0]), fields["actual_histogram"][0],
    )


def run_shard(executable, kernel, residue):
    generator = subprocess.Popen(
        [executable, "-q", "-c", "10", f"{residue}/16"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    audit = subprocess.Popen(
        [kernel], stdin=generator.stdout, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, env=environment,
    )
    generator.stdout.close()
    output, audit_stderr = audit.communicate()
    generator_stderr = generator.stderr.read()
    generator_status = generator.wait()
    require(generator_status == 0,
            f"gentourng shard {residue}: " + generator_stderr.decode(errors="replace"))
    require(audit.returncode == 0,
            f"kernel shard {residue}: " + audit_stderr.decode(errors="replace"))
    return parse_kernel_output(residue, output.decode())


def full_replay(workers):
    gentourng = find_program("THM4137_GENTOURNG", ("gentourng", "nauty-gentourng"))
    require(file_digest(gentourng) == EXPECTED_GENTOURNG_SHA256,
            "audited gentourng binary hash")
    all_stream = stream_count_digest(gentourng, False)
    strong_stream = stream_count_digest(gentourng, True)
    require(all_stream == (9733056, EXPECTED_ALL_RAW_SHA256), "all-class raw stream")
    require(strong_stream == (9355949, EXPECTED_STRONG_RAW_SHA256),
            "strong-class raw stream")
    with TemporaryDirectory(prefix="thm4137-") as directory:
        kernel = compile_kernel(transformed_kernel(), directory)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            rows = tuple(pool.map(
                lambda residue: run_shard(gentourng, kernel, residue), range(16)
            ))
    require(rows == SHARDS, "full replay differs from frozen shard certificate")
    return rows


def main():
    parser = ArgumentParser()
    parser.add_argument("--full", action="store_true",
                        help="replay all 9,355,949 strong classes")
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    arguments = parser.parse_args()
    require(1 <= arguments.workers <= 16, "workers in 1..16")

    transformed_kernel()  # Fast mode also pins the exact load-bearing source.
    rows = full_replay(arguments.workers) if arguments.full else SHARDS
    summary = certificate_summary(rows)
    controls = exact_controls()
    certificate_hash = digest(tuple(shard_dict(row) for row in rows))
    if EXPECTED_CERTIFICATE_SHA256 is not None:
        require(certificate_hash == EXPECTED_CERTIFICATE_SHA256,
                "frozen aggregate certificate digest")
    ledger = {
        "theorem": "THM-4137",
        "universe": "all strong tournament isomorphism classes of order ten",
        "generator": {
            "version": "Homebrew nauty 2.9.3",
            "binary_sha256": EXPECTED_GENTOURNG_SHA256,
            "all_command": "gentourng -q 10",
            "strong_command": "gentourng -q -c 10",
            "shard_command": "gentourng -q -c 10 res/16, res=0..15",
            "all_raw_sha256": EXPECTED_ALL_RAW_SHA256,
            "strong_raw_sha256": EXPECTED_STRONG_RAW_SHA256,
        },
        "kernel": {
            "base_sha256": EXPECTED_BASE_CPP_SHA256,
            "transformed_sha256": EXPECTED_TRANSFORMED_CPP_SHA256,
            "carrier": "contracted good/reversed block DP + exact full cut/layer replay",
            "profile_compatibility_header": "thm4131-n9-independent-profile-v1; payload is transformed N=10 (45-bit label, 9 layers, 1024 responses)",
        },
        "summary": summary,
        "worst_packets": WORST_PACKETS,
        "minimum_margin_packets": MIN_MARGIN_PACKETS,
        "controls": controls,
        "certificate_sha256": certificate_hash,
        "scope": (
            "FINITE-EXACT order ten rational and exact-coset support-floor centrality; "
            "actual maxima remain separate; THM-4133 refutes all-order centrality at twelve"
        ),
    }
    semantic = digest(ledger)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, "frozen semantic digest")

    print("status=PASS")
    print("theorem=THM-4137 complete strong order-ten centrality")
    print("verification=fast frozen 16-shard certificate; --full regenerates all streams and exact profiles")
    print("generator=Homebrew nauty 2.9.3; gentourng -q -c 10 res/16 for res=0..15")
    print(f"generator_binary_sha256={EXPECTED_GENTOURNG_SHA256}")
    print(f"class_counts=all:{summary['class_counts'][0]};strong:{summary['class_counts'][1]}")
    print(f"raw_sha256=all:{EXPECTED_ALL_RAW_SHA256};strong:{EXPECTED_STRONG_RAW_SHA256}")
    print(f"shard_certificate_sha256={certificate_hash}")
    print(f"shard_rows={tuple((row[0], row[1], row[2], row[3]) for row in rows)}")
    print("failures=rational:0;coset:0")
    print("rational_histogram=((0,),9355949);coset_histogram=((0,),9355949);coset_reorders=0")
    print("minimum_strict_coset_margin=24;minimum_margin_multiplicity=2")
    print(f"minimum_margin_packets={MIN_MARGIN_PACKETS}")
    print("worst_rho=34499248/37325237;worst_multiplicity=2")
    print(f"worst_packets={WORST_PACKETS}")
    print(f"theta_zero={summary['theta_zero']}")
    print(f"actual_noncentral_only={summary['actual_noncentral_only']}")
    print(f"actual_histogram={summary['actual_histogram']}")
    print(f"controls={controls}")
    print("replay=python3 -B 04-computation/tournament_strong_centrality_order_ten_thm4137.py --full --workers 8")
    print("scope=FINITE-EXACT n=10 support floors; actual maxima separate; n=11 open; all-order refuted at n=12")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
