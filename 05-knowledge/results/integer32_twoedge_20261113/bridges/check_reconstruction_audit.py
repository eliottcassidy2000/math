#!/usr/bin/env python3
"""Independent local encoding/cut/nearest-word oracles against the frozen producer."""
import hashlib
import importlib.util
import itertools as it
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RECON = HERE.parent / "reconstruction"
PINS = {
    "PROOF.md": "f2ba72de4e691ff25f46a1849109a4db30fa0a6917c31a3363b0190b5998647c",
    "INHERITED_AUDIT.md": "089421c85674284d820f12b2b22f116aaee29a05775d0f139e5ff58754ed448a",
    "check.py": "2f6991f470496b2b802fe8a20760fd459562148678b9015be2fab1272cc8ab05",
    "normal.json": "bdd9fbc57328a637c1153cde5bea2c7c6f29d650d5744b6eb1d758055dbfcab1",
    "optimized.json": "bdd9fbc57328a637c1153cde5bea2c7c6f29d650d5744b6eb1d758055dbfcab1",
}
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def main():
    for name, digest in PINS.items():
        need(hashlib.sha256((RECON / name).read_bytes()).hexdigest() == digest, "review pin " + name)
    spec = importlib.util.spec_from_file_location("frozen_reconstruction", RECON / "check.py")
    producer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(producer)
    rows = []
    for n in range(3, 7):
        edges = list(it.combinations(range(n), 2))
        wedges = [(v, u, w) for v in range(n) for u, w in it.combinations([z for z in range(n) if z != v], 2)]
        wi = {w: k for k, w in enumerate(wedges)}
        images = {}
        mindistance = len(wedges)
        base = None
        for x in range(0, 1 << len(edges), 2):
            out = [set() for _ in range(n)]
            deg = [0] * n
            for k, (u, v) in enumerate(edges):
                if x >> k & 1:
                    out[v].add(u)
                    deg[u] += 1
                    deg[v] += 1
                else:
                    out[u].add(v)
            # Independent directed incidence, no producer formula/base subtraction.
            beta = sum(int((u in out[v]) != (w in out[v])) << k for k, (v, u, w) in enumerate(wedges))
            need(beta == producer.encode(n, x), "direct incidence encoding")
            need(beta not in images, "converse quotient injective")
            images[beta] = x
            if base is None:
                base = beta
            cut = sum(d * (n - 1 - d) for d in deg)
            need(cut == (beta ^ base).bit_count(), "independent degree cut")
            if x:
                mindistance = min(mindistance, cut)
        need(mindistance == 2 * n - 4, "independent exact distance")
        # Formal observer fields at n<=4: raw star/triangle equations and exhaustive
        # nearest-word comparison, without using a producer recognition method.
        formal_count = nearest_count = 0
        if n <= 4:
            for beta in range(1 << len(wedges)):
                def bit(v, u, w):
                    return (beta >> wi[v, min(u, w), max(u, w)]) & 1
                stars = all((bit(v, a, b) ^ bit(v, a, c) ^ bit(v, b, c)) == 0
                            for v in range(n) for a, b, c in it.combinations([z for z in range(n) if z != v], 3))
                triangles = all((bit(a, b, c) ^ bit(b, a, c) ^ bit(c, a, b)) == 1
                                for a, b, c in it.combinations(range(n), 3))
                need((stars and triangles) == (beta in images), "independent complete parity image")
                need(producer.decode(n, beta) == images.get(beta), "producer exact decoder vs complete image")
                close = [x for word, x in images.items() if (word ^ beta).bit_count() <= n - 3]
                need(len(close) <= 1, "independent unique radius")
                expected = close[0] if close else None
                need(producer.repair(n, beta) == expected, "repair vs exhaustive nearest word")
                formal_count += 1
                nearest_count += bool(close)
        # n5: all guaranteed reference errors, independently counted/decoded against
        # the entire 512-word image. n6: a declared deterministic 256-pattern prefix.
        controls = 0
        if n >= 5:
            errors = [0]
            for k in range(1, n - 2):
                errors.extend(sum(1 << i for i in s) for s in it.combinations(range(len(wedges)), k))
            if n == 6:
                errors = errors[:256]
            for error in errors:
                observed = base ^ error
                close = [x for word, x in images.items() if (word ^ observed).bit_count() <= n - 3]
                need(close == [0], "full-image nearest control")
                need(producer.repair(n, observed) == 0, "producer repair vs independent oracle")
                controls += 1
        rows.append({"n": n, "image_size": len(images), "minimum_distance": mindistance,
                     "all_formal_fields": formal_count, "formal_fields_in_radius": nearest_count,
                     "reference_error_controls": controls})
    print(json.dumps({"status": "FINITE-EXACT independent audit controls PASS", "reviewed_sha256": PINS,
                      "scope": "independent direct-incidence encoding and degree-cut on all converse states n3..6; complete formal-field/nearest-word checks n3..4; all466 reference errors n5; first256 weight-ordered reference errors n6",
                      "rows": rows, "checks": CHECKS,
                      "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
