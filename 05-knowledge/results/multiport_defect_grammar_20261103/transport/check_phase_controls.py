"""Small readout audit only; does not alter the stable forest checker/proof."""
from pathlib import Path
import argparse
import hashlib
import json
from check_transport import EndpointEvaluator, literal_birth, endpoint_dp, require

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent


def phi343(v):
    h, x, y = v
    return [(h + 4*x + 6*y) % 7, 2*x % 7, (3*x + 4*y) % 7]


def chart(golden):
    h, x, y = golden
    v = [a % 7 for a in golden]
    a, b, c = (x + 4*(h-y)) % 7, y % 7, (h-y) % 7
    require([(b+c) % 7, (a-4*c) % 7, b] == v, "phase inverse")
    orbit = [v, phi343(v), phi343(phi343(v))]
    require(phi343(orbit[-1]) == v, "cubic codec closure")
    return {"observer_mod7": v, "F343_vector_codec": orbit,
            "F49_charge": [a, b, c], "F49_frobenius": [(a+4*b) % 7, -b % 7, c]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE / "phase_controls.json")
    args = parser.parse_args()
    hostile_path = RESULTS / "return6325_completion_20261103/halving_response.json"
    phase_path = RESULTS / "zero_ready_strong_switch_20261103/history_phase/validation_new.json"
    hostile = json.loads(hostile_path.read_text())["hostile"]
    hostile_rows = []
    for index, rows in enumerate(hostile["rows"]):
        e = EndpointEvaluator(rows)
        h, f, s = e.h[-1], e.F(8), e.s[8][-1]
        response = h + f - s
        direct = endpoint_dp(literal_birth(rows, {8}, {8}, False))[0][-1]
        require(h == 462915 and response == direct == hostile["responses"][index], "literal cut256")
        hostile_rows.append({"rows": rows, "H": h, "E8_2": f, "start8": s,
                             "cut": 256, "response": response})
    saved = json.loads(phase_path.read_text())
    records, cache = [], {}
    for state in saved["states"]:
        if state["H"] != 170807:
            continue
        key = tuple(state["rows"])
        if key not in cache:
            cache[key] = EndpointEvaluator(key)
        e = cache[key]
        root = state["names"].index(state["root_name"])
        end = e.e[root][-1]
        s2 = sum(e.s[root][m]*e.h[e.full ^ m] for m in range(e.full+1))
        require((s2-end) % 2 == 0, "integral golden readout")
        golden = [e.h[-1], (s2-end)//2, (s2+end)//2]
        actual = chart(golden)
        require(golden == state["golden"], "literal golden source")
        for field, value in actual.items():
            require(value == state["charts"][field], "literal field chart " + field)
        joints = []
        for other_name in (6, 50, 39):
            if other_name == state["root_name"] or other_name not in state["ports"]:
                continue
            other = state["names"].index(other_name)
            p, q = state["ports"].index(state["root_name"]), state["ports"].index(other_name)
            pair = [e.J(root, other), e.J(other, root)]
            require(pair == [state["all_ports"]["ordered_W2"][p][q],
                             state["all_ports"]["ordered_W2"][q][p]], "retained joint ports")
            joints.append({"ports": [state["root_name"], other_name], "ordered_pair": pair})
        records.append({"tag": state["tag"], "root_name": state["root_name"],
                        "H": e.h[-1], "end_root": end, "S_root_2": s2,
                        "golden": golden, "charts": actual, "joint_controls": joints})
    require(len(records) == 2, "both actual source gauges")
    result = {"status": "FINITE-EXACT readout compatibility; no physical Frobenius transport inferred",
              "hostile_81_27": hostile_rows, "literal_170807_gauges": records,
              "source_sha256": {str(p.relative_to(RESULTS)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in (hostile_path, phase_path)},
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS: literal cut256 controls; both170807 root39/root50 integer and F49/F343 decoders; retained joint ports")


if __name__ == "__main__":
    main()
