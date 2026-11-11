"""Independent literal audit of the seven-persistent-vertex141 reset.

No alternative-lane constructor, DP, attachment, or isomorphism imports.
"""
from pathlib import Path
import hashlib
import json

from check_five_cores import check, paths_dp, ready_word, relabel, strong
from replay141 import attach_out


def old_source(k):
    rows = ready_word(["11"]*k+["01","11","10"])
    return relabel(rows,tuple(range(3,k+6))+(1,2,0))


def induced(rows, labels):
    return tuple(sum(((rows[v] >> w) & 1) << j for j,w in enumerate(labels))
                 for v in labels)


def main():
    reports = []
    for k in range(1,7):
        source = old_source(k)
        retained = [v for v in range(len(source)) if v != k+4]
        deleted = induced(source,retained)
        out = 7 << (k+1)
        first = attach_out(deleted,out)
        reset = attach_out(first,out)
        birth1 = attach_out(reset,1 << (k+2))
        birth2 = attach_out(birth1,1 << len(reset))
        graphs = [source,deleted,first,reset,birth1,birth2]
        counts = [paths_dp(rows) for rows in graphs]
        d = 1 << k
        expected_head = [34*d+5,7*d+5,15*d+9,31*d+17]
        check(counts[:4] == expected_head, "uniform count head mismatch")
        check(counts[-1] == 3*counts[3]+2, "converse ready continuation mismatch")
        for rows in graphs[1:]:
            check(induced(rows,list(range(len(retained)))) == deleted,
                  "persistent old arc changed")
        check(all(strong(rows) for rows in graphs), "intermediate strongness failed")
        reports.append({"k":k,"retained_original_vertices":retained,"count_trace":counts,
                        "outmasks":[out,out,1 << (k+2),1 << len(reset)],
                        "old_count_restored":counts[3] == counts[0],
                        "rows_at_every_stage":[list(t) for t in graphs]})
    chosen = reports[1]
    mapping = [1,2,3,6,7,8,0,4,5]
    u = (254,508,504,240,480,448,384,256,9)
    reset = chosen["rows_at_every_stage"][3]
    check(all(((reset[v] >> w)&1) == ((u[mapping[v]] >> mapping[w])&1)
              for v in range(9) for w in range(9)), "persistent U isomorphism failed")
    check(chosen["count_trace"] == [141,33,69,141,215,425], "141 trace mismatch")
    report = {"status":"PASS independent literal audit; all-k formula also checked algebraically in audit note",
              "reviewed_scope":"RESET sections1-2; exhaustive broader searches and neutral component not independently rerun here",
              "literal_family_controls":reports,"current_order_to_U_labels":mapping,
              "seven_persistent_vertex_ledger":dict(zip([0,1,2,3,4,5,7],mapping[:7])),
              "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name("reset_audit.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":"PASS","traces":[r["count_trace"] for r in reports],
                      "persistent_vertex_map":report["seven_persistent_vertex_ledger"]},indent=2))


if __name__ == "__main__":
    main()
