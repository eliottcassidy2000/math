"""Extract already-counted start defects; no new graph query or search."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
path = HERE/"receipt.json"
source = json.loads(path.read_text())
records = []
for record in source["records"]:
    target = record["target_sidecar"]
    defects = {name: value["H_X_Y"][1]+value["H_X_Y"][2]-target["H"]
               for name, value in target["ports"].items()}
    if not all(value >= 0 for value in defects.values()):
        raise RuntimeError("Start defects must be nonnegative")
    records.append({"family": record["family"], "t": record["t"], "H": target["H"],
                    "sigma_S2_minus_H": defects,
                    "sigma1_ports": [int(v) for v,d in defects.items() if d == 1],
                    "sigma2_ports": [int(v) for v,d in defects.items() if d == 2],
                    "sigma0_ports": [int(v) for v,d in defects.items() if d == 0]})
report = {"status": "FINITE-EXACT extraction from existing endpoint counts",
          "formula": "sigma_v = S_v(2)-H = X_v+Y_v-H",
          "scope": "Only ports99,0,1,2,3,4; no claim about other old/promoted/new ports.",
          "receipt_sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "records": records}
(HERE/"START_DEFECTS.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps(records,indent=2))
