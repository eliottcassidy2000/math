"""Run every producer twice; enforce identical receipts under -O."""
import hashlib
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
jobs = [("verify.py", "verification.json"),
        ("geometry_checks.py", "geometry_verification.json"),
        ("path_kernel.py", "path_kernel_verification.json")]
for script, receipt in jobs:
    subprocess.run([sys.executable, "-X", "utf8", "-B", str(root/script)],
                   check=True, capture_output=True)
    normal = (root/receipt).read_bytes()
    subprocess.run([sys.executable, "-X", "utf8", "-O", "-B", str(root/script)],
                   check=True, capture_output=True)
    optimized = (root/receipt).read_bytes()
    if normal != optimized:
        raise RuntimeError(f"Optimization changed {receipt}")
    print(f"PASS {script}: normal/-O identical; SHA256 {hashlib.sha256(normal).hexdigest()}")
