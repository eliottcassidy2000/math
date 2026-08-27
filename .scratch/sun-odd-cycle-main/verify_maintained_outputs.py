from pathlib import Path
import subprocess
import sys


PAIRS = [
    ("04-computation/sun_2468_reflection_orbit_tower_thm4246.py", "05-knowledge/results/sun_2468_reflection_orbit_tower_thm4246.out"),
    ("04-computation/sun_2468_reflection_orbit_tower_independent_audit_thm4246.py", "05-knowledge/results/sun_2468_reflection_orbit_tower_independent_audit_thm4246.out"),
    ("04-computation/jc23_w0_involution_projection_histogram_thm4247.py", "05-knowledge/results/jc23_w0_involution_projection_histogram_thm4247.out"),
    ("04-computation/jc23_w0_hidden_degree12_attachment_audit_thm4247.py", "05-knowledge/results/jc23_w0_hidden_degree12_attachment_audit_thm4247.out"),
    ("04-computation/jc23_w0_hidden_degree12_orbit0_symbolic_audit_thm4247.py", "05-knowledge/results/jc23_w0_hidden_degree12_orbit0_symbolic_audit_thm4247.out"),
    ("04-computation/lrc14_thm4231_4238_4242_residual_odd_cycle_core_20260826.py", "05-knowledge/results/lrc14_thm4231_4238_4242_residual_odd_cycle_core_20260826.out"),
    ("04-computation/lrc14_thm4231_4238_4242_residual_odd_cycle_core_independent_20260826.py", "05-knowledge/results/lrc14_thm4231_4238_4242_residual_odd_cycle_core_independent_20260826.out"),
]


for script, output in PAIRS:
    actual = subprocess.check_output([sys.executable, "-B", script]).replace(b"\r\n", b"\n")
    expected = Path(output).read_bytes().replace(b"\r\n", b"\n")
    assert actual == expected, (script, len(actual), len(expected))
    optimized = subprocess.check_output([sys.executable, "-B", "-O", script]).replace(b"\r\n", b"\n")
    assert optimized == expected, (script, "optimized")
    print(f"{script}: normal/optimized/frozen PASS")
