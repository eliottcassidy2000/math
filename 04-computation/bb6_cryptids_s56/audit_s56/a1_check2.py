import sys, os, re, random, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(os.path.dirname(HERE), 'session_copy'))
from tmcheck import check, HALT
import a1_check as M   # re-runs nothing heavy? (module body runs) -> avoid: copy minimal pieces instead
