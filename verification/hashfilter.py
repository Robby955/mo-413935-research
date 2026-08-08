"""Keep only the lines whose stable hash falls in one residue class.

Used to split a level of the tower into passes so that no single sort has to see
the whole raw stream.  The union over p = 0..P-1 is the whole stream.
"""

import hashlib
import sys

P = int(sys.argv[1])
p = int(sys.argv[2])
out = sys.stdout
for line in sys.stdin:
    s = line.strip()
    if not s:
        continue
    if int.from_bytes(hashlib.blake2b(s.encode(), digest_size=8).digest(), "big") % P == p:
        out.write(s + "\n")
