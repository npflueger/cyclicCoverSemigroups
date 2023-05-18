from genSemigroups import *
import sys

maxg = int(sys.argv[1])

print(f"This script will generate all semigroups arising from cyclic covers with N = 3, a-vector selected from the feasible set, and with genus <= {maxg}.")

for g in range(maxg+1):
    print(f"\n**** Now lisiting example in genus {g}\n")
    count = listByGenus(g)
    if count == 0: print("(None found)")
    else: print(f"{count} semigroups found in genus {g}")
