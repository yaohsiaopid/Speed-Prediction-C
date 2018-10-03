import sys
flnm = sys.argv[1]
f = open(flnm)
i = 0
import re
for line in f:
    if i >= 6:
        nm = re.search(r'(\w*)', line)[0]
        if len(nm) != 0:
            print()
            print(nm)
        nums = re.findall(r"[-+]?\d*\.\d+|\d+", line)
        for n in nums:
            print('%s; ' % n, end='')
        
    i += 1

