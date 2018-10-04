import sys
# flnm = sys.argv[1]
f = open('cout')
g = open('edcout', 'w')
i = 0
import re
for line in f:
    l = re.sub(r'0+;',';', line)
    ln = re.sub(r'(\.;)', ';', l)
    g.write(ln)

f.close()
g.close()


