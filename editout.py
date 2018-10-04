import sys
flnm = sys.argv[1]
f = open(flnm)
g = open('ed%s'%flnm , 'w')
i = 0
import re
for line in f:
    l = re.sub(r'0+;',';', line)
    ln = re.sub(r'(\.;)', '; ', l)
    ln = re.sub(r'(\.\.\.\n)', ' ', ln)
    ln = ' '.join(ln.split())
    if '=' in ln: 
        g.write('\n')
    else:
        g.write(' ')
    g.write(ln)

f.close()
g.close()


