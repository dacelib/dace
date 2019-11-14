import sys
with open(sys.argv[1]) as f:
    for line in f:
       line = line.rstrip()
       if line[0:2] == 'CD':
          line = '  ' + line[2:].ljust(74,' ') + 'CD'
       print( line )
