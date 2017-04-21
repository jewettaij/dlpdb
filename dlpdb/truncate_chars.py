#!/usr/bin/env python

"""
This script throws away the leading and trailing characters on each line.
It truncates the first and last a,b characters from each line, where a,b 
are arguments from sys.argv.  (In retrospect, I should have used awk.)
"""

import sys

def main():
    if (len(sys.argv) <= 1):
        sys.stderr.write('Error: expected the number of truncations as an argument\n')
    truncate_a = int(sys.argv[1])
    if (len(sys.argv) >= 3):
        truncate_b = int(sys.argv[2])
    else:
        truncate_b = truncate_a

    for line in sys.stdin:
        line = line.strip()
        if (len(line) > (truncate_a + truncate_b)):
            truncated_line = line[truncate_a:len(line)-truncate_b]
        else:
            truncated_line = ''
        sys.stdout.write(truncated_line+'\n')


if __name__ == "__main__":
    main()
