#!/usr/bin/env python

"""
Truncates the first and last a,b tokens from each line,
where a,b are arguments from sys.argv.
(The truncate_char.py on the other hand, truncates individual characters)
Strings are split into "tokens" use white-space as a delimiter.
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
        tokens = line.split()
        if (len(line) > (truncate_a + truncate_b)):
            truncated_tokens = tokens[truncate_a:len(tokens)-truncate_b]
        else:
            truncated_tokens = []
        for i in range(0,len(truncated_tokens)):
            sys.stdout.write(truncated_tokens[i])
            if i+1 < len(truncated_tokens):
                sys.stdout.write(' ')
        sys.stdout.write('\n')


if __name__ == "__main__":
    main()
