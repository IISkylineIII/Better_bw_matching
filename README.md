# Better_bw_matching

# Description
Python implementation of the Better Burrows-Wheeler Matching algorithm.
It efficiently finds the number of occurrences of given patterns in a string, using the Burrows-Wheeler Transform (BWT), a fundamental technique in data compression and bioinformatics.

```

def preprocess_bwt(bwt):
    # First column is the sorted BWT
    first_col = ''.join(sorted(bwt))

    # Count table to store counts of characters up to each position in BWT
    count = {}
    for char in set(bwt):
        count[char] = [0] * (len(bwt) + 1)

    # Fill the count table
    for i in range(1, len(bwt) + 1):
        char = bwt[i-1]
        for c in count:
            count[c][i] = count[c][i-1] + (1 if c == char else 0)

    # Occurrence map to store the first occurrence of each character in the first column
    occ = {}
    for i, char in enumerate(first_col):
        if char not in occ:
            occ[char] = i

    return first_col, count, occ

def count_occurrences(pattern, bwt, first_col, count, occ):
    top = 0
    bottom = len(bwt) - 1

    while top <= bottom and pattern:
        symbol = pattern[-1]
        pattern = pattern[:-1]

        if symbol in bwt[top:bottom+1]:
            top = occ[symbol] + count[symbol][top]
            bottom = occ[symbol] + count[symbol][bottom + 1] - 1
        else:
            return 0
    return bottom - top + 1

def better_bw_matching(bwt, patterns):
    first_col, count, occ = preprocess_bwt(bwt)
    return [count_occurrences(pattern, bwt, first_col, count, occ) for pattern in patterns]

# Input and output handling
bwt = "GGCGCCGC$TAGTCACACACGCCGTA"
patterns = ["ACC", "CCG", "CAG"]

# Call the function and print results
result = better_bw_matching(bwt, patterns)
print(' '.join(map(str, result)))

```

# Example Output

1 2 1

# Applications
* Efficient pattern matching using the Burrows-Wheeler Transform.
* Useful in bioinformatics for genome string matching.
* A fundamental component of FM-indexes for compressed full-text indexing.

# License
This project is licensed under the MIT License.

