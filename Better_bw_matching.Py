def preprocess_bwt(bwt):

first_col = ''.join(sorted(bwt))

    count = {}
    for char in set(bwt):
        count[char] = [0] * (len(bwt) + 1)

    # Fill the count table
    for i in range(1, len(bwt) + 1):
        char = bwt[i-1]
        for c in count:
            count[c][i] = count[c][i-1] + (1 if c == char else 0)

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

bwt = "GGCGCCGC$TAGTCACACACGCCGTA"
patterns = ["ACC", "CCG", "CAG"]

result = better_bw_matching(bwt, patterns)
print(' '.join(map(str, result)))
