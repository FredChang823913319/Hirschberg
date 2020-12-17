from typing import Dict, List


UP = (-1, 0)
LEFT = (0, -1)
TOPLEFT = (-1, -1)
ORIGIN = (0, 0)

class Hirschberg():
    def __init__(self, scores: Dict, keys: List):
        # score is a dictionary that record scores for match, mismatch, gap penalty
        # keys are all the characters that will appear in the sequences
        self.scores = scores
        self.keys = keys

    def delta_score(self):
        delta = {}
        l = len(self.keys)
        for i in range(l):
            char1 = self.keys[i]
            if char1 == '-':
                delta['-'] = {k: v for (k, v) in zip(self.keys, [self.scores['gap'] for _ in range(l)])}
                continue
            delta[char1] = {}
            for j in range(l):
                char2 = self.keys[j]
                if char2 == '-':
                    delta[char1][char2] = self.scores['gap']
                elif i == j:
                    delta[char1][char2] = self.scores['match']
                else:
                    delta[char1][char2] = self.scores['mismatch']
        return delta

    @staticmethod
    def prefix(v, w, delta):
        n, m = len(v), len(w)
        M = [[0 for _ in range(m+1)] for _ in range(n+1)]
        for i in range(n+1):
            for j in range(m+1):
                if i == 0 and j == 0:
                    M[i][j] = 0
                elif i == 0:
                    M[i][j] = M[i][j - 1] + delta['-'][w[j - 1]]
                elif j == 0:
                    M[i][j] = M[i - 1][j] + delta[v[i - 1]]['-']
                else:
                    M[i][j] = max(M[i][j - 1] + delta['-'][w[j - 1]],
                                  M[i - 1][j] + delta[v[i - 1]]['-'],
                                  M[i-1][j-1] + delta[v[i - 1]][w[j-1]])
            # Since we only need two row to calculate, we clear the prev row from memory
            if i != 0:
                M[i-1] = []
        return M[n]

    @staticmethod
    def suffix(v, w, delta):
        n, m = len(v), len(w)
        M = [[0 for _ in range(m+1)] for _ in range(n+1)]
        # we need to reverse every edge in calculating the suffix
        # so we calculate score from n-i and m-j
        for i in range(n+1):
            for j in range(m+1):
                if i == 0 and j == 0:
                    M[i][j] = 0
                elif i == 0:
                    M[i][j] = M[i][j - 1] + delta['-'][w[m - j]]
                elif j == 0:
                    M[i][j] = M[i - 1][j] + delta[v[n - i]]['-']
                else:
                    M[i][j] = max(M[i][j - 1] + delta['-'][w[m - j]],
                                  M[i - 1][j] + delta[v[n - i]]['-'],
                                  M[i-1][j-1] + delta[v[n - i]][w[m - j]])
            # Since we only need two row to calculate, we clear the prev row from memory
            if i != 0:
                M[i-1] = []
        return M[n]

    @staticmethod
    def traceback_global(v, w, pointers):
        i, j = len(v), len(w)
        new_v = []
        new_w = []
        while True:
            di, dj = pointers[i][j]
            if (di, dj) == LEFT:
                new_v.append('-')
                new_w.append(w[j - 1])
            elif (di, dj) == UP:
                new_v.append(v[i - 1])
                new_w.append('-')
            elif (di, dj) == TOPLEFT:
                new_v.append(v[i - 1])
                new_w.append(w[j - 1])
            i, j = i + di, j + dj
            if i <= 0 and j <= 0:
                break
        return ''.join(new_v[::-1]), ''.join(new_w[::-1])

    def base_case(self, v, w, delta):
        M = [[0 for _ in range(len(w) + 1)] for _ in range(len(v) + 1)]
        pointers = [[ORIGIN for _ in range(len(w) + 1)] for _ in range(len(v) + 1)]
        for i in range(len(v) + 1):
            for j in range(len(w) + 1):
                if i == 0 and j == 0:
                    M[i][j] = 0
                elif i == 0:
                    M[i][j] = M[i][j - 1] + delta['-'][w[j - 1]]
                    pointers[i][j] = LEFT
                elif j == 0:
                    M[i][j] = M[i - 1][j] + delta[v[i - 1]]['-']
                    pointers[i][j] = UP
                else:
                    best_sub = max([(LEFT, M[i][j - 1] + delta['-'][w[j - 1]]),
                                    (UP, M[i - 1][j] + delta[v[i - 1]]['-']),
                                    (TOPLEFT, M[i - 1][j - 1] + delta[v[i - 1]][w[j - 1]])], key=lambda x: x[1])
                    pointers[i][j] = best_sub[0]
                    M[i][j] = best_sub[1]
            if i != 0:
                M[i-1] = []
        score = M[len(v)][len(w)]
        alignment_v, alignment_w = self.traceback_global(v, w, pointers)
        return score, alignment_v, alignment_w

    def hirschberg(self, v, w, delta):
        n, m = len(v), len(w)
        # in cases where any one of the sequence has only one character, we can directly use needleman-wunsch
        # as it gurantees linear space (since the table size would be 1 * m or 1 * n)
        if n <= 1 or m <= 1:
            return self.base_case(v, w, delta)
        else:
            prefix = self.prefix(v[:n//2], w, delta)
            suffix = self.suffix(v[n//2:], w, delta)
            index = [prefix[i] + suffix[m - i] for i in range(m + 1)]
            index = index.index(max(index))
            prefix = self.hirschberg(v[:n//2], w[:index], delta)
            suffix = self.hirschberg(v[n//2:], w[index:], delta)
            score = prefix[0] + suffix[0]
            alignment_v = prefix[1] + suffix[1]
            alignment_w = prefix[2] + suffix[2]
            return score, alignment_v, alignment_w

    def align(self, v, w):
        delta = self.delta_score()
        score, alignment_v, alignment_w = self.hirschberg(v, w, delta)
        return score, (alignment_v, alignment_w)


if __name__ == '__main__':
    keys = ['A', 'C', 'T', 'G', '-']
    scores = {'match': 1,
              'mismatch': -1,
              'gap': -1}
    nw = Hirschberg(scores, keys)
    score, alignment = nw.align('TAGATA', 'GTAGGCTTAAGGTTA')
    print(score)
    print(alignment[0])
    print(alignment[1])