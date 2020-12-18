from typing import Dict, List


UP = (-1, 0)
LEFT = (0, -1)
TOPLEFT = (-1, -1)
ORIGIN = (0, 0)


class NeedlemanWunsch():
    def __init__(self, scores: Dict = None, keys: List = None, delta: Dict = None):
        # score is a dictionary that record scores for match, mismatch, gap penalty
        # keys are all the characters that will appear in the sequences
        self.scores = scores
        self.keys = keys
        self.delta = delta

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

    def align(self, v: str, w: str):
        if self.delta == None:
            delta = self.delta_score()
        else:
            delta = self.delta
        M = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]
        pointers = [[ORIGIN for j in range(len(w) + 1)] for i in range(len(v) + 1)]
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
        score = M[len(v)][len(w)]
        alignment_v, alignment_w = self.traceback_global(v, w, pointers)
        return score, (alignment_v, alignment_w)


if __name__ == '__main__':
    keys = ['A', 'C', 'T', 'G', '-']
    scores = {'match': 1,
              'mismatch': -1,
              'gap': -1}
    nw = NeedlemanWunsch(scores, keys)
    score, alignment = nw.align('TAGATA', 'GTAGGCTTAAGGTTA')
    print(score)
    print(alignment[0])
    print(alignment[1])
