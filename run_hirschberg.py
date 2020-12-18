import os
import time
from hirschberg import Hirschberg

def read_delta(path):
    delta = {}
    with open(path, 'r') as f:
        lines = f.readlines()
        keys = lines[0].split()
        keys[-1] = '-'
        for i, line in enumerate(lines[1:]):
            delta[keys[i]] = {k: int(v) for (k, v) in zip(keys, line.split()[1:])}
    return delta, keys

def main(args):
    path = None
    if args.path is not None:
        path = args.path
        if not os.path.isdir(path):
            os.mkdir(path)

    file_name = args.file_name

    alignment1, alignment2 = None, None
    if args.alignment1 is not None and args.alignment2 is not None:
        alignment1 = args.alignment1
        alignment2 = args.alignment2
        if os.path.isfile(alignment1):
            alignment = ''
            with open(alignment1, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    alignment += line.strip()
            alignment1 = alignment
        if os.path.isfile(alignment2):
            alignment = ''
            with open(alignment2, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    alignment += line.strip()
            alignment2 = alignment
    else:
        print("Missing one or more sequences to align with!")
        exit(0)

    delta = None
    score = None
    keys = None
    if args.delta is not None:
        delta, keys = read_delta(args.delta)
    else:
        if args.match is not None and args.mismatch is not None and args.gap is not None:
            score = {'match': args.match,
                     'mismatch': args.mismatch,
                     'gap': args.gap}
        if args.keys is not None:
            keys = args.keys.split(',')
        else:
            print("Symbols will be used by the sequences are missing!")
            exit(0)

    hs = Hirschberg(score, keys, delta)
    start_time = time.time()
    score, alignments = hs.align(alignment1, alignment2)
    end_time = time.time()
    elapsed_time = end_time - start_time
    if path is None:
        print("Best Alignment Score:", score)
        print("Sequence 1: ", alignments[0])
        print("Sequence 2: ", alignments[1])
        print("Alignment is done in %.4f seconds!" % elapsed_time)
    else:
        with open(os.path.join(path, file_name), 'w+') as f:
            f.write("Best Alignment Score: %s \n" % str(score))
            f.write("Sequence 1: %s \n" % alignments[0])
            f.write("Sequence 2: %s \n" % alignments[1])
        print("Alignment is done in %.4f seconds!" % elapsed_time)
        print("Result saved at %s" % (path + file_name))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, default=None, help='writing directory for the output alignment; '
                                                               'if want output in terminal, leave this blank')
    parser.add_argument('--file_name', type=str, default='hirschberg_output.txt', help='file to write final output')
    parser.add_argument('--keys', type=str, default=None, help='symbols in the sequences')
    parser.add_argument('--alignment1', type=str, default=None, help='first sequence for alignment;'
                                                                     'can be either an alignment or '
                                                                     'a txt file containing the alignment')
    parser.add_argument('--alignment2', type=str, default=None, help='second sequence for alignment;'
                                                                     'can be either an alignment or '
                                                                     'a txt file containing the alignment')
    parser.add_argument('--delta', type=str, default=None, help='file path to delta scoring function; '
                                                                'leaving this blank will cause the program to generate a naive one based on provided match, mismatch, gap scores')
    parser.add_argument('--match', type=int, default=1, help='score for matched symbol')
    parser.add_argument('--mismatch', type=int, default=-1, help='score for mismatched symbor')
    parser.add_argument('--gap', type=int, default=-1, help='gap penalty for scoring matrix')
    args = parser.parse_args()

    main(args)
