import random
import time
import itertools


def read_file(file):
    f = open(file, "r")
    listseq = []
    lines = f.read().split()
    t = int(lines[0])
    n = int(lines[1])
    l = int(lines[2])
    for i in range(3, len(lines), 1):
        listseq.append(lines[i])
    f.close()
    return [t, n, l, listseq]


res = read_file("rawDNA.txt")
t = res[0]
n = res[1]
l = res[2]
seqlist = res[3]


def generate_random(t, n):
    l = int(input("please Enter size of the pattern you want"))
    randomlist = []
    for i in range(t):
        seq = ""
        for j in range(n):
            seq += random.choice("ACGT")
        randomlist.append(seq)
    return randomlist


def MinHammingDistance(pattern, sequence, verbose=False, k=None):
    pattern = pattern.upper()
    sequence = sequence.upper()
    min_distance = len(pattern)
    min_i = -1
    for i in range(len(sequence) - (len(pattern) if k == None else k) + 1):
        distance = sum(
            [
                1
                for j in range(len(pattern))
                if pattern[j] != sequence[i : i + len(pattern)][j]
            ]
        )

        min_distance = min(min_distance, distance)
        if distance == min_distance:
            min_i = i

    if verbose:
        print(
            f"Align Seq: {sequence}, Index: {min_i}, Dist: {min_distance}, Compared Pattern: {pattern}, Best Align: {sequence[min_i:min_i+len(pattern)]}"
        )
    return min_distance


def GenerateArray(k):
    bases = ["A", "C", "G", "T"]
    array = bases
    for n in range(k - 1):
        array = [i + j for i in array for j in bases]
    return array


def FindMedianString(k, dna):
    pattern = GenerateArray(k)
    distance_of_pattern_dna = {}
    min_string = len(dna) * len(pattern)
    for i in pattern:
        sum_distance = 0
        for j in range(len(dna)):
            sum_distance += MinHammingDistance(i, dna[j])
        distance_of_pattern_dna[i] = sum_distance
        if sum_distance < min_string:
            min_string = sum_distance
    for t in distance_of_pattern_dna.keys():
        if distance_of_pattern_dna[t] == min_string:
            return t


def branchAndBoundMedianStringReq(DNA, k, pattern, min_dist, best_pattern):
    ham_dists_sum = 0
    for seq in DNA:
        ham_dists_sum += MinHammingDistance(pattern, seq, k=k)

    if min_dist < ham_dists_sum:
        return min_dist, best_pattern

    if len(pattern) == k:
        return ham_dists_sum, pattern

    for c in "ACGT":
        node_min_dist, node_best_pattern = branchAndBoundMedianStringReq(
            DNA, k, pattern + c, min_dist, best_pattern
        )

        if node_min_dist < min_dist:
            min_dist = node_min_dist
            best_pattern = node_best_pattern

    return min_dist, best_pattern


def branchAndBoundMedianString(DNA, k):
    return branchAndBoundMedianStringReq(DNA, k, "", float("inf"), ["A"] * l)[1]


print("choose a IF you want to run algorithm  from file")
print("choose b IF you want to run algorithm  ON Randomly generated sequences")
choice = input("Please Choose a or b  ")
# choice = "a"

l = 0
seqlist = []

if choice == "a":
    file = input("please Enter path of the file :")
    # file = "rawDNA.txt"
    result = read_file(file)
    # t = result[0]
    # n = result[1]
    l = result[2]
    seqlist = result[3]
elif choice == "b":
    seqlist = generate_random(t, n)
    l = input("Please enter l: ")


print("\n=======================================================")
print("Calculating median string with branch and bound...")
start = time.time()
mStr = branchAndBoundMedianString(seqlist, l)
print("Branch & Bound Median String:", mStr)
end = time.time()
print("Branch & Bound Time:", end - start)
print()
for seq in seqlist:
    MinHammingDistance(mStr, seq, True)
print("=======================================================\n")

print("\n=======================================================")
print("Calculating median string with brute force...")
start = time.time()
mStr = FindMedianString(l, seqlist)
print("Brute Force Median String:", mStr)
end = time.time()
print("Brute Force Time:", end - start)
print()
for seq in seqlist:
    MinHammingDistance(mStr, seq, True)
print("=======================================================\n")

print("\nConsensus String: ", mStr)
