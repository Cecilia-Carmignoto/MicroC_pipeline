
import multiprocessing as mp
import itertools as it
import json
import os
import sys
import time
import unittest

class MatchKind:
    FULL_MATCH_FORWARD_MIDDLE = 1
    FULL_MATCH_REVERSE_MIDDLE = 2
    PARTIAL_MATCH_FORWARD_END = 3
    PARTIAL_MATCH_REVERSE_END = 4
    PARTIAL_MATCH_FORWARD_BEGINNING = 5
    PARTIAL_MATCH_REVERSE_BEGINNING = 6

    @staticmethod
    def str(x):
        return {
            MatchKind.FULL_MATCH_FORWARD_MIDDLE: "Full match forward middle",
            MatchKind.FULL_MATCH_REVERSE_MIDDLE: "Full match reverse middle",
            MatchKind.PARTIAL_MATCH_FORWARD_END: "Partial match forward end",
            MatchKind.PARTIAL_MATCH_REVERSE_END: "Partial match reverse end",
            MatchKind.PARTIAL_MATCH_FORWARD_BEGINNING: "Partial match forward beginning",
            MatchKind.PARTIAL_MATCH_REVERSE_BEGINNING: "Partial match reverse beginning"
        }[x]

class Match:
    def __init__(self, kind, length):
        self.kind = kind
        self.length = length

    def isfull(self):
        return self.kind in (
            MatchKind.FULL_MATCH_FORWARD_MIDDLE,
            MatchKind.FULL_MATCH_REVERSE_MIDDLE
        )

    def ispartial(self):
        return not self.isfull()
    
    def __repr__(self):
        return f"Match(kind={MatchKind.str(self.kind)}, length={self.length})"
    
    def __hash__(self):
        return hash((self.kind, self.length))
    
    def __eq__(self, other):
        return self.kind == other.kind and self.length == other.length

def match_partial_bridge_end(seq, bridge):
    for i in range(len(bridge), 1, -1):
        if seq.endswith(bridge[:i]):
            return i
    return 0

def _match_partial_bridge_beginning(seq, bridge):
    for i in range(len(bridge), 1, -1):
        if seq.startswith(bridge[-i:]):
            return i
    return 0
    
def _check_bridge(seq, bridge):
    possible_matches = set()
    rev_bridge = bridge[::-1]
    if bridge in seq:
        possible_matches.add(
            Match(MatchKind.FULL_MATCH_FORWARD_MIDDLE, len(bridge)))
    if rev_bridge in seq:
        possible_matches.add(
            Match(MatchKind.FULL_MATCH_REVERSE_MIDDLE, len(bridge)))
    match_fwd_end_count = match_partial_bridge_end(seq, bridge)
    if match_fwd_end_count > 0:
        possible_matches.add(
            Match(MatchKind.PARTIAL_MATCH_FORWARD_END, match_fwd_end_count))
    match_rev_end_count = match_partial_bridge_end(seq, rev_bridge)
    if match_rev_end_count > 0:
        possible_matches.add(
            Match(MatchKind.PARTIAL_MATCH_REVERSE_END, match_rev_end_count))
    match_fwd_beginning_count = _match_partial_bridge_beginning(seq, bridge)
    if match_fwd_beginning_count > 0:
        possible_matches.add(
            Match(MatchKind.PARTIAL_MATCH_FORWARD_BEGINNING, match_fwd_beginning_count))
    match_rev_beginning_count = _match_partial_bridge_beginning(seq, rev_bridge)
    if match_rev_beginning_count > 0:
        possible_matches.add(
            Match(MatchKind.PARTIAL_MATCH_REVERSE_BEGINNING, match_rev_beginning_count))

    return possible_matches

class PairKind:
    FULL = 1
    PARTIAL = 2
    NO_MATCH = 3

def bridge_count(seq1, seq2, bridge):
    max_pair = (0, 0)

    matches1 = _check_bridge(seq1, bridge)
    matches2 = _check_bridge(seq2, bridge)
    max_matches_1 = max([match.length for match in matches1], default=0)
    max_matches_2 = max([match.length for match in matches2], default=0)
    if any(match.isfull() for match in matches1):
        return (len(bridge), max_matches_2), PairKind.FULL
    if any(match.isfull() for match in matches2):
        return (max_matches_1, len(bridge)), PairKind.FULL
    if not matches1 or not matches2:
        return max_pair, PairKind.NO_MATCH
    
    # here we know that there are, and there are only, partial matches in both sequences
    # need to check for consistent matches, keeping the overall maximum
    max_match = 0
    for match1, match2 in it.product(matches1, matches2):
        match = min(match1.length, match2.length)
        if (match1.kind == MatchKind.PARTIAL_MATCH_FORWARD_END and match2.kind == MatchKind.PARTIAL_MATCH_FORWARD_BEGINNING) or \
            (match1.kind == MatchKind.PARTIAL_MATCH_REVERSE_END and match2.kind == MatchKind.PARTIAL_MATCH_REVERSE_BEGINNING) or \
            (match1.kind == MatchKind.PARTIAL_MATCH_FORWARD_BEGINNING and match2.kind == MatchKind.PARTIAL_MATCH_FORWARD_END) or \
            (match1.kind == MatchKind.PARTIAL_MATCH_REVERSE_BEGINNING and match2.kind == MatchKind.PARTIAL_MATCH_REVERSE_END):
            if match > max_match:
                max_match = match
                max_pair = (match1.length, match2.length)
    if max_match == 0:
        return (max_matches_1, max_matches_2), PairKind.NO_MATCH
    return max_pair, PairKind.PARTIAL

def process_chunk(chunk):
    fastq1_path, fastq2_path, bridge, start, end = chunk
    print(f"Chunk started: {start:_} - {end:_}")
    full, null = 0, 0
    full_matches, partial_matches = [], []
    with open(fastq1_path, 'r') as fastq1, open(fastq2_path, 'r') as fastq2:
        fastq1.seek(start)
        fastq2.seek(start)

        full, null = 0, 0
        full_matches, partial_matches, null_matches = [], [], []
        tot = 0
        start_time = time.time()
        lap_time = start_time
        while True:
            header1 = fastq1.readline()
            header2 = fastq2.readline()
            id1 = header1.split(':')[0][1:]
            id2 = header2.split(':')[0][1:]

            # ids should match between the two files
            if id1 != id2:
                print(f"Error: IDs do not match: {id1} != {id2}")
                exit(1)
            
            if not header1 or not header2:
                print(f"Chunk {start:_} - {end:_} EOF")
                exit(1)

            seq1 = fastq1.readline()[:-1] # remove newline
            seq2 = fastq2.readline()[:-1] # remove newline
            tuple_lengths, pair_kind = bridge_count(seq1, seq2, bridge)

            if pair_kind == PairKind.FULL:
                full += 1
                full_matches.append({
                    # 'id': id1,
                    # 'seq1': seq1,
                    # 'seq2': seq2,
                    'lengths': tuple_lengths
                })
            elif pair_kind == PairKind.PARTIAL:
                partial_matches.append({
                    # 'id': id1,
                    # 'seq1': seq1,
                    # 'seq2': seq2,
                    'lengths': tuple_lengths
                })
            else:
                null += 1
                null_matches.append({
                    # 'id': id1,
                    # 'seq1': seq1,
                    # 'seq2': seq2,
                    'lengths': tuple_lengths
                })
            
            for _ in range(2):
                fastq1.readline()
                fastq2.readline()
            
            if fastq1.tell() >= end:
                break

            tot += 1
            if tot % 1_000_000 == 0:
                lap_time = time.time() - lap_time
                elapsed_time = time.time() - start_time
                print(f"Already analyzed {tot:_} lines on both files in {elapsed_time:.2f} seconds (lap time: {lap_time:.2f} seconds)")
                lap_time = time.time()
                print(f"Matches (Full: {full:_}, No matches: {null:_})")
    
    print(f"Chunk finished: {start:_} - {end:_} in {time.time() - start_time:.2f} seconds")
    return full, null, full_matches, partial_matches, null_matches

SEQUENCE_LENGTH = sum([51, 101, 2, 101]) # 255 bytes, divided each line (4 lines in a sequence)
def get_file_chunks(fastq1_path, fastq2_path, num_chunks):
    size = os.path.getsize(fastq1_path)
    if size != os.path.getsize(fastq2_path):
        print("Error: Files have different sizes")
        exit(1)
    number_of_sequences = size // SEQUENCE_LENGTH
    chunk_size = number_of_sequences // num_chunks
    print(f"Total size: {size:_} bytes, Chunk size: {chunk_size:_} sequences ({chunk_size * SEQUENCE_LENGTH:_} bytes)")
    chunks = []
    for i in range(num_chunks):
        start = i * chunk_size * SEQUENCE_LENGTH
        end = (i + 1) * chunk_size * SEQUENCE_LENGTH
        chunks.append((start, end))
    return chunks

# python3 bridge_count.py mC_wt_R1_48h_L1_1.fq mC_wt_R1_48h_L1_2.fq AGGTTCGTCCATCGATCGATGGACGAACCT
def parallel_main():
    fastq1_path, fastq2_path, bridge = sys.argv[1:4]
    np = mp.cpu_count() * 2
    print(f"Using {np} processes")
    chunks = get_file_chunks(fastq1_path, fastq2_path, np)
    with mp.Pool(np) as pool:
        results = pool.map(process_chunk, [(fastq1_path, fastq2_path, bridge, start, end) for start, end in chunks])
    print("All processes finished")

    full, null = 0, 0
    full_matches, partial_matches, null_matches = [], [], []
    for full_chunk, null_chunk, full_matches_chunk, partial_matches_chunk, null_matches_chunk in results:
        full += full_chunk
        null += null_chunk
        full_matches.extend(full_matches_chunk)
        partial_matches.extend(partial_matches_chunk)
        null_matches.extend(null_matches_chunk)

    print(f"Full matches: {full:_}")
    print(f"Partial matches: {len(partial_matches):_}")
    print(f"No matches: {null:_}")
    print()

    with open("output.json", 'w') as output:
        json.dump({
            'full': full,
            'null': null,
            'full_matches': full_matches,
            'partial_matches': partial_matches,
            'null_matches': null_matches
        }, output)
    print("Output written to output.json")

# python3 bridge_count.py mC_wt_R1_48h_L1_1.fq mC_wt_R1_48h_L1_2.fq AGGTTCGTCCATCGATCGATGGACGAACCT
def main():
    fastq1_path, fastq2_path, bridge = sys.argv[1:4]

    tot = 0
    full, null = 0, 0
    full_matches, partial_matches = [], []
    start_time = time.time()
    lap_time = start_time
    with open(fastq1_path, 'r') as fastq1, open(fastq2_path, 'r') as fastq2:
        while True:
            header1 = fastq1.readline()
            header2 = fastq2.readline()
            id1 = header1.split(':')[0][1:]
            id2 = header2.split(':')[0][1:]
            if id1 != id2:
                print(f"Error: IDs do not match: {id1} != {id2}")
                break

            if not header1 or not header2:
                break

            seq1 = fastq1.readline()[:-1] # remove newline
            seq2 = fastq2.readline()[:-1] # remove newline
            tuple_lengths, pair_kind = bridge_count(seq1, seq2, bridge)

            if pair_kind == PairKind.FULL:
                full += 1
                full_matches.append({
                    'id': id1,
                    'seq1': seq1,
                    'seq2': seq2,
                    'lengths': tuple_lengths
                })
            elif pair_kind == PairKind.PARTIAL:
                partial_matches.append({
                    'id': id1,
                    'seq1': seq1,
                    'seq2': seq2,
                    'lengths': tuple_lengths
                })
            else:
                null += 1

            for _ in range(2):
                fastq1.readline()
                fastq2.readline()

            tot += 1
            if tot % 1_000_000 == 0:
                lap_time = time.time() - lap_time
                elapsed_time = time.time() - start_time
                print(f"Already analyzed {tot:_} lines on both files in {elapsed_time:.2f} seconds (lap time: {lap_time:.2f} seconds)")
                lap_time = time.time()
                print(f"Matches (Full: {full:_}, No matches: {null:_})")
                print()

    print(f"Total of {tot:_} lines analyzed in {time.time() - start_time:.2f} seconds")
    print(f"Full matches: {full:_}")
    print(f"No matches: {null:_}")
    print()

    with open("output.json", 'w') as output:
        json.dump({
            'full': full,
            'null': null,
            'full_matches': full_matches,
            'partial_matches': partial_matches
        }, output)
    print("Output written to output.json")

# python3 -m unittest bridge_count.py  
class TestCheckBridge(unittest.TestCase):
    def test_no_match(self):
        seq = "CCACCTTATTAACATTCC"
        bridge = "ACT"
        matches = _check_bridge(seq, bridge)
        self.assertEqual(matches, set())

    def test_no_full_match(self):
        seq = "ACCTTATTAACATT"
        bridge = "ACT"
        matches = _check_bridge(seq, bridge)
        full_matches = any(match.isfull() for match in matches)
        self.assertFalse(full_matches)

    def test_full_match(self):
        seq = "ACCTTATTAACTT"
        bridge = "ACT"
        matches = _check_bridge(seq, bridge)
        full_matches = any(match.isfull() for match in matches)
        self.assertTrue(full_matches)

    def test_full_reverse_match(self):
        seq = "ACCTTATTCAATT"
        bridge = "ACT"
        matches = _check_bridge(seq, bridge)
        full_matches = any(match.isfull() for match in matches)
        self.assertTrue(full_matches)

    def test_partial_match_reverse_end(self):
        seq = "CCATCTTATTAACATTCA"
        bridge = "ACT"
        matches = _check_bridge(seq, bridge)
        self.assertIn(Match(MatchKind.PARTIAL_MATCH_REVERSE_END, 3), matches)

    def test_partial_match_forward_end(self):
        seq = "CCACCTTATTAACATTCAC"
        bridge = "ACT"
        matches = _check_bridge(seq, bridge)
        self.assertIn(Match(MatchKind.PARTIAL_MATCH_FORWARD_END, 2), matches)
        self.assertIn(Match(MatchKind.FULL_MATCH_REVERSE_MIDDLE, 3), matches)

if __name__ == "__main__":
    parallel_main()