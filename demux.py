import csv
from Bio import pairwise2
from Bio import SeqIO
from Bio.Seq import Seq
from Levenshtein import distance
import gzip
import sys

def find_best_match(sequence, reference_sequences):
    best_id = None
    best_distance = float('inf')
    for ref_id, ref_seq in reference_sequences.items():
        dist = distance(sequence, ref_seq)
        if dist <= 2 and dist < best_distance:
            best_id = ref_id
            best_distaece = dist
    return best_id

def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())

#fastq_file_r2 = "R1.fastq.gz"
#fastq_file_r1 = "R2.fastq.gz"
fastq_file_r2 = sys.argv[1]
fastq_file_r1 = sys.argv[2]
target_sequence_1 = "GTTGTAGCTCCCTTTCTCATTTCGG"
target_sequence_2 = "TTTAAGGGGCATCGTTTATTTTTTG"
ref_file = "demuxRef.csv"

# User-defined extraction coordinates
sp_start = 24
sp_end = 0  # Exclusive, adjust to desired length
bc_start = 0
bc_end_8 = 13  # Exclusive, adjust to desired length

output_list = []
no_match_list = []

# Load reference sequences from demuxRef.csv
reference_sequences = {}
with open(ref_file, "r") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        reference_sequences[row["id"]] = row["sp"]

reference_sequences_bc = {}
with open(ref_file, "r") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        reference_sequences_bc[row["id"]] = row["bc"]

# Initialize counts dictionary to keep track of the number of matched reads for each ID
counts = {id: 0 for id in reference_sequences}

with gzip.open(fastq_file_r1, "rt") as hr1, gzip.open(fastq_file_r2, "rt") as hr2:
    for record_r1, record_r2 in zip(SeqIO.parse(hr1, "fastq"), SeqIO.parse(hr2, "fastq")):
        merged_sequence = str(record_r1.seq) + reverse_complement(str(record_r2.seq))
        merged_sequence = reverse_complement(merged_sequence)
        quality_scores = record_r1.letter_annotations["phred_quality"] + record_r2.reverse_complement().letter_annotations["phred_quality"]
        average_quality = sum(quality_scores) / len(quality_scores)

        if average_quality >= 26:
            # Align target_sequence_1 to the merged sequence
            alignments_1 = pairwise2.align.localms(merged_sequence, target_sequence_1, 2, -1, -2, -2, one_alignment_only=True)
            # Align target_sequence_2 to the merged sequence
            alignments_2 = pairwise2.align.localms(merged_sequence, target_sequence_2, 2, -1, -2, -2, one_alignment_only=True)

            if alignments_1 and alignments_2:
                best_alignment_1 = alignments_1[0]
                best_alignment_2 = alignments_2[0]
                start_1 = best_alignment_1[3]
                end_2 = best_alignment_2[4]

                # Extract sp and bc sequences from the aligned regions
                sp = merged_sequence[start_1 - sp_start : start_1 - sp_end]
                bc = merged_sequence[end_2 + bc_start : end_2 + bc_end_8]

                # Find the best matches for sp and bc using reference sequences
                best_match_sp_id = find_best_match(sp, reference_sequences)
                best_match_bc_id = find_best_match(bc, reference_sequences_bc)

                if best_match_sp_id and best_match_bc_id:
                    output_list.append((best_match_sp_id, best_match_bc_id))
                    counts[best_match_sp_id] += 1

                    if best_match_sp_id == best_match_bc_id:
                        with open(f"{best_match_sp_id}_matched.fastq", "a") as f:
                            SeqIO.write(record_r1, f, "fastq")
                            
                else:
                    no_match_list.append(merged_sequence)
            else:
                no_match_list.append(merged_sequence)

# Save output_list to a CSV file
with open("output_list.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["best_match_sp_id", "best_match_bc_id"])
    writer.writerows(output_list)

# Save no_match_list to a CSV file
with open("no_match_list.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Sequence"])
    writer.writerows([[sequence] for sequence in no_match_list])

print("Output lists saved to CSV files: output_list.csv and no_match_list.csv")
