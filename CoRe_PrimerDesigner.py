# -*- coding: utf-8 -*-

"""
@author: Javier Abellán García ; contact: abellangarciajavier@gmail.com

This script filters the best conserved regions of an alignment for primer-design. It can be used in PCR, RT-PCR, qPCR, RT-qPCR and NGS library development, though NGS would require modifications on the size_length input to use for overlapping regions.

    Instructions:

Step 1. Perform a Multiple Sequence Alignment (MSA) on your sequences of interest and get a FASTA file ONLY with the consensus sequence for your alignment.

Step 2. Get the identity values for each nucleotide (nt) of the consensus sequence. For example
    2.1. Download "Geneious Prime: free version".
    2.2. Import the alignment file or use Geneious Prime to do the MSA.
    2.3. Select the alignment file -> 'export' -> 'export other' -> 'Graphs'
    2.4. On 'Graphs to export' only select the 'identity' checkbox, and on 'export options' format has to be csv and the checkbox 'exclude columns where the consensus sequence is a gap'. Click 'OK' and save the file to your preferred directory.

The identity_value of a nt is the percentage of the aligned sequences that have the same nt as the consensus sequence in a given position (expressed as the probability of a random sequence from the alignment having the same nt as the consensus sequence). The Identity.csv file contains the identity values for every nt of the consensus sequence in order.

Step 3. Open the Identity.csv file in Excel or any csv editor and remove the 'position' column and the 'title' row. The file should only have one column with no title, and each cell of the column should only contain a number between 0.00 and 1.00 (with '.' separated decimals).

Step 4. Modify this script as needed. Without changes, the script will use a consensus sequence and it's identity file to output a .csv file with all the possible primers (forward and reverse) for the consensus sequence, each with an individual identity score that indicates the conservation of their complementary region on the consensus sequence. Every segment of the code is carefully commented to make understanding what can be adapted to specific cases easier, but feel free to contact the author if you need help understanding or using it.

At the time of publishing this

Step 5. Running the code.
    5.1. Open the Terminal of your computer. If you're not computer-savy, and you're feeling overwhelmed, you can try looking up tutorials on the basics of using scripts for your operating system (macOS, Windows...). Make sure you have python installed on your computer, or you can use a python virtual environment (make sure to understand what this means before using the Terminal on your computer, research is key).
    5.2.Make sure to have the 'BioPython' library installed, otherwise the script won't be able to import and read the fasta files. You can install BioPython using pip by typing "pip install biopython" in a new line of your Terminal.
    5.3. You are ready to run the code. Either from a normal terminal window or opening one in the directory you'll be using (the latter recommended for simplicity), type "python3 CoRe_PrimerDesigner.py". The script should run smoothly if all the steps above have been followed.
    5.4. Tune the parameters in "SECTION 0: User-configurable parameters" to your desired goals.
    5.5. Answer the prompts. The terminal window will ask you to input some information, namely confirming the sequence type.
    5.6. If everything goes well, a message will appear saying that the file was created successfully.

The results will be two newly created files, both with the same name as the consensus sequence fasta file but one followed by "_primers.csv" and another followed by "_taqman.csv".


! Note: this script DOES NOT evaluate the experimental properties of each oligo, just identifies possible targets and uses theoretical properties to score the candidates. Further optimization is needed to make sure which of the output primers (if any) are good options for each project. Neither the author, nor any of the contributors, can guarantee the results of this script, and understanding every step of the process it is performing is strongly recommended.

Acknowledgements:
    -Bandara, Chamara (2021) 10.13140/RG.2.2.27207.21924/3. For the idea and the base for this Script. In a way, this script is a refactoring of the aforementioned, but it includes new capabilities such as the creation of reverse primers, a clearer prompt structure for inputs and the adaptation for sequences longer than 1000bp, which required an in-depth rewriting.
    -Carlota Monedero Herranz and Manuel Gámez Jurado. For their invaluable help in troubleshooting.

Don't forget to cite the code if you use it (either as is or modified) in you applications.
"""


import os
from Bio import SeqIO


# SECTION 0: User-Configurable Parameters.

"""FILES"""

# Write the filenames for consensus and identity files and for the working directory between quotation marks:
directory = ""
consensus_file = "Consensus.fa"  # Type the consensus filename (including .fa extension). "example_con.fa"
identity_file = "Identity.csv"  # Type the identity filename (including .csv extension). "example_id.csv"



"""PRIMER PROPERTIES"""

# Range of primer sizes (bp):
min_primer_length = 18
max_primer_length = 22


# Range of Tm values (ºC):
min_tm = 50.0
max_tm = 55.0


# Range of GC% values (%):
min_gc = 40.0
max_gc = 60.0


# Window size for GC-distribution analysis (bp):
window_size = 5


# GC-Clamp parameters
use_gc_clamp = True  # Set to False to disable GC-clamp prioritization.
clamp_length = 2  # Length of the GC-clamp (number of nt 3'->5').


# Weight of every parameter for primer quality assessment. The sum of all parameters should be 100 (it's a percentage).
gc_weight = 40  # Weight (%) assigned to %GC deviation from range in primer scoring.
self_comp_weight = 20  # Weight (%) assigned to self complementarity in primer scoring.
hairpin_weight = 20  # Weight (%) assigned to hairpins in primer scoring.
dist_weight = 5  # Weight (%) assigned to GC distribution in primer scoring.
gc_clamp_weight = 15  # Weight (%) assigned to GC-clamp in primer scoring. If GC-clamp prioritization is disabled this value will not be considered for the percentage.


# Parameters for Primer Ranking
cons_weight = 90
quality_weight = 10



"""PRIMER-PAIRS PROPERTIES"""

# Parameters for Product-size Range
range_choice = "Y"  # If you want to specify a range for PCR product size type "Y"
if range_choice == "Y":
    min_size = 150
    max_size = 250
else:
    min_size = max_size = None


# Parameters for Pair Formation
pair_tm_deviation = 3  # Maximum Tm difference allowed among primers from the same pair (ºC)
max_pair_complementarity = 3.0  # Maximum acceptable complementarity score between primers from the same pair


# Parameters for pair quality assessment
pair_comp_weight = 50
tm_diff_weight = 30
size_weight = 20


# Parameters for Pair Ranking
identity_weight = 90
pair_quality_weight = 10



"""TAQMAN PROBE PROPERTIES"""

num_primer_pairs = 10  # Number of primer pairs to analyze (top -> down)
num_probes_per_pair = 5  # Number of probes to propose per primer pair
min_distance_between_probes = 10  # Minimum distance between probes in nt. A low value will give very similar probes

min_probe_size = 15  # Minimum probe size
max_probe_size = 30  # Maximum probe size

ideal_tm_increase = 10  # Difference of Tm between the average of the primer-pairs and the TaqMan probe



# SECTION 1: File input.

def get_directory(directory):
    """Use current directory."""
    direction = directory
    return direction if direction else os.getcwd()

def get_file_path(directory, file):
    """Use the file name and return its path if exists."""
    filename = file
    file_path = os.path.join(directory, filename)
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Error: File '{file_path}' not found.")
    return file_path

def read_fasta(file_path):
    """Read a FASTA file and returns the sequence."""
    with open(file_path, "r") as consensus_seq:
        fasta_sequences = SeqIO.parse(consensus_seq, "fasta")
        for fasta in fasta_sequences:
            return str(fasta.seq)
    raise ValueError("Error: No sequences found in the FASTA file.")

def read_identity_scores(file_path):
    """Read identity scores and return them as a list."""
    with open(file_path, "r") as identity_file:
        return [float(line.strip()) for line in identity_file.readlines()]


default_fasta_dir = get_directory(directory)
fasta_file_path = get_file_path(default_fasta_dir, consensus_file)
identity_file_path = get_file_path(default_fasta_dir, identity_file)



# SECTION 2: Input validation

"""Check if primer length values are positive integers"""
if not (isinstance(min_primer_length, int) and isinstance(max_primer_length, int) and min_primer_length > 0 and max_primer_length > 0):
    raise ValueError("Error: 'min_primer_length' and 'max_primer_length' must be positive integers.")


"""Check if min_tm, max_tm, min_gc y max_gc are logical values"""
if not (min_tm < max_tm):
    raise ValueError("Error: 'min_tm' must be lower than 'max_tm'.")
if not (min_gc < max_gc):
    raise ValueError("Error: 'min_gc' must be lower than 'max_gc'.")


"""Check if percentages add to 100: """
# Primer quality percentages
primer_quality_check = gc_weight + self_comp_weight + hairpin_weight + dist_weight + gc_clamp_weight if use_gc_clamp else 0
if primer_quality_check != 100:
    raise ValueError(f"Error: Weights of primer quality properties (gc_weight, self_comp_weight, hairpin_weight, dist_weight and gc_clamp_weight) should add to 100. Currently they add to {primer_quality_check}. Please, make sure this variables are set accordingly.")
# Primer ranking percentages
primer_ranking_check = cons_weight + quality_weight
if primer_ranking_check != 100:
    raise ValueError(f"Error: Weights of primer ranking properties (cons_weight and quality_weight) should add to 100. Currently they add to {primer_ranking_check}. Please, make sure this variables are set accordingly.")
# Pair quality percentages
pair_quality_check = pair_comp_weight + tm_diff_weight + size_weight
if pair_quality_check != 100:
    raise ValueError(f"Error: Weights of primer quality properties (pair_comp_weight, tm_diff_weight and size_weight) should add to 100. Currently they add to {pair_quality_check}. Please, make sure this variables are set accordingly.")
# Pair ranking percentages
pair_ranking_check = identity_weight + pair_quality_weight
if pair_ranking_check != 100:
    raise ValueError(f"Error: Weights of primer ranking properties (cons_weight and quality_weight) should add to 100. Currently they add to {pair_ranking_check}. Please, make sure this variables are set accordingly.")


"""Check if the probe size range is within acceptable limits"""
if min_probe_size < 15 or max_probe_size > 30:
    user_confirmation = input(
        f"The specified probe size range ({min_probe_size}-{max_probe_size} nts) is outside the normal parameters for TaqMan probes. Do you want to continue? (Y/N): ").upper()
    if user_confirmation != "Y":
        print("Exiting script.")
        exit()


sequence = read_fasta(fasta_file_path)


def is_dna(seq):
    """Check if the sequence is DNA."""
    return 'T' in seq and 'U' not in seq

def is_rna(seq):
    """Check if the sequence is RNA."""
    return 'U' in seq and 'T' not in seq

def get_sequence_type(seq):
    """Determine if the sequence is DNA or RNA."""
    if is_dna(seq) and not is_rna(seq):
        return "DNA"
    elif is_rna(seq) and not is_dna(seq):
        return "RNA"
    else:
        raise ValueError("Error: Sequence is neither clearly DNA nor RNA.")


sequence_type = get_sequence_type(sequence)


# Sequence type confirmation
confirmation = input(f"The sequence appears to be {sequence_type}. Is this correct? (Y/N): ").lower()
if confirmation != "y":
    raise ValueError("Error: Please check the input sequence. Exiting.")


def reverse_complement(seq):
    """Return the reverse complement of a DNA/RNA sequence."""
    if sequence_type == "DNA":
        complement = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K',
            'S': 'S', 'W': 'W', 'B': 'V', 'D': 'H',
            'H': 'D', 'V': 'B', 'N': 'N'
        }
        return ''.join(complement.get(base, base) for base in reversed(seq))

    elif sequence_type == "RNA":
        complement = {
            'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
            'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K',
            'S': 'S', 'W': 'W', 'B': 'V', 'D': 'H',
            'H': 'D', 'V': 'B', 'N': 'N'
        }
        return ''.join(complement.get(base, base) for base in reversed(seq))


"""Setting the variables for fw/rv sequences and identity scores"""
rev_comp_sequence = reverse_complement(sequence)
identity_scores = read_identity_scores(identity_file_path)
reversed_scores = identity_scores[::-1]


# Debugging lines. Ensure the sequence is being read appropriately.
print(f"Original sequence: {sequence[:30]}...{sequence[-30:]}")
print(f"Reverse complement sequence: {rev_comp_sequence[:30]}...{rev_comp_sequence[-30:]}")
print(f"{os.path.basename(fasta_file_path)} sequence is {len(sequence)}bp long")



# SECTION 3: Primer design

def melting(seq, sequence_type):
    """Calculate Tm and GC content for a given fragment."""
    primer_len = len(seq)
    count_a = seq.count("A")
    count_t = seq.count("T") if sequence_type == "DNA" else 0
    count_g = seq.count("G")
    count_c = seq.count("C")
    count_u = seq.count("U") if sequence_type == "RNA" else 0
    count_r = seq.count("R") / 2
    count_y = seq.count("Y") / 2
    count_k = seq.count("K") / 2
    count_m = seq.count("M") / 2
    count_s = seq.count("S") / 2
    count_w = seq.count("W") / 2
    count_b = seq.count("B") / 3
    count_d = seq.count("D") / 3
    count_h = seq.count("H") / 3
    count_v = seq.count("V") / 3
    count_n = seq.count("N") / 4

    count_gc = count_g + count_c + count_s + count_b + count_v
    count_at = count_a + (count_t if sequence_type == "DNA" else 0) + count_w + count_d + count_h + (count_u if sequence_type == "RNA" else 0)
    count_other = count_r + count_y + count_k + count_m + count_n

    # Calculate GC-content (%).
    gc = (count_gc / primer_len) * 100

    # Calculate Tm.
    if primer_len >= 13:
        tm = 64.9 + 41 * (count_gc - 16.4) / (count_gc + count_at + count_other)
    else:
        tm = (2 * count_at) + (4 * count_gc)

    return str("%.2f" % tm), str("%.2f" % gc)


def calculate_gc_content(seq, window_size):
    """Calculate the GC-content in sliding-windows throughout the sequence of each primer."""
    gc_content = []
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i + window_size]
        gc_count = window.count('G') + window.count('C')
        gc_percentage = (gc_count / window_size) * 100
        gc_content.append(gc_percentage)
    return gc_content


def evaluate_gc_distribution(gc_distribution, window_size):
    """Evaluate the uniformity of the GC distribution."""
    iqr_values = []
    for i in range(len(gc_distribution) - window_size + 1):
        window = gc_distribution[i:i + window_size]
        q1 = sorted(window)[len(window) // 4]
        q3 = sorted(window)[len(window) * 3 // 4]
        iqr = q3 - q1
        iqr_values.append(iqr)

    if len(iqr_values) == 0:
        return float('inf'), float('inf')
    normalized_iqr = sum(iqr_values) / len(iqr_values) / window_size

    max_variation_gc = [0, 1] * (window_size // 2)
    if window_size % 2 != 0:
        max_variation_gc.append(0)
    max_iqr_values = []
    for i in range(len(gc_distribution) - window_size + 1):
        max_variation_gc = [(100 * x / (window_size - 1)) for x in range(window_size)]
        window = max_variation_gc
        q1 = sorted(window)[len(window) // 4]
        q3 = sorted(window)[len(window) * 3 // 4]
        iqr = q3 - q1
        max_iqr_values.append(iqr)
    max_iqr = max(max_iqr_values) / window_size

    return normalized_iqr, max_iqr


def has_gc_clamp(seq, clamp_length):
    """Checks if the primer has a GC-clamp."""
    clamp_score = 0
    if seq[-clamp_length:].count('G') + seq[-clamp_length:].count('C') >= clamp_length:
        clamp_score = gc_clamp_weight
    return clamp_score


def calculate_dimer_stability(seq):
    """Calculate the stability of auto-dimers."""
    max_stability = 0
    length = len(seq)
    for i in range(length):
        for j in range(i + 1, length):
            stability = 0
            for k in range(min(length - j, j - i)):
                if seq[i + k] == reverse_complement(seq[j + k]):
                    stability += 1
                elif seq[i + k] == 'N' or seq[j + k] == 'N':
                    stability += 0.25
                else:
                    stability -= 0.5
            max_stability = max(max_stability, stability)

    # Normalize to percentage (0% - 100%)
    if max_stability == 0:
        return 0.0
    else:
        percent_stability = max_stability / length
        return percent_stability


def penalize_hairpins(primer):
    """Penalize the possible formation of hairpins loops."""
    hairpin_penalty = 0
    max_hairpin_penalty = 0
    length = len(primer)
    for i in range(length - 2):
        for j in range(i + 3, length):
            # Max penalty for each potential hairpin
            max_hairpin_penalty += 2
            if primer[i] == 'C' and primer[i + 1] == 'G':
                max_hairpin_penalty += 1  # Extra penalty for CG pairs

            # Actual penalty for the primer
            if primer[i] == reverse_complement(primer[j]) and primer[i + 1] == reverse_complement(primer[j - 1]):
                hairpin_penalty += 2
            if primer[i] == 'C' and primer[i + 1] == 'G':
                hairpin_penalty += 1  # Extra penalty for CG pairs

    return hairpin_penalty, max_hairpin_penalty


def calculate_gc_deviation(gc):
    gc_deviation = 0
    fragment_gc = float(gc)
    if min_gc <= fragment_gc <= max_gc:
        gc_deviation = gc_weight
    elif float(gc) <= min_gc or max_gc <= float(gc):
        gc_deviation = abs(float(gc) - ((max_gc + min_gc) / 2))

    return gc_deviation


def design_primers(seq, scores, min_primer_length, max_primer_length, sequence_type, use_gc_clamp):
    """Create and score primers."""
    primers = []
    # Design every possible primer in size range.
    for primer_size in range(min_primer_length, max_primer_length + 1):
        for i in range(len(scores) - primer_size + 1):
            fragment = seq[i:i + primer_size]
            tm, gc = melting(fragment, sequence_type)
            if min_tm <= float(tm) <= max_tm:
                score = sum(scores[i:i + primer_size])
                start = i + 1
                end = i + primer_size
                abs_score = score
                rel_score = (score / primer_size) * 100

                # Score every primer.
                self_comp = calculate_dimer_stability(fragment)
                self_comp_score = (1 - self_comp) * self_comp_weight

                hairpin_penalty = penalize_hairpins(fragment)
                hairpin_penalty_score = (hairpin_penalty[0] / hairpin_penalty[1])
                hairpin_score = (1 - hairpin_penalty_score) * hairpin_weight

                gc_score = calculate_gc_deviation(gc)

                gc_cont = calculate_gc_content(fragment, window_size)
                gc_dist = evaluate_gc_distribution(gc_cont, window_size)
                gc_dist_score = (1 - (gc_dist[0] / gc_dist[1])) * dist_weight

                gc_clamp_score = has_gc_clamp(fragment, clamp_length) if use_gc_clamp else 0
                check_gc_clamp = "Yes" if gc_clamp_score != 0 else "No"

                # Calculate primer quality.
                primer_quality = gc_score + self_comp_score + gc_dist_score + gc_clamp_score + hairpin_score

                # Calculate primer score (ranking position).
                primer_score = (rel_score * (cons_weight / 100)) + (primer_quality * (quality_weight / 100))


                primers.append([
                    round(primer_score, 2), round(abs_score, 2), len(fragment), round(rel_score, 2), round(primer_quality, 2), start, fragment, end, tm, gc, round(gc_score, 2), round(self_comp_score, 2), check_gc_clamp, round(hairpin_score, 2), round(gc_dist_score, 2)
                ])
    return primers


forward_primers = design_primers(sequence, identity_scores, min_primer_length, max_primer_length, sequence_type, use_gc_clamp)
reverse_primers = design_primers(rev_comp_sequence, reversed_scores, min_primer_length, max_primer_length, sequence_type, use_gc_clamp)


forward_primers.sort(key=lambda x: (x[0], x[3], x[4]), reverse=True)
reverse_primers.sort(key=lambda x: (x[0], x[3], x[4]), reverse=True)


# SECTION 4: Design primer pairs.


def calculate_pair_comp(fwd_primer, rev_primer):
    """Penalize the formation of primer-dimers."""
    max_score = 0
    rev_primer_rc = reverse_complement(rev_primer)
    length_fwd = len(fwd_primer)
    length_rev = len(rev_primer_rc)
    for i in range(length_fwd):
        for j in range(length_rev):
            score = 0
            for k in range(min(length_fwd - i, length_rev - j)):
                if sequence_type == "DNA":
                    if fwd_primer[i + k] == rev_primer_rc[j + k]:
                        score += 1
                    elif fwd_primer[i + k] in "ATGCN" and rev_primer_rc[j + k] in "ATGCN":
                        score -= 0.5
                elif sequence_type == "RNA":
                    if fwd_primer[i + k] == rev_primer_rc[j + k]:
                        score += 1
                    elif fwd_primer[i + k] in "AUGCN" and rev_primer_rc[j + k] in "AUGCN":
                        score -= 0.5
            if score > max_score:
                max_score = score
    if max_score == 0:
        return 0.0
    else:
        complementarity = (max_score / max(length_fwd, length_rev))
        return complementarity


def get_primer_pairs(forward_primers, reverse_primers, min_size, max_size):
    """Generate primer pairs within the product size range."""
    primer_pairs = []
    for i, f_primer in enumerate(forward_primers):
        fw_tm = f_primer[8]
        for j, r_primer in enumerate(reverse_primers):
            rv_tm = r_primer[8]
            if abs(float(fw_tm) - float(rv_tm)) > pair_tm_deviation:
                continue
            if r_primer[7] > f_primer[5]:
                product_size = r_primer[7] - f_primer[5] + 1
                if min_size <= product_size <= max_size:
                    if r_primer[7] > f_primer[5] + 1:

                        # Pair properties
                        combined_identity = ((f_primer[3] + r_primer[3]) / 2)

                        amplicon_range = f"{f_primer[5]}->{r_primer[7]}"

                        pair_comp = calculate_pair_comp(f_primer[6], r_primer[6])
                        pair_comp_score = (1 - pair_comp) * pair_comp_weight

                        tm_diff = abs(float(f_primer[8]) - float(r_primer[8]))
                        tm_diff_score = (1 - tm_diff / pair_tm_deviation) * tm_diff_weight

                        size_score = (1 - (abs(product_size - min_size) / (max_size - min_size))) * size_weight

                        # Calculate pair quality.
                        pair_quality = pair_comp_score + tm_diff_score + size_score

                        # Calculate pair score.
                        pair_score = (combined_identity * (identity_weight / 100)) + (pair_quality * (pair_quality_weight / 100))

                        primer_pairs.append([
                            "", round(pair_score, 2), round(combined_identity, 2), round(pair_quality, 2),
                            round(pair_comp_score, 2), round(tm_diff, 2), round(size_score, 2),
                            len(f_primer[6]), fw_tm, f_primer[3], f_primer[6],
                            len(r_primer[6]), rv_tm, r_primer[3], r_primer[6],
                            product_size, amplicon_range, i + 1, j + 1  # Fw and Rv primer index
                        ])

    primer_pairs.sort(key=lambda x: (x[1], x[2], x[3]), reverse=True)

    # Assign pair names based on ranking.
    for idx, pair in enumerate(primer_pairs, start=1):
        pair_name = f"pair{idx}_fw{pair[-2]}_rv{pair[-1]}"
        pair[0] = pair_name
    # Remove unwanted primer index.
    for pair in primer_pairs:
        pair.pop(-1)
        pair.pop(-1)

    return primer_pairs


primer_pairs = get_primer_pairs(forward_primers, reverse_primers, min_size, max_size)



# SECTION 5: Write "primer and pair" output.

def write_output(file_path, fwd_primers, rev_primers, primer_pairs, fasta_filename):
    """Write the primers to a CSV file."""

    with open(file_path, "w") as output_file:

        # Write headers
        output_file.write(
            "Fw_Name\tTotal_Score (%)\tCons. (abs)\tCons. (%)\tQuality (%)\tStart\tFw_Seq\tEnd\tTm (ºC)\tGC (%)\tGC_Score\tSelf-Comp.\tCG_Clamp\tHairpin_score\tGC-Dist.\t\t\t" +
            "Rv_Name\tTotal_Score (%)\tCons. (abs)\tCons. (%)\tQuality (%)\tStart\tRv_Seq\tEnd\tTm (ºC)\tGC (%)\tGC_Score\tSelf-Comp.\tCG_Clamp\tHairpin_score\tGC-Dist.\t\t\t" +
            "Pair_Name\tTotal_Score (%)\tPair_Cons. (%)\tQuality (%)\tPair-Comp\tPair_Tm (ºC avg)\tTm_difference (ºC)\tSize_score\tFw_size\tFw_Tm (ºC)\tFw_Cons. (%)\tFw_Seq\tRv_size\tRv_Tm (ºC)\tRv_Cons. (%)\tRv_Seq\tLength\tAmplicon_Range\n"
        )

        # Write data for primers and pairs
        for idx, (f_primer, r_primer, pair) in enumerate(zip(fwd_primers, rev_primers, primer_pairs), start=1):

            f_primer_name = f"{fasta_filename}_fw{idx}"
            f_abs_size = f"{f_primer[1]}/{f_primer[2]}\t"

            r_primer_name = f"{fasta_filename}_rv{idx}"
            r_abs_size = f"{r_primer[1]}/{r_primer[2]}\t"

            pair_name = f"{pair[0]}"
            pair_avg_tm = f"{round(((float(pair[8]) + float(pair[12])) / 2), 2)}"

            output_file.write(
                f"{f_primer_name}\t" + f"{f_primer[0]}\t" + f"{f_abs_size}" + "\t".join(map(str, f_primer[3:15])) + "\t\t\t" +
                f"{r_primer_name}\t" + f"{r_primer[0]}\t" + f"{r_abs_size}" + "\t".join(map(str, r_primer[3:15])) + "\t\t\t" +
                f"{pair_name}\t" + f"{pair[1]}\t" + f"{pair[2]}\t" + f"{pair[3]}\t" + f"{pair[4]}\t" + f"{pair_avg_tm}\t" + "\t".join(map(str, pair[5:17])) + "\n"
            )


output_file_path = os.path.splitext(fasta_file_path)[0] + "_primers.csv"
write_output(output_file_path, forward_primers, reverse_primers, primer_pairs, os.path.basename(fasta_file_path).split('.')[0])

print(f"Primers created successfully and written to {os.path.basename(output_file_path)}")



# SECTION 6: TaqMan probes design.

def calculate_probe_conservation(probe_sequence, identity_scores, start_index):
    """Calculate conservation for any probe."""
    probe_len = len(probe_sequence)
    conservation_score = sum(identity_scores[start_index:start_index + probe_len])
    return conservation_score


def taqman_size(sequence, min_probe_size, max_probe_size, ideal_tm):
    """Finds the best probe with Tm closest to the ideal Tm."""
    best_probe_size = None
    closest_tm_diff = float('inf')
    for probe_size in range(min_probe_size, max_probe_size + 1):
        for i in range(len(sequence) - probe_size + 1):
            probe_sequence = sequence[i:i + probe_size]
            tm, gc = melting(probe_sequence, sequence_type)
            tm_diff = abs(float(tm) - ideal_tm)
            if tm_diff < closest_tm_diff:
                closest_tm_diff = tm_diff
                best_probe_size = probe_size

    return best_probe_size


def design_taqman_probes_for_pairs(sequence, identity_scores, primer_pairs, min_size, max_size, num_pairs, max_primer_length, min_distance_between_probes, min_conservation=0):
    """Designs TaqMan probes based on a limited number of primer pairs and desired product size range."""
    taqman_probes = []
    probe_counter = 1  # Initialize probe_counter outside the loop
    for idx, pair in enumerate(primer_pairs[:num_pairs]):
        product_length = int(pair[15])
        start_index = int(pair[16].split("->")[0])
        end_index = int(pair[16].split("->")[1])
        amplicon_sequence = sequence[start_index:end_index]
        amplicon_length = len(amplicon_sequence)
        if min_size <= product_length <= max_size:
            probe_name_base = f"{pair[0]}_probe"

            # Calculate the average Tm for the primer pair.
            forward_tm = float(pair[8])
            reverse_tm = float(pair[12])
            pair_avg_tm = (forward_tm + reverse_tm) / 2

            # Exclude primer regions.
            adjusted_amplicon_sequence = amplicon_sequence[max_primer_length:-max_primer_length]
            adjusted_identity_scores = identity_scores[start_index + max_primer_length:end_index - max_primer_length]
            adjusted_amplicon_length = len(adjusted_amplicon_sequence)

            # Determine the ideal Tm for the TaqMan probe.
            ideal_tm = pair_avg_tm + ideal_tm_increase

            # Find the best probe size.
            best_probe_size = taqman_size(adjusted_amplicon_sequence, min_probe_size, max_probe_size, ideal_tm)

            if best_probe_size is None:
                continue

            best_probes = []

            # Generate potential probe sequences.
            i = 0
            while len(best_probes) < num_probes_per_pair and i + best_probe_size <= adjusted_amplicon_length:
                probe_sequence = adjusted_amplicon_sequence[i:i + best_probe_size]
                abs_conservation = calculate_probe_conservation(probe_sequence, adjusted_identity_scores, i)

                # Select the probe with the highest conservation.
                if abs_conservation >= min_conservation:
                    tm, gc = melting(probe_sequence, sequence_type)

                    # Check if probe is too close to others previously selected.
                    is_valid_probe = True
                    for selected_probe in best_probes:
                        if abs(selected_probe[1] - (start_index + max_primer_length + i)) <= min_distance_between_probes:
                            is_valid_probe = False
                            break
                    if is_valid_probe:
                        best_probes.append((
                            probe_sequence, start_index + max_primer_length + i,
                            round(abs_conservation, 2), tm, gc, best_probe_size,
                            pair_avg_tm
                        ))

                i += 1

            # Sort probes and select the best num_probes_per_pair probes
            best_probes.sort(key=lambda x: (x[2], x[3]), reverse=True)
            best_probes = best_probes[:num_probes_per_pair]

            for probe_sequence, start_index, abs_conservations, tm, gc, probe_size, pair_avg_tm in best_probes:
                probe_name = f"{probe_name_base}_{probe_counter}"
                taqman_probes.append([
                    probe_name, abs_conservations, round((abs_conservations / probe_size) * 100, 2), start_index + 1, probe_sequence, start_index + probe_size, tm, pair_avg_tm, gc, probe_size,
                ])
                probe_counter += 1

            # Blank row to separate groups of possible probes for different primer pairs
            taqman_probes.append([""] * 10)

    return taqman_probes


taqman_probes = design_taqman_probes_for_pairs(sequence, identity_scores, primer_pairs, min_size, max_size, num_primer_pairs, max_primer_length, min_distance_between_probes)



# SECTION 7: Write TaqMan output

output_file_path_taqman = os.path.splitext(fasta_file_path)[0] + "_TaqMan.csv"

with open(output_file_path_taqman, "w") as output_taqman_file:
    # Write the header
    output_taqman_file.write(
        "Probe_Name\tCons. (abs)\tCons. (%)\tStart\tProbe_Sequence\tEnd\tProbe Tm (ºC)\tPair_Tm (ºC avg)\tGC%\n"
    )

    # Write data for TaqMan probes
    for probe in taqman_probes:
        # Add the size of the probe after the absolute conservation value
        abs_conservation_with_size = f"{probe[1]}/{probe[9]}\t"

        output_taqman_file.write(
            f"{probe[0]}\t" + f"{abs_conservation_with_size}" + "\t".join(map(str, probe[2:9])) + "\n"
        )

    print(f"TaqMan probes successfully designed and written to {os.path.basename(output_file_path_taqman)}.")
