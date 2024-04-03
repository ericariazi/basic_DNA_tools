# Basic DNA tools
"""
The Basic DNA Tools script carries out common tasks performed by molecular biologists
and bioinformaticians for manipulating DNA sequences, such as constructing the 
reverse complement, or generating all possible polypeptides from a DNA sequence.

The script runs in two modes, Screen Mode and File Mode. 

In Screen Mode, the script takes a single DNA sequence and optionally a name. By default, output is sent stdout. But results can be optionally sent to output CSV file if the path is provided.

In File Mode, the script takes an input FASTA or FASTQ file as a command-line argument, and optionally a user designated CSV file.

FASTA and FASTQ files must be in the proper format. FASTA file extensions must be ".fasta", ".fas", ".fa", ".fna", or ".ffn". FASTQ file extensions must be ".fastq".
"""

import sys
import argparse
import csv
import os.path

class DNA:
    """Represents a DNA sequence with optional name and quality attributes."""

    def __init__(self, sequence, name="", quality=""):
        """
        Initialize a new DNA instance.
        
        :param sequence: str, DNA nucleotide sequence
        :param name: str, optional, name of the DNA sequence
        :param quality: str, optional, quality scores of the DNA sequence
        """   
        self.seq = sequence # This will also set self._seq and validate the sequence
        self.name = name
        self.quality = quality
    
    def __str__(self):
        """
        Returns the string representation of the DNA object.

        The output format varies depending on whether the DNA object has a name and/or quality scores defined.
        """
        if self.name and self.quality: # If both name and quality scores are available, include them in the output
            return f"DNA name: {self.name}\nDNA sequence: {self.seq}\nLength: {self.length}\nQuality scores: {self.quality}"
        elif self.name and not self.quality: # If only the name is available, include it and the sequence in the output
            return f"DNA name: {self.name}\nDNA sequence: {self.seq}\nLength: {self.length}"
        else: # If neither name nor quality scores are available, only include the sequence in the output
            return f"DNA sequence: {self.seq}\nLength: {self.length}"

    @property
    def seq(self):
        return self._seq
    
    @seq.setter
    def seq(self, sequence):
        sequence = sequence.upper()
        for base in sequence:
            if base not in ["A", "T", "G", "C"]:
                raise ValueError(f"Invalid base {base} in sequence")
        self._seq = sequence
        # Now that the nucleotides are validated, update the length
        self.length = len(sequence)
    

def main():
    """
    Entry point of the script. Parses command-line arguments to run in either screen mode or file mode.
    
    In screen mode, a DNA sequence is inputted directly via command-line arguments along with optional sequence name.
    In file mode, DNA sequences are read from a specified FASTA or FASTQ file. Outputs can be directed to a CSV file.
    """
    parser = argparse.ArgumentParser(description="Provides molecular biology functions for DNA molecules")
    # Create a mutually exclusive groups of args such that inputting the DNA sequence occurs via the command-line or a file
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-s", "--seq", help="Screen mode. Input 1 DNA sequence. Default output is standard out if no output path provided.)")
    group.add_argument("-f", "--file", help="File mode. Path to input FASTA or FASTQ file. The filename must have the appropriate extension. Output is written to a CSV file.")
    parser.add_argument("-n", "--name", help = "Optional name of inputted sequence in screen mode. Not used in File mode.")
    parser.add_argument("-o", "--output", help="Path to output CSV file. If not provided while in File mode, output is written to 'output.csv'.")
    args = parser.parse_args()
    # Check if --name is provided without --seq
    if args.name and not args.seq:
        parser.error("--name argument is only valid when provided in screen mode with --seq")
    if args.seq:  # Screen mode
        run_screen_mode(args.seq, args.name)
    elif args.file:  # File mode
        run_file_mode(args.file, args.output)


def run_screen_mode(sequence, name, output_path=None):
    """
    Handles processing of a DNA sequence provided via command-line arguments in screen mode.
    
    :param sequence: str, DNA nucleotide sequence
    :param name: str, optional name of the DNA sequence
    :param output_path: str, optional path to write the output to a CSV file. Prints to stdout if None.
    """
    print(f"\nRunning screen mode\n")
    try:
        dna = DNA(sequence, name)
        if dna:
            if output_path:
                write_functions(dna, output_path)
            else:
                print_functions(dna)
    except ValueError as e:
        sys.exit(e)


def run_file_mode(file, output_path=None):
    """
    Processes DNA sequences from a specified FASTA or FASTQ file in file mode.
    
    :param file: str, path to the input FASTA or FASTQ file
    :param output_path: str, optional path to the output CSV file. Uses 'output.csv' as default if not provided.
    """
    print(f"\nRunning file mode")  
    if not os.path.isfile(file):
        sys.exit("File not found.")
    # Logic to determine if an output path was provided, and if not, set the default output path 
    if not output_path:
        output_path = "./output.csv"
    # Input file extension check
    fasta_extensions = (".fasta", ".fas", ".fa", ".fna", ".ffn")
    fastq_extension = (".fastq")
    if not file.endswith(fasta_extensions) and not file.endswith(fastq_extension):
        sys.exit("Filename does not end in a FASTA or FASTQ file extension.")  
    if file.endswith(fasta_extensions):
        file_type = "FASTA"
    elif file.endswith(fastq_extension):
        file_type = "FASTQ"
    with open(file, "r") as f:
        record_count = 0  # Number of processed sequence records 
        while True:
            # Skip any empty lines
            header = ""
            while not header:
                line = f.readline() # Do not strip to allow for checking if this is the End of File
                if not line:  # Check for End of file wihtout stripping, will need to break out of 2 loops
                    break  # Exit first loop (inner)
                header = line.strip()
            if not line: # Need to check again to exit second loop (outer)
                break
            name = header[1:] # Discard symbol from header but keep the rest
            # Input file format validation 
            if file_type == "FASTA":
                if not header.startswith(">"):
                    sys.exit("Incorrect format. FASTA record does not start with '>'")
                seq = process_FASTA_record(f)
                quality = ""  # Default value in case of FASTA file        
            elif file_type == "FASTQ":
                if not header.startswith("@"):
                    sys.exit("Incorrect format. FASTQ record does not start with '@'")
                seq, quality = process_FASTQ_record(f, name)
            # Process the DNA sequence
            try:
                dna = DNA(seq, name, quality)
                write_functions(dna, output_path)
                record_count += 1  # Increment record count after successful processing
            except ValueError:
                sys.exit(f"Invalid base found while processing record # {record_count}, sequence name: {name}")  
        print(f"Processed {record_count} records.")  # Output the number of processed records


def process_FASTA_record(f):
    """
    Extracts a DNA sequence from a FASTA file format record.
    
    :param f: file object, input file opened in read mode
    :return: str, extracted DNA sequence
    """ 
    seq = ""  # Initialize empty string for accumulating sequence data.
    # Read and accumulate the sequence until the next symbol separator, empty line or end of file is found
    while True:
        position = f.tell()  # Get current file location
        line = f.readline().strip()
        if line.startswith(">") or not line:  # There is no more sequence data to accumulate for this record.
            f.seek(position)  # Rewind if the line starts a new record or EOF
            break
        else: # Accumulate sequence lines
            seq += line
    return seq


def process_FASTQ_record(f, name):
    """
    Extracts a DNA sequence and its corresponding quality scores from a FASTQ file format record.
    
    :param f: file object, input file opened in read mode
    :param name: str, name of the sequence being processed (used for error messages)
    :return: tuple, containing the extracted DNA sequence and quality scores
    """  
    seq, quality = "", ""  # Initialize empty strings for accumulating sequence and quality data.
    # Read and accumulate the sequence until the '+' separator is found
    while True:
        line = f.readline().strip()
        if line.startswith("+"):  # Accumulate sequence lines until "+" encountered. This is a comment or blank line
            break
        else: 
            seq += line
    # Immediately after "+" line, start accumulating quality scores until the lengths of the q scores and sequence are the same
    while len(quality) < len(seq):  
        line = f.readline().strip()
        quality += line
    #validate sequence and quality lengths match
    if len(seq) != len(quality):
        sys.exit(f"Length mismatch between sequence and quality scores in record: {name}")
    return seq, quality


def get_complement(dna):
    """
    Constructs the complement of a DNA sequence using a dictionary to map and replace each nucleotide 
    in the input DNA sequence with its complementary base.

    :param dna: An instance of the DNA class containing the DNA sequence.
    :type dna: DNA 
    :return: The constructed complement of the DNA sequence.
    :rtype: str
    """
    # use a dictionary to map complementary base replacements
    complement = {"G": "C", "A": "T", "T": "A", "C": "G"}
    # use list comprehension to loop through the sequence,look up base in complement dict for replacement, and concatenate list
    return  "".join(complement[base] for base in dna.seq)


def get_reverse_complement(dna):
    """
    Generates the reverse complement of a given DNA sequence. This involves creating the complement of the sequence and then reversing it.

    :param dna: An instance of the DNA class containing the DNA sequence.
    :type dna: DNA 
    :return: The reverse complement of the DNA sequence.
    :rtype: str
    """
    # Utilizes get_complement to create the complement and then reverses it
    return get_complement(dna)[::-1]  


def calculate_GC_pct(dna):
    """
    Computes the GC content percentage of a given DNA sequence. 

    :param dna: An instance of the DNA class containing the DNA sequence to analyze.
    :type dna: DNA 
    :return: The percentage of G and C nucleotides in the DNA sequence as a float.
    :rtype: float
    """
    G_count = dna.seq.count("G")
    C_count = dna.seq.count("C")
    return ((G_count + C_count)/len(dna.seq)) * 100

 
def calculate_AT_pct(dna):
    """
    Computes the AT content percentage of a given DNA sequence. 

    :param dna: An instance of the DNA class containing the DNA sequence to analyze.
    :type dna: DNA 
    :return: The percentage of A and T nucleotides in the DNA sequence as a float.
    :rtype: float
    """
    A_count = dna.seq.count("A")
    T_count = dna.seq.count("T")
    return ((A_count + T_count)/len(dna.seq)) * 100


def methyl_conversion(dna):
    """
    Simulates the conversion of unmethylated cytosine (C) to thymine (T) in a DNA sequence, representing 
    the process occurring during bisulfite treatment or enzymatic methyl conversion reactions. This conversion
    is typically used in DNA methylation studies to identify unmethylated cytosines.

    :param dna: The DNA sequence to be converted.
    :type dna: DNA
    :return: A string representing the DNA sequence after conversion, with converted cytosines represented as lowercase 't'.
    :rtype: str
    """
    return dna.seq.replace("C", "t")


def transcribe(dna):
    """
    Simulates the transcription of a DNA sequence into an RNA sequence by replacing T with U.

    :param dna: The DNA sequence to be transcribed.
    :type dna: DNA
    :return: The RNA sequence as a string resulting from the transcription process.
    :rtype: str
    """
    return dna.seq.replace("T", "U")


def translate(dna,letters=1):
    """
    Simulates the translation of a DNA sequence into polypeptide sequences for all 6 reading frames,
    including three frames from the original sequence and three from its reverse complement. 
    This is achieved by mapping DNA codons to amino acids according to a translation table.

    :param dna: The DNA sequence to be transcribed.
    :type dna: DNA
    :param letters: Specifies the amino acid representation format: 1 for single-letter codes (default),
                    and 3 for three-letter codes.
    :type letters: int
    :return: A dictionary containing all six reading frames with amino acid sequences resulting from the translation.
             Keys indicate the strand and reading frame, while values are the corresponding amino acid sequences.
    :rtype: dict
    """
    codon_table = {
        "TTT":("Phe","F"), "TTC":("Phe","F"), "TTA":("Leu","L"), "TTG":("Leu","L"),
        "CTT":("Leu","L"), "CTC":("Leu","L"), "CTA":("Leu","L"), "CTG":("Leu","L"),
        "ATT":("Ile","L"), "ATC":("Ile","L"), "ATA":("Ile","L"), "ATG":("Met","M"),
        "GTT":("Val","V"), "GTC":("Val","V"), "GTA":("Val","V"), "GTG":("Val","V"),

        "TCT":("Ser","S"), "TCC":("Ser","S"), "TCA":("Ser","S"), "TCG":("Ser","S"),
        "CCT":("Pro","P"), "CCC":("Pro","P"), "CCA":("Pro","P"), "CCG":("Pro","P"),
        "ACT":("Thr","T"), "ACC":("Thr","T"), "ACA":("Thr","T"), "ACG":("Thr","T"),
        "GCT":("Ala","A"), "GCC":("Ala","A"), "GCA":("Ala","A"), "GCG":("Ala","A"),

        "TAT":("Tyr","Y"), "TAC":("Tyr","Y"), "TAA":("*","*"), "TAG":("*","*"),
        "CAT":("His","H"), "CAC":("His","H"), "CAA":("Gln","Q"), "CAG":("Gln","Q"),
        "AAT":("Asn","N"), "AAC":("Asn","N"), "AAA":("Lys","K"), "AAG":("Lys","K"),
        "GAT":("Asp","D"), "GAC":("Asp","D"), "GAA":("Glu","E"), "GAG":("Glu","E"),

        "TGT":("Cys","C"), "TGC":("Cys","C"), "TGA":("*","*"), "TGG":("Trp","W"),
        "CGT":("Arg","R"), "CGC":("Arg","R"), "CGA":("Arg","R"), "CGG":("Arg","R"),
        "AGT":("Ser","S"), "AGC":("Ser","S"), "AGA":("Arg","R"), "AGG":("Arg","R"),
        "GGT":("Gly","G"), "GGC":("Gly","G"), "GGA":("Gly","G"), "GGG":("Gly","G"),
    }
    proteins = {}
    reverse_complement = get_reverse_complement(dna)
    for j in range(0,6):  # Loop through top strand (0<j<3) and then bottom strand (3<j<6)
        protein = []  # collect all 6 reading frames in protein dictionary
        for i in range(0, len(dna.seq), 3): # Step through sequence 3 bases at a time
            if j < 3: # Loop through top strand
                codon = dna.seq[i+j:i+3+j] # i indexes for reading frame and j for strand
            else: # Loop through bottom strand
                codon = reverse_complement[i+j-3:i+3+j-3]  #search the bottom strand (reverse_complement)
            if codon in codon_table:
                # Decide whether to append the full name or the single-letter code based on `letters`
                if letters == 3:
                    protein.append(codon_table[codon][0]) # Three-letter code
                    if codon_table[codon][0] == "*":
                        break
                elif letters == 1:
                    protein.append(codon_table[codon][1]) # Single-letter code
                    if codon_table[codon][0] == "*":
                        break
        if letters == 3:
            if j < 3:
                proteins["Top strand, Reading frame " + str(j+1) + " (3-letter code)"] = "".join(protein)
            else:
                proteins["Bottom strand, Reading frame " + str(j+1) + " (3-letter code)"] = "".join(protein)
        elif letters == 1:
            if j < 3: 
                proteins["Top strand, Reading frame " + str(j+1) + " (1-letter code)"] = "".join(protein)
            else:
                proteins["Bottom strand, Reading frame " + str(j+1) + " (1-letter code)"] = "".join(protein)
    return proteins
    

def print_functions(dna):
    """
    Prints various properties and translations of the given DNA object to stdout.
    
    :param dna: DNA, a DNA object for which properties are calculated and printed.
    """
    print(dna)
    print(f"Complement: {get_complement(dna)}")   
    print(f"Reverse complement: {get_reverse_complement(dna)}")
    print(f"\nGC content: {calculate_GC_pct(dna):.2f}")
    print(f"AT content: {calculate_AT_pct(dna):.2f}")
    print(f"\nmethylation conversion: {methyl_conversion(dna)}")
    print(f"RNA transcription: {transcribe(dna)}")
    print(f"\n3-letter code translation:")
    proteins_3_letter_code = (translate(dna, letters=3))
    protein_printer(proteins_3_letter_code)
    print(f"\n1-letter code translation:")
    proteins_1_letter_code = (translate(dna))
    protein_printer(proteins_1_letter_code)


def protein_printer(proteins_dict):
    """
    Prints the protein sequences contained in the provided dictionary, 
    indicating the corresponding reading frames and code types.
    
    :param proteins_dict: dict, dictionary containing protein sequences with reading frame and code type as keys
    """
    for reading_frame in proteins_dict:
        print(f"{reading_frame}: {proteins_dict[reading_frame]}")


def write_functions(dna, file_name):
    """
    Writes various properties and analyses of a DNA sequence to a CSV file.

    This function takes a DNA object and a file name as input. It then calculates various properties of the DNA sequence,
    including its length, complement, reverse complement, GC content percentage, AT content percentage, result of methylation
    conversion, and RNA transcription. It also includes the translation of the DNA sequence into polypeptide sequences
    for all 6 reading frames, both in 3-letter and 1-letter amino acid codes. These properties are written to the specified
    CSV file, with a header row defining each column if the file does not already exist.

    :param dna: The DNA object containing the sequence and optionally name and quality scores.
    :type dna: DNA
    :param file_name: Path to the CSV file where the DNA sequence properties and analyses are to be written.
    :type file_name: str
    """
    header = [
        "Name", "Sequence", "Quality Scores", "Length", "Complement", "Reverse Complement",
        "GC %", "AT %", "Methylation Conversion", "RNA transcription",
        "Top strand, Reading frame 1 (3-letter code)", "Top strand, Reading frame 2 (3-letter code)", "Top strand, Reading frame 3 (3-letter code)",
        "Bottom strand, Reading frame 4 (3-letter code)", "Bottom strand, Reading frame 5 (3-letter code)", "Bottom strand, Reading frame 6 (3-letter code)",
        "Top strand, Reading frame 1 (1-letter code)", "Top strand, Reading frame 2 (1-letter code)", "Top strand, Reading frame 3 (1-letter code)",
        "Bottom strand, Reading frame 4 (1-letter code)", "Bottom strand, Reading frame 5 (1-letter code)", "Bottom strand, Reading frame 6 (1-letter code)"
    ]
    file_exists = os.path.isfile(file_name) and os.path.getsize(file_name) > 0
    with open(file_name, "a") as o:
        writer = csv.DictWriter(o, fieldnames=header)
        if not file_exists:
            writer.writeheader()        
        row_data = {
            "Name": dna.name,
            "Sequence": dna.seq,
            "Quality Scores": dna.quality,
            "Length": dna.length,
            "Complement": get_complement(dna),
            "Reverse Complement": get_reverse_complement(dna),
            "GC %": f"{calculate_GC_pct(dna):.2f}",
            "AT %": f"{calculate_AT_pct(dna):.2f}",
            "Methylation Conversion": methyl_conversion(dna),
            "RNA transcription": transcribe(dna),
            }
        row_data.update(translate(dna, letters=3))
        row_data.update(translate(dna))
        writer.writerow(row_data)


if __name__ == "__main__":
    main()