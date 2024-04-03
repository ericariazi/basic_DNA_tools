# Basic DNA Tools

## Video Demo: <https://youtu.be/pLDzBPYO4zo>

## Description:
The Basic DNA Tools script carries out common tasks performed by molecular biologists and bioinformaticians for manipulating DNA sequences, such as constructing the reverse complement, or generating all possible polypeptides from a DNA sequence.

The script runs in two modes, Screen Mode and File Mode.

In Screen Mode, the script takes a single DNA sequence and optionally a name. By default, output is sent stdout. But results can be optionally sent to output CSV file if the path is provided.

In File Mode, the script takes an input FASTA or FASTQ file as a command-line argument, and optionally a user designated CSV file.

FASTA and FASTQ files must be in the proper format. FASTA file extensions must be ".fasta", ".fas", ".fa", ".fna", or ".ffn". FASTQ file extensions must be ".fastq".

All script functions are performed on a given inputted DNA sequence.

### DNA class
The script implements a custom class termed "DNA". DNA objects take a sequence, and upon initialization, validate that the sequence only contains the 4 nucleotides 'A', 'T', 'G', and 'C'. It then stores as attributes the DNA sequence, length, and optionally its name and quality scores. A DNA sequence's name is provided as an optional argument or extracted from the header in FASTA/FASTQ files. A DNA sequnce's quality scores are extracted from FASTQ files.

The DNA class includes a string representation of DNA objects. The output format varies depending on whether the DNA object has a name and/or quality scores defined.

#### DNA class attributes:
`DNA.seq` = DNA nucleotide sequence

`DNA.name` = Name or description of DNA sequence

`DNA.length` = Length of DNA sequence

`DNA.quality` = DNA nucleotide quality scores

### Functions
`main()`:  Entry point of the script. Parses command-line arguments to run in either screen mode or file mode.

`get_complement()`:  Constructs the complement of a DNA sequence using a dictionary to map and replace each nucleotide in the input DNA sequence with its complementary base.

`get_reverse_complement()`:  Generates the reverse complement of a given DNA sequence. This involves creating the complement of the sequence and then reversing it.

`calculate_GC_pct()`: Computes the GC content percentage of a given DNA sequence.

`calculate_AT_pct()`:  Computes the AT content percentage of a given DNA sequence.

`methyl_conversion()`: Simulates the conversion of unmethylated cytosine (C) to thymine (T) in a DNA sequence, representing the process occurring during bisulfite treatment or enzymatic methyl conversion reactions. This conversion is typically used in DNA methylation studies to identify unmethylated cytosines.

`transcribe()`:  Simulates the transcription of a DNA sequence into an RNA sequence by replacing T with U.

`translate()`:  Simulates the translation of a DNA sequence into polypeptide sequences for all 6 reading frames, including three frames from the original sequence and three from its reverse complement. This is achieved by mapping DNA codons to amino acids according to a translation table.

`print_functions()`: Prints various properties and translations of the given DNA object to stdout.

`protein_printer()`:  Prints the protein sequences contained in the provided dictionary, indicating the corresponding reading frames and code types.

`write_functions()`:  Writes various properties and analyses of a DNA sequence to a CSV file.
This function takes a DNA object and a file name as input. It then calculates various properties of the DNA sequence, including its length, complement, reverse complement, GC content percentage, AT content percentage, result of methylation conversion, and RNA transcription. It also includes the translation of the DNA sequence into polypeptide sequences for all 6 reading frames, both in 3-letter and 1-letter amino acid codes. These properties are written to the specified CSV file, with a header row defining each column if the file does not already exist.

## Usage
The script can be run in two modes:

### Screen Mode
Directly input a DNA sequence via command-line arguments:

`python basic_DNA_tools.py -s agttcacatgactgaactctg -n Sample_Sequence`

Example output to the screen:

```
Running screen mode

DNA name: Sample_Sequence
DNA sequence: AGTTCACATGACTGAACTCTG
Length: 21
Complement: TCAAGTGTACTGACTTGAGAC
Reverse complement: CAGAGTTCAGTCATGTGAACT

GC content: 42.86
AT content: 57.14

methylation conversion: AGTTtAtATGAtTGAAtTtTG
RNA transcription: AGUUCACAUGACUGAACUCUG

3-letter code translation:
Top strand, Reading frame 1 (3-letter code): SerSerHisAsp*
Top strand, Reading frame 2 (3-letter code): ValHisMetThrGluLeu
Top strand, Reading frame 3 (3-letter code): PheThr*
Bottom strand, Reading frame 4 (3-letter code): GlnSerSerValMet*
Bottom strand, Reading frame 5 (3-letter code): ArgValGlnSerCysGlu
Bottom strand, Reading frame 6 (3-letter code): GluPheSerHisValAsn

1-letter code translation:
Top strand, Reading frame 1 (1-letter code): SSHD*
Top strand, Reading frame 2 (1-letter code): VHMTEL
Top strand, Reading frame 3 (1-letter code): FT*
Bottom strand, Reading frame 4 (1-letter code): QSSVM*
Bottom strand, Reading frame 5 (1-letter code): RVQSCE
Bottom strand, Reading frame 6 (1-letter code): EFSHVN
```

### File Mode

Analyze sequences stored in FASTA or FASTQ files:

`python basic_DNA_tools.py -f ./test_seq.fa`

`python basic_DNA_tools.py -f ./test_seq.fastq`

You can also specify an output path for the results:

`python basic_DNA_tools.py -f ./test_seq.fa -o ./test_FASTA_output.csv`

`python basic_DNA_tools.py -f ./test_seq.fastq -o ./test_FASTq_output.csv`

Example output to a CSV file:
```
Name,Sequence,Quality Scores,Length,Complement,Reverse Complement,GC %,AT %,Methylation Conversion,RNA transcription,"Top strand, Reading frame 1 (3-letter code)","Top strand, Reading frame 2 (3-letter code)","Top strand, Reading frame 3 (3-letter code)","Bottom strand, Reading frame 4 (3-letter code)","Bottom strand, Reading frame 5 (3-letter code)","Bottom strand, Reading frame 6 (3-letter code)","Top strand, Reading frame 1 (1-letter code)","Top strand, Reading frame 2 (1-letter code)","Top strand, Reading frame 3 (1-letter code)","Bottom strand, Reading frame 4 (1-letter code)","Bottom strand, Reading frame 5 (1-letter code)","Bottom strand, Reading frame 6 (1-letter code)"
test sequence 1,AGTTCACATGACTGAACTCTG,,21,TCAAGTGTACTGACTTGAGAC,CAGAGTTCAGTCATGTGAACT,42.86,57.14,AGTTtAtATGAtTGAAtTtTG,AGUUCACAUGACUGAACUCUG,SerSerHisAsp*,ValHisMetThrGluLeu,PheThr*,GlnSerSerValMet*,ArgValGlnSerCysGlu,GluPheSerHisValAsn,SSHD*,VHMTEL,FT*,QSSVM*,RVQSCE,EFSHVN
test sequence 2,ATGAGTAACAATGAGTTCACATGACTGAACTCTG,,34,TACTCATTGTTACTCAAGTGTACTGACTTGAGAC,CAGAGTTCAGTCATGTGAACTCATTGTTACTCAT,38.24,61.76,ATGAGTAAtAATGAGTTtAtATGAtTGAAtTtTG,AUGAGUAACAAUGAGUUCACAUGACUGAACUCUG,MetSerAsnAsnGluPheThr*,*,Glu*,GlnSerSerValMet*,ArgValGlnSerCysGluLeuIleValThrHis,GluPheSerHisValAsnSerLeuLeuLeu,MSNNEFT*,*,E*,QSSVM*,RVQSCELLVTH,EFSHVNSLLL
test sequence 3,ATCGACCGTTCACATGACTGAACTCTG,,27,TAGCTGGCAAGTGTACTGACTTGAGAC,CAGAGTTCAGTCATGTGAACGGTCGAT,48.15,51.85,ATtGAttGTTtAtATGAtTGAAtTtTG,AUCGACCGUUCACAUGACUGAACUCUG,IleAspArgSerHisAsp*,SerThrValHisMetThrGluLeu,ArgProPheThr*,GlnSerSerValMet*,ArgValGlnSerCysGluArgSer,GluPheSerHisValAsnGlyArg,LDRSHD*,STVHMTEL,RPFT*,QSSVM*,RVQSCERS,EFSHVNGR
test sequence 4,ATGAGTAACAATGAGTTCACATGACTGAACTCTG,,34,TACTCATTGTTACTCAAGTGTACTGACTTGAGAC,CAGAGTTCAGTCATGTGAACTCATTGTTACTCAT,38.24,61.76,ATGAGTAAtAATGAGTTtAtATGAtTGAAtTtTG,AUGAGUAACAAUGAGUUCACAUGACUGAACUCUG,MetSerAsnAsnGluPheThr*,*,Glu*,GlnSerSerValMet*,ArgValGlnSerCysGluLeuIleValThrHis,GluPheSerHisValAsnSerLeuLeuLeu,MSNNEFT*,*,E*,QSSVM*,RVQSCELLVTH,EFSHVNSLLL
```

## Additional Options
-n, --name:  Optional name for the DNA sequence (only in screen mode).

-o, --output: Path to the output CSV file (mainly for file mode).
