import basic_DNA_tools as tools
import pytest


def test_valid_DNA_object():
    """
    Test the creation of a DNA object with a valid sequence.
    Asserts that the sequence is correctly stored and uppercase.
    """
    seq = "agttcacatgactgaactctg"
    name = "test_seq"
    dna = tools.DNA(seq, name)
    assert dna.seq == "AGTTCACATGACTGAACTCTG"
    assert dna.name == "test_seq"


def test_invalid_DNA_object():
    seq = "agttcacatgnactgaactctg"
    name = "test_seq"
    with pytest.raises(ValueError):
        tools.DNA(seq, name) == ValueError


def test_run_screen_mode_without_output_path():
    seq = "agttcacatgactgaactctg"
    name = "test_seq"
    # Assume function works if no exceptions are raised
    try:
        tools.run_screen_mode(seq, name)
    except Exception as e:
        pytest.fail(f"run_screen_mode without output path raised an exception: {e}")


def test_run_screen_mode_with_output_file():
    seq = "agttcacatgactgaactctg"
    name = "test_seq"
    output = "test_screen_mode.csv"
    try:
        tools.run_screen_mode(seq, name, output)
    except Exception as e:
        pytest.fail(f"run_screen_mode with output path raised an exception: {e}")


def test_run_valid_FASTA_file_mode():
    file = "test_seq.fa"
    output = "test_valid_FASTA_file.csv"
    try:
        tools.run_file_mode(file, output)
    except Exception as e:
        pytest.fail(f"run_file_mode with valid FASTA file raised an exception: {e}")


def test_run_invalid_FASTA_extension():
    file = "test_seq.seq"
    output = "test_invalid_FASTA_extension.csv"
    with pytest.raises(SystemExit):
        tools.run_file_mode(file, output)


def test_run_invalid_FASTA_format():
    file = "invalid_seq.fa"
    output = "test_invalid_FASTA_format.csv"
    with pytest.raises(SystemExit):
        tools.run_file_mode(file, output)


def test_run_valid_FASTQ_file_mode():
    file = "test_seq.fastq"
    output = "test_valid_FASTQ_file.csv"
    try:
        tools.run_file_mode(file, output)
    except Exception as e:
        pytest.fail(f"run_file_mode with valid FASTQ file raised an exception: {e}")
    

def test_run_invalid_FASTQ_format():
    file = "invalid_seq.fastq"
    output = "test_invalid_FASTQ_format.csv"
    with pytest.raises(SystemExit):
        tools.run_file_mode(file, output)


def test_get_complement():
    seq = "agttcacatgactgaactctg"
    name = "test_seq"
    dna = tools.DNA(seq, name)
    assert tools.get_complement(dna) == "TCAAGTGTACTGACTTGAGAC"


def test_get_reverse_complement():
    seq = "agttcacatgactgaactctg"
    name = "test_seq"
    dna = tools.DNA(seq, name)
    assert tools.get_reverse_complement(dna) == "CAGAGTTCAGTCATGTGAACT"


def test_calculate_GC_pct():
    seq = "gcatcatagcttagtatcga"
    name = "test_seq"
    dna = tools.DNA(seq, name)
    assert tools.calculate_GC_pct(dna) == 40.00


def test_calculate_AT_pct():
    seq = "gcatcatagcttagtatcga"
    name = "test_seq"
    dna = tools.DNA(seq, name)
    assert tools.calculate_AT_pct(dna) == 60.00


def test_methyl_conversion():
    seq = "gcatcatagcttagtatcga"
    name = "test_seq"
    dna = tools.DNA(seq, name)
    assert tools.methyl_conversion(dna) == "GtATtATAGtTTAGTATtGA"


def test_transcribe():
    seq = "gcatcatagcttagtatcga"
    name = "test_seq"
    dna = tools.DNA(seq, name)
    assert tools.transcribe(dna) == "GCAUCAUAGCUUAGUAUCGA"


def test_translate_3_letter_code():
    seq = "AGTTCACATGACTGAACTCTG"
    name = "test_seq"
    dna = tools.DNA(seq, name)
    proteins = {
        "Top strand, Reading frame 1 (3-letter code)": "SerSerHisAsp*",
        "Top strand, Reading frame 2 (3-letter code)": "ValHisMetThrGluLeu",
        "Top strand, Reading frame 3 (3-letter code)": "PheThr*",
        "Bottom strand, Reading frame 4 (3-letter code)": "GlnSerSerValMet*",
        "Bottom strand, Reading frame 5 (3-letter code)": "ArgValGlnSerCysGlu",
        "Bottom strand, Reading frame 6 (3-letter code)": "GluPheSerHisValAsn",      
    }
    assert tools.translate(dna, letters=3) == proteins


def test_translate_1_letter_code():
    seq = "AGTTCACATGACTGAACTCTG"
    name = "test_seq"
    dna = tools.DNA(seq, name)
    proteins = {
        "Top strand, Reading frame 1 (1-letter code)": "SSHD*",
        "Top strand, Reading frame 2 (1-letter code)": "VHMTEL",
        "Top strand, Reading frame 3 (1-letter code)": "FT*",
        "Bottom strand, Reading frame 4 (1-letter code)": "QSSVM*",
        "Bottom strand, Reading frame 5 (1-letter code)": "RVQSCE",
        "Bottom strand, Reading frame 6 (1-letter code)": "EFSHVN",
    }
    assert tools.translate(dna) == proteins


def test_print_functions():
    seq = "AGTTCACATGACTGAACTCTG"
    name = "test_seq"
    dna = tools.DNA(seq, name)
    try:
        tools.print_functions(dna)
    except Exception as e:
        pytest.fail(f"print functions raised an exception {e}")


def test_protein_printer():
    proteins = {
        "Top strand, Reading frame 1 (1-letter code)": "SSHD*",
        "Top strand, Reading frame 2 (1-letter code)": "VHMTEL",
        "Top strand, Reading frame 3 (1-letter code)": "FT*",
        "Bottom strand, Reading frame 4 (1-letter code)": "QSSVM*",
        "Bottom strand, Reading frame 5 (1-letter code)": "RVQSCE",
        "Bottom strand, Reading frame 6 (1-letter code)": "EFSHVN",
    }
    try:
        tools.protein_printer(proteins)
    except Exception as e:
        pytest.fail(f"protein printer raised an exception {e}")


def test_write_functions():
    seq = "AGTTCACATGACTGAACTCTG"
    name = "test_seq"
    dna = tools.DNA(seq, name)
    try:
        tools.write_functions(dna, file_name="test_write_functions.csv")
    except Exception as e:
        pytest.fail(f"print functions raised an exception {e}")