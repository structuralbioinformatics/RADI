#!/usr/bin/env python2.7
import os, sys, re
import optparse
import subprocess

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("./%prog -f <famsa_dir> -i <input_file> -m <mmseqs_dir> -u <uniref_dir> [-n <nr_db> -r <db>] [--dummy=<dummy_dir> -o <output_dir> -s <max_sequences> -t <threads>")

    parser.add_option("-f", action="store", type="string", dest="famsa_dir", help="Full path to FAMSA bin directory (i.e. where \"famsa\" is located; e.g. $FAMSA_PATH/bin)", metavar="<famsa_dir>")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in FASTA format)", metavar="<input_file>")
    parser.add_option("-m", action="store", type="string", dest="mmseqs_dir", help="Full path to MMseqs2 bin directory (i.e. where \"mmseqs\" is located; e.g. $MMSEQS2_PATH/bin)", metavar="<mmseqs_dir>")
    parser.add_option("-u", action="store", type="string", dest="uniref_dir", help="UniRef DB directory (i.e. where uniref.sh was executed)", metavar="<uniref_dir>")

    group = optparse.OptionGroup(parser, "UniRef DB options", "The script will build a profile of the query by searching the non-redundant DB (default = \"uniref50\") with MMseqs2 for 4 iterations, after which it will switch to searching the redundant DB (default = \"uniref100\") with that profile.")
    group.add_option("-n", action="store", default="uniref50", type="string", dest="nr_db", help="Non-redundant UniRef DB (\"uniref50\", \"uniref90\" or \"uniref100\"; default = uniref50)", metavar="<nr_db>")
    group.add_option("-r", action="store", default="uniref100", type="string", dest="db", help="Redundant UniRef DB (\"uniref50\", \"uniref90\" or \"uniref100\"; default = uniref100)", metavar="<db>")
    parser.add_option_group(group)

    group = optparse.OptionGroup(parser, "Non-mandatory options")
    group.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="<dummy_dir>")
    group.add_option("-o", action="store", default="./", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="<output_dir>")
    group.add_option("-s", action="store", default=100000, type="int", dest="max_sequences", help="Max. number of sequences (default = 100000)", metavar="<max_sequences>")
    group.add_option("-t", "--threads", default=1, action="store", type="int", dest="threads", help="Total number of cores to be used for the computation (default = 1)", metavar="<threads>")
    parser.add_option_group(group)

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def parse_file(file_name):
    """
    This function parses any file and yields lines one by one.
    
    @input:
    file_name {string}
    @return:
    line {string}

    """
 
    if os.path.exists(file_name):
        # Initialize #
        f = None
        # Open file handle #
        try: f = open(file_name, "rt")
        except: raise ValueError("Could not open file %s" % file_name)
        # For each line... #
        for line in f:
            yield line.strip("\n")
        f.close()
    else:
        raise ValueError("File %s does not exist!" % file_name)

def parse_fasta_file(file_name, clean=True):
    """
    This function parses any FASTA file and yields sequences one by one
    in the form header, sequence.

    @input:
    file_name {string}
    clean {boolean} if true, converts \W to X
    @return:
    line {list} header, sequence

    """

    # Initialize #
    header = ""
    sequence = ""
    # For each line... #
    for line in parse_file(file_name):
        if len(line) == 0: continue
        if line.startswith("#"): continue
        if line.startswith(">"):
            if sequence != "":
                if clean:
                    sequence = re.sub("\W|\d", "X", sequence)
                yield header, sequence
            m = re.search("^>(.+)", line)
            header = m.group(1)
            sequence = ""
        else:
            sequence += line.upper()
    if clean:
        sequence = re.sub("\W|\d", "X", sequence)

    yield header, sequence

def write(file_name=None, content=""):
    """
    This function writes any {content} to a file or to stdout if no
    file is provided. If the file already exists, it pushed the {content}
    at the bottom of the file.

    @input:
    file_name {string}
    content {string}

    """
    if file_name is not None:
        try:
            f = open(file_name, "a")
            f.write("%s\n" % content)
        except:
            raise ValueError("Could create file %s" % file_name)
    else:
        sys.stdout.write("%s\n" % content)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Initialize #
    nr_db = os.path.join(os.path.abspath(options.uniref_dir), "%s.fa.db" % options.nr_db)
    db = os.path.join(os.path.abspath(options.uniref_dir), "%s.fa.db" % options.db)
    # Create dummy dir #
    dummy_dir = os.path.join(os.path.abspath(options.dummy_dir), "%s.%s" % (os.path.basename(__file__), os.getpid()))
    if not os.path.exists(dummy_dir):
        os.makedirs(dummy_dir)
    # Create output dir #
    if not os.path.exists(os.path.abspath(options.output_dir)):
        os.makedirs(os.path.abspath(options.output_dir))

    #----------#
    # MMseqs2  #
    #----------#

    # Skip if nr query db already exists #
    nr_query_db = os.path.join(os.path.abspath(options.output_dir), "query.%s.db" % options.nr_db)
    if not os.path.exists(nr_query_db):
        # Create DB #
        process = subprocess.check_output(["mmseqs", "createdb", os.path.abspath(options.input_file), nr_query_db])
    # Skip if nr alignment file already exists #
    nr_alignment_file = os.path.join(os.path.abspath(options.output_dir), "query.%s.ali" % options.nr_db)
    if not os.path.exists(nr_alignment_file):
        # Search DB #
        process = subprocess.check_output(["mmseqs", "search", nr_query_db, nr_db, nr_alignment_file, dummy_dir, "--threads", options.threads, "-s", "7.5", "--max-seq-id", "1.0", "--num-iterations", "4"])
    # Skip if redundant query db already exists #
    query_db = os.path.join(os.path.abspath(options.output_dir), "query.%s.db" % options.db)
    if not os.path.exists(query_db):
        # Create DB #
        process = subprocess.check_output(["mmseqs", "result2profile", nr_query_db, nr_db, nr_alignment_file, query_db])
    # Skip if alignment file already exists #
    alignment_file = os.path.join(os.path.abspath(options.output_dir), "query.%s.ali" % options.db)
    if not os.path.exists(alignment_file):
        # Search DB #
        process = subprocess.check_output(["mmseqs", "search", query_db, db, alignment_file, dummy_dir, "--max-seqs", str(options.max_sequences), "--threads", options.threads, "-s", "7.5", "--max-seq-id", "1.0"])
    # Skip if nr sequences file already exists #
    nr_sequences_file = os.path.join(os.path.abspath(options.output_dir), "query.%s.fa" % options.nr_db)
    if not os.path.exists(nr_sequences_file):
        # Get FASTA sequences #
        process = subprocess.check_output(["mmseqs", "createseqfiledb", nr_db, nr_alignment_file, nr_sequences_file])
    # Skip if redundant sequences file already exists #
    sequences_file = os.path.join(os.path.abspath(options.output_dir), "query.%s.fa" % options.db)
    if not os.path.exists(sequences_file):
        # Get FASTA sequences #
        process = subprocess.check_output(["mmseqs", "createseqfiledb", db, alignment_file, sequences_file])

    #----------#
    # FAMSA    #
    #----------#

    # Skip if FAMSA input file already exists #
    famsa_in_file = os.path.join(os.path.abspath(options.output_dir), "famsa.in.fa")
    if not os.path.exists(famsa_in_file):
        # Initialize #
        sequences = {}
        # For header, sequence... #
        for header, sequence in parse_fasta_file(os.path.abspath(options.input_file)):
            # Add sequence #
            sequences.setdefault(sequence, header)
        # For header, sequence... #
        for header, sequence in parse_fasta_file(nr_sequences_file):
            # Skip if enough sequences #
            if len(sequences) == options.max_sequences: break
            # Skip if sequence already exists #
            if sequence in sequences: continue
            sequences.setdefault(sequence, header)
        # For header, sequence... #
        for header, sequence in parse_fasta_file(sequences_file):
            # Skip if enough sequences #
            if len(sequences) == options.max_sequences: break
            # Skip if sequence already exists #
            if sequence in sequences: continue
            sequences.setdefault(sequence, header)
        # For each sequence... #
        for sequence in sequences:
            # Write #
            write(famsa_in_file, ">%s\n%s" % (sequences[sequence], sequence))
    # Skip if FAMSA output file already exists #
    famsa_out_file = os.path.join(os.path.abspath(options.output_dir), "famsa.out.fa")
    if not os.path.exists(famsa_out_file):
        # Create MSA #
        process = subprocess.check_output(["famsa", "-t", "32", famsa_in_file, famsa_out_file])

    #----------#
    # MSAs     #
    #----------#

    # Skip if MSA file already exists #
    famsa_msa_file = os.path.join(os.path.abspath(options.output_dir), "msa.famsa.fa")
    if not os.path.exists(famsa_msa_file):
        # Initialize #
        uniq = set()
        headers = []
        sequences = []
        # For header, sequence... #
        for header, sequence in parse_fasta_file(famsa_out_file, clean=False):
            # Add to lists #
            headers.append(header)
            sequences.append(list(sequence))
        # Transpose sequences #
        sequences = zip(*sequences)
        # For each position... #
        for i in reversed(range(len(sequences))):
            # If query position is a gap... #
            if sequences[i][0] == "-": sequences.pop(i)
        # Transpose sequences #
        sequences = zip(*sequences)
        # For each sequence... #
        for i in range(len(headers)):
            # If sequence is unique... #
            sequence = "".join(sequences[i])
            if sequence not in uniq:
                # Write #
                write(famsa_msa_file, ">%s\n%s" % (headers[i], sequence))
                # Sequence is unique #
                uniq.add(sequence)
