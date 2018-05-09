import os, sys, re
import optparse
import subprocess

# i.e. directory where uniref.sh was exec
uniref_path = "/home/ofornes/scratch/RADI/uniref"

#-------------#
# Parsers     #
#-------------#

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

#-------------#
# Write       #
#-------------#

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
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python2.7 %prog -i <input_file> [--dummy=dummy_dir -n <nr_database> -r <db> -o <output_dir> -s <max_sequences> -u <uniref_dir>]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="<dummy_dir>")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in FASTA format)", metavar="<input_file>")
    parser.add_option("-n", action="store", default="uniref50", type="string", dest="nr_db", help="Non-redundant database (\"uniref50\", \"uniref90\" or \"uniref100\"; default=uniref50)", metavar="<nr_db>")
    parser.add_option("-o", action="store", default="./", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="<output_dir>")
    parser.add_option("-r", action="store", default="uniref100", type="string", dest="db", help="Redundant database (\"uniref50\", \"uniref90\" or \"uniref100\"; default=uniref100)", metavar="<db>")
    parser.add_option("-s", action="store", default=1000000, type="int", dest="max_sequences", help="Max. number of sequences (default=1000000)", metavar="<max_sequences>")
    parser.add_option("-u", action="store", type="string", dest="uniref_dir", help="UniRef directory (i.e. where uniref.sh was exec)", metavar="<uniref_dir>")

    (options, args) = parser.parse_args()

    if options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Initialize #
    nr_db = os.path.join(uniref_path, "%s.db" % options.nr_db)
    db = os.path.join(uniref_path, "%s.db" % options.db)
    dummy_dir = os.path.abspath(options.dummy_dir)
    # Create output dir #
    if not os.path.exists(os.path.abspath(options.output_dir)):
        os.makedirs(os.path.abspath(options.output_dir))
    # Create dummy dir #
    if not os.path.exists(dummy_dir):
        os.makedirs(dummy_dir)

    #----------#
    # MMseqs2  #
    #----------#
        
    # Skip if redundant query db already exists #
    query_db = os.path.join(os.path.abspath(options.output_dir), "query.%s.db" % options.db)
    if not os.path.exists(query_db):
        # Skip if nr query db already exists #
        nr_query_db = os.path.join(os.path.abspath(options.output_dir), "query.%s.db" % options.nr_db)
        if not os.path.exists(query_db):
            # Create DB #
            process = subprocess.check_output(["mmseqs", "createdb", os.path.abspath(options.input_file), nr_query_db])
        # Skip if alignment file already exists #
        alignment_file = os.path.join(os.path.abspath(options.output_dir), "query.%s.ali" % options.nr_db)
        if not os.path.exists(alignment_file):
            # Search DB #
            process = subprocess.check_output(["mmseqs", "search", nr_query_db, nr_db, alignment_file, dummy_dir, "--threads", "32", "-s", "7.5", "--num-iterations", "4"])
        # Create DB #
        process = subprocess.check_output(["mmseqs", "result2profile", nr_query_db, nr_db, alignment_file, query_db])
    # Skip if alignment file already exists #
    alignment_file = os.path.join(os.path.abspath(options.output_dir), "query.%s.ali" % options.db)
    if not os.path.exists(alignment_file):
        # Search DB #
        process = subprocess.check_output(["mmseqs", "search", query_db, db, alignment_file, dummy_dir, "--max-seqs", str(options.max_sequences), "--threads", "32", "-s", "7.5", "--max-seq-id", "1.0"])
    # Skip if nr sequences file already exists #
    nr_sequences_file = os.path.join(os.path.abspath(options.output_dir), "query.%s.fa" % options.nr_db)
    if not os.path.exists(nr_sequences_file):
        # Get FASTA sequences #
        process = subprocess.check_output(["mmseqs", "createseqfiledb", nr_db, alignment_file, nr_sequences_file])
    # Skip if redundant sequences file already exists #
    sequences_file = os.path.join(os.path.abspath(options.output_dir), "query.%s.fa" % options.db)
    if not os.path.exists(sequences_file):
        # Get FASTA sequences #
        process = subprocess.check_output(["mmseqs", "createseqfiledb", db, alignment_file, sequences_file])

    #----------#
    # ClustalO #
    #----------#

    # Skip if ClustalO input file already exists #
    clustalo_in_file = os.path.join(os.path.abspath(options.output_dir), "clustalo.in.fa")
    if not os.path.exists(clustalo_in_file):
        # Initialize #
        msa = []
        sequences = {}
        # For header, sequence... #
        for header, sequence in parse_fasta_file(os.path.abspath(options.input_file)):
            # Add to MSA #
            msa.append((header, sequence))
            sequences.setdefault(sequence, header)
            break
        # For header, sequence... #
        for header, sequence in parse_fasta_file(nr_sequences_file):
            # Skip if sequence already exists #
            if sequence in sequences: continue
            # Add to MSA #
            msa.append((header, sequence))
            sequences.setdefault(sequence, header)
        # For header, sequence... #
        for header, sequence in msa:
            # Write #
            write(clustalo_in_file, ">%s\n%s" % (header, sequence))
    # Skip if ClustalO output file already exists #
    clustalo_out_file = os.path.join(os.path.abspath(options.output_dir), "clustalo.out.fa")
    if not os.path.exists(clustalo_out_file):
        # Create MSA #
        process = subprocess.check_output(["clustalo", "-i", clustalo_in_file, "-o", clustalo_out_file,"--threads", "32"])

    #----------#
    # HMMER    #
    #----------#

    # Skip if HMMalign input file already exists #
    hmmalign_in_file = os.path.join(os.path.abspath(options.output_dir), "hmmalign.in.fa")
    if not os.path.exists(hmmer_in_file):
        # Initialize #
        msa = []
        sequences = {}
        # For header, sequence... #
        for header, sequence in parse_fasta_file(sequences_file):
            # Skip if sequence already exists #
            if sequence in sequences: continue
            # Add to MSA #
            msa.append((header, sequence))
            sequences.setdefault(sequence, header)
        # For header, sequence... #
        for header, sequence in msa:
            # Write #
            write(hmmalign_in_file, ">%s\n%s" % (header, sequence))
    # Skip if HMMalign input hmm already exists #
    hmmalign_in_hmm = os.path.join(os.path.abspath(options.output_dir), "hmmalign.in.hmm")
    if not os.path.exists(hmmalign_in_hmm):
        # Create HMM #
        process = subprocess.check_output(["hmmbuild", hmmalign_in_hmm, clustalo_out_file]) 
    # Skip if HMMER output file already exists #
    hmmalign_out_hmm = os.path.join(os.path.abspath(options.output_dir), "hmmalign.out.a2m")
    if not os.path.exists(hmmalign_out_hmm):
        # Create MSA #
        process = subprocess.check_output(["hmmalign", "--mapali", clustalo_out_file, "-o", hmmalign_out_hmm, "--outformat", "A2M", "--trim", hmmalign_in_hmm, hmmalign_in_file])

#    msa_file = os.path.join(os.path.abspath(options.output_dir), "msa.fa")
#    if not os.path.exists(msa_file):
#        # Initialize #
#        headers = []
#        sequences = []
#        # For header, sequence... #
#        for header, sequence in parse_fasta_file(famsa_out_file, clean=False):
#            # Add to lists #
#            headers.append(header)
#            sequences.append(list(sequence))
#        # Transpose sequences #
#        sequences = zip(*sequences)
#        # For each position... #
#        for i in reversed(range(len(sequences))):
#            # If query position is a gap... #
#            if sequences[i][0] == "-": sequences.pop(i)
#        # Transpose sequences #
#        sequences = zip(*sequences)
#        # For each sequence... #
#        for i in range(len(headers)):
#            # Write #
#            write(msa_file, ">%s\n%s" % (headers[i], "".join(sequences[i])))
