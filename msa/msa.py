import os, sys, re
import optparse
import subprocess

#-------------#
# Config      #
#-------------#
account = None
cluster_name = "hydra.prib.upf.edu"
memory = "512Gb"
platform = "slurm"
processors = 32
queue = None
user_email = "oriol@cmmt.ubc.ca"
walltime = "168:00:00"
# e.g. conda install -c bioconda clustalo
clustalo_path = "/home/ofornes/.anaconda2.7/bin"
# e.g. conda install -c bioconda mmseqs2
mmseqs_path = "/home/ofornes/.anaconda2.7/bin"
# i.e. dir where uniref.sh was exec
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

    parser = optparse.OptionParser("python2.7 %prog -i <input_file> [--dummy=dummy_dir -j <max_iterations> -n <nr_database> -r <redundant_db> -o <output_dir> -s <max_sequences> -v]")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="<dummy_dir>")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (in FASTA format)", metavar="<input_file>")
    parser.add_option("-j", action="store", default=8, type="int", dest="max_iterations", help="Max. number of iterations (default=8)", metavar="<max_iterations>")
    parser.add_option("-n", action="store", default="uniref50", type="string", dest="nr_db", help="Non-redundant database (\"uniref50\", \"uniref90\" or \"uniref100\"; default=uniref50)", metavar="<nr_db>")
    parser.add_option("-o", action="store", default="./", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="<output_dir>")
    parser.add_option("-r", action="store", default="uniref100", type="string", dest="redundant_db", help="Redundant database (\"uniref50\", \"uniref90\" or \"uniref100\"; default=uniref100)", metavar="<redundant_db>")
    parser.add_option("-s", action="store", default=25000, type="int", dest="max_sequences", help="Max. number of sequences (default=25000)", metavar="<max_sequences>")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")

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
    msa = {}
    nr_db = os.path.join(uniref_path, options.nr_db)
    redundant_db = os.path.join(uniref_path, options.redundant_db)
    dummy_dir = os.path.abspath(dummy_dir)
    
    # Create output dir #
    if not os.path.exists(os.path.abspath(options.output_dir)):
            os.makedirs(os.path.abspath(options.output_dir))
    
    # Skip if query file already exists #
    query_file = os.path.join(os.path.abspath(options.output_dir), "query.fa")
    if not os.path.exists(query_file):
        # For header, sequence... #
        for header, sequence in parse_fasta_file(os.path.abspath(options.input_file)):
            # Write #
            write(query_file, ">%s\n%s" % (header, sequence))
            # Add to MSA #
            msa.setdefault(sequence, header)
    
    # Skip if query db already exists #
    query_db = os.path.join(os.path.abspath(options.output_dir), "query.%s.it_1.db" % options.nr_db)
    if not os.path.exists(query_db):
        # Create DB #
        process = subprocess.check_output([os.path.join(mmseqs_path, "mmseqs"), "createdb", query_file, query_db])
    
    # For each iteration... #
    for i in range(8):
        # Skip if alignment file already exists #
        alignment_file = os.path.join(os.path.abspath(options.output_dir), "query.%s.it_%s.ali" % (options.nr_db, i + 1))
        if not os.path.exists(alignment_file):
            # Search DB #
            process = subprocess.check_output([os.path.join(mmseqs_path, "mmseqs"), "search", query_db, nr_db, alignment_file, dummy_dir, "--max-seqs", "300", "--split-memory-limit", "512000000000", "--threads", "32", "-s", "7.5"])
        # If alignment has enough sequences or this is the last iteration... #
        if (len([j for j in parse_file(alignment_file)]) == 300) or (i + 1 == 8):
            next_query_db = os.path.join(os.path.abspath(options.output_dir), "query.%s.db" % options.redundant_db)
        # ... Else... #
        else:
            next_query_db = os.path.join(os.path.abspath(options.output_dir), "query.%s.it_%s.db" % (options.nr_db, i + 2))
        # Skip if next query DB already exists #
        if not os.path.exists(next_query_db):
            # Create DB #
            process = subprocess.check_output([os.path.join(mmseqs_path, "mmseqs"), "result2profile", query_db, nr_db, alignment_file, next_query_db])
            query_db = next_query_db
        # End for loop if switched to redundant db #
        if options.redundant_db in query_db: break

    # Skip if alignment file already exists #
    alignment_file = os.path.join(os.path.abspath(options.output_dir), "query.%s.ali" % options.redundant_db)
    if not os.path.exists(alignment_file):
        # Search DB #
        process = subprocess.check_output([os.path.join(mmseqs_path, "mmseqs"), "search", query_db, redundant_db, alignment_file, dummy_dir, "--max-seqs", str(options.max_sequences * 2), "--split-memory-limit", "512000000000", "--threads", "32", "-s", "7.5", "--max-seq-id", "0.999"])
    
    # Skip if sequences file already exists #
    sequences_file = os.path.join(os.path.abspath(options.output_dir), "query.%s.fa" % options.redundant_db)
    if not os.path.exists(sequences_file):
        # Get FASTA sequences #
        process = subprocess.check_output([os.path.join(mmseqs_path, "mmseqs"), "createseqfiledb", redundant_db, alignment_file, sequences_file])

        #cat ./examples/1ATG_A/query.fa ./examples/1ATG_A/query_uniref100.ali.fa > ./examples/1ATG_A/query_clustalo_in.fa
#clustalo -i ./examples/1ATG_A/query_clustalo_in.fa -o ./examples/1ATG_A/query_clustalo_out.fa
#
#    [ofornes@cdr462 RADI]$ clustalo -i ./examples/1ATG_A/query_clustalo_in.fa -o ./examples/1ATG_A/query_clustalo_out.fa --threads=32 -v
#    Using 32 threads
#    Read 10001 sequences (type: Protein) from ./examples/1ATG_A/query_clustalo_in.fa
#    Using 176 seeds (chosen with constant stride from length sorted seqs) for mBed (from a total of 10001 sequences)
#    Calculating pairwise ktuple-distances...
#    Ktuple-distance calculation progress done. CPU time: 144.20u 0.09s 00:02:24.29 Elapsed: 00:02:24
#    mBed created 181 cluster/s (with a minimum of 1 and a soft maximum of 100 sequences each)
#    Distance calculation within sub-clusters done. CPU time: 40.49u 0.01s 00:00:40.50 Elapsed: 00:00:40
#    Guide-tree computation (mBed) done.
#    Progressive alignment progress: 99 % (9900 out of 10000)
#    Progressive alignment progress done. CPU time: 390.10u 161.68s 00:09:11.78 Elapsed: 00:09:12
#    Alignment written to ./examples/1ATG_A/query_clustalo_out.fa