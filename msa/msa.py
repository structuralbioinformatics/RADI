import os, sys, re
import ConfigParser
import optparse
import subprocess


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
        if gz:
            try: f = gzip.open(file_name, "rt")
            except: raise ValueError("Could not open file %s" % file_name)
        else:
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
    for line in parse_file(file_name, gz):
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

def write(file_name=None, content=None):
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
    parser.add_option("-o", action="store", default="./", type="string", dest="output_file", help="Output directory (default = ./)", metavar="<output_dir>")
    parser.add_option("-r", action="store", default="uniref100", type="string", dest="redundant_db", help="Redundant database (\"uniref50\", \"uniref90\" or \"uniref100\"; default=uniref100)", metavar="<redundant_db>")
    parser.add_option("-s", action="store", default=10000, type="int", dest="max_sequences", help="Max. number of sequences (default=10000)", metavar="<max_sequences>")
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

    # Create output dir #
    if not os.path.exists(os.path.abspath(options.output_dir)):
            os.makedirs(os.path.abspath(options.output_dir))
    
    # Initialize #
    query_file = os.path.join(os.path.abspath(options.output_dir), "query.fa")
    query_db = os.path.join(os.path.abspath(options.output_dir), "query.it_0.db")
    nr_db
    # For header, sequence... #
    for header, sequence in parse_fasta_file(os.path.abspath(options.input_file)):
        # Write #
        write(query_file, ">%s\n%s" % (header, sequence))
    # Create DB #
    process = subprocess.check_output(["mmseqs", "createdb", query_file, query_db])
    
    
#mmseqs search ./examples/1ATG_A/query.db ./uniref/uniref50 ./examples/1ATG_A/query_uniref50.ali /home/ofornes/scratch/tmp/ --split-memory-limit 512000000000 --threads 32 -s 7.5
#mmseqs result2profile ./examples/1ATG_A/query.db ./uniref/uniref50 ./examples/1ATG_A/query_uniref50.ali ./examples/1ATG_A/query_uniref50.db
#mmseqs search ./examples/1ATG_A/query_uniref50.db ./uniref/uniref100 ./examples/1ATG_A/query_uniref100.ali /home/ofornes/scratch/tmp/ --max-seqs 10000 --split-memory-limit 512000000000 --threads 32 -s 7.5 --max-seq-id 0.999
#mmseqs createseqfiledb ./uniref/uniref100 ./examples/1ATG_A/query_uniref100.ali ./examples/1ATG_A/query_uniref100.ali.fa
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