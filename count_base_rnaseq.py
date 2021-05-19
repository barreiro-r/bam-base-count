# coding=utf8
import pysam
from optparse import OptionParser

def get_pileupred_value(pileupread):
	'''(pileup) -> str'''
    if pileupread.is_del:
        value = "is_del"
    elif pileupread.is_refskip:
        value = "is_refskip"
    else:
        value = pileupread.alignment.query_sequence[pileupread.query_position]
    return value

def cont_all_reads(pileupcolumn):
	'''(pileup) -> dict'''
    values = {
        "A":0,
        "T":0,
        "C":0,
        "G":0,
        "N":0,
        "is_del":0,
        "is_refskip":0
    }

    for pileupread in pileupcolumn.pileups:
        values[get_pileupred_value(pileupread)] += 1

    return values

def print_table(positions_bases, my_chr):
	'''(int,str) -> None'''
    head = ["chr","pos","A","T","C","G","N","is_del","is_refskip"]
    print("\t".join(head))
    for my_position in sorted(list(positions_bases.keys())):
        #me deixa
        row_list = [
            my_chr,
            my_position,
            positions_bases[my_position]["A"],
            positions_bases[my_position]["T"],
            positions_bases[my_position]["C"],
            positions_bases[my_position]["G"],
            positions_bases[my_position]["N"],
            positions_bases[my_position]["is_del"],
            positions_bases[my_position]["is_refskip"]]
        print("\t".join([ str(x) for x in row_list]))

def main():
		#setting parse
    usage  = "usage: %prog [options] position1 position2 ..."
    parser = OptionParser(usage = usage)
    parser.add_option("-f", "--file", dest="bamfile",
                      help=".bam file for counting", metavar="FILE")
    parser.add_option("-c", "--chr", dest="chr",
                      help="chromosome e.g. chr1", type = 'string')

    (options, args) = parser.parse_args()

    #load bam and set other parameters
    samfile   = pysam.AlignmentFile(options.bamfile, "rb")
    my_chr    = options.chr
    positions = [ int(x) for x in args ]

    positions_bases = {}

    #populate position_bases
    for pileupcolumn in samfile.pileup(reference = my_chr, start=min(positions), stop=max(positions)+1, min_base_quality = 0):
        if pileupcolumn.reference_pos in positions and pileupcolumn.reference_name == my_chr:
            positions_bases[pileupcolumn.reference_pos] = cont_all_reads(pileupcolumn)

    #print output
    print_table(positions_bases,my_chr)

if __name__ == "__main__":
    main()

