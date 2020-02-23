#Bed file:

#chr1:117223140-117223856 3 7
#chr1:117223140-117223856 5 9

#Fasta file:
#>chr1:117223140-117223856
#CGCGTGGGCTAGGGGCTAGCCCC

#Desired output:
#>chr1:117223140-117223856
#CGTGG
#>chr1:117223140-117223856
#TGGGC


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

# read names and postions from bed file
positions = defaultdict(list)
with open('pos_forward.bed') as f:
    for line in f:
        name, start, stop = line.split()
        positions[name].append((int(start), int(stop)))

# parse faste file and turn into dictionary
records = SeqIO.to_dict(SeqIO.parse(open('hg19.ucsc.fa'), 'fasta'))

# search for short sequences
short_seq_records = []
for name in positions:
    for (start, stop) in positions[name]:
        long_seq_record = records[name]
        long_seq = long_seq_record.seq
        alphabet = long_seq.alphabet
        short_seq = str(long_seq)[start-1:stop]
        short_seq = short_seq.upper()
        #short_seq_record = SeqRecord(Seq(short_seq, alphabet), id=name, description='')
        short_seq_record = SeqRecord(Seq(short_seq, alphabet), id=name, description= str(start) + ':' + str(stop))

        short_seq_records.append(short_seq_record)

# write to file
with open('output_seq_forward.fasta', 'w') as f:
    SeqIO.write(short_seq_records, f, 'fasta')


