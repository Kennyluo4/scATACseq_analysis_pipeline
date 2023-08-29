import pysam

file = 'wgEncodeUwRepliSeqBjS3AlnRep1.bam'

samf = pysam.AlignmentFile(file, 'rb')

for reads in samf:
    print(reads)
    
samf.check_index()