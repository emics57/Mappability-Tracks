import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import matplotlib.patches as mplpatches
from operator import itemgetter
import argparse
import pysam

parser = argparse.ArgumentParser()
parser.add_argument('--bam','-b',default="None",type=str,action='store',help='Path to BAM file goes here')
parser.add_argument('--out','-o',default="None",type=str,action='store',help='Path to BED output goes here')

bamInput=parser.parse_args().bam

# convert BAM file into a Dataframe with read ID, SAMflag, mapped chr name, CIGAR string,
# alignment score, NM score (edit distance)
save = pysam.set_verbosity(0)
bam_data = []
alignments = pysam.AlignmentFile(bamInput,"rb")
pysam.set_verbosity(save)

for read in alignments:
    readName = read.query_name # read id
    derivedName = readName.split(':')[0] # chromosome the read was derived from
    flag = read.flag # SAM flag
    cigarSize = read.infer_query_length() # read size inferred from cigar string
    rname = alignments.get_reference_name(read.reference_id) # mapped chromosome name
    ref_start = read.reference_start # mapped coordinate
    as_score = read.get_tag("AS") if read.has_tag("AS") else None
    nm_score = read.get_tag("NM") if read.has_tag("NM") else None
    bam_data.append([readName,derivedName,flag,rname,cigarSize,ref_start,as_score,nm_score])
df = pd.DataFrame(bam_data)
chr = df.loc[0,1]
readSize=df.loc[0,4]
output=parser.parse_args().out

# obtain uniquely mapped reads
value_counts = df[0].value_counts()
unique_mapped_values = value_counts[value_counts == 1].index.tolist()
unique_mapped_df = df[df[0].isin(unique_mapped_values)& (df[1] != '4')]
unique_mapped_df[5] = unique_mapped_df[5].astype(int)
# mapped coordinates of uniquely mapped reads
uniqueList = list(zip(unique_mapped_df[5], unique_mapped_df[4]))
uniqueList = sorted(uniqueList, key=lambda x: x[1])

# plots just the multimapped reads with three or less top alignments
# Extract values that appear more than once
multi_mapped_values = value_counts[value_counts > 1].index.tolist()
# Filter the DataFrame to get the multi-mapped entries
multi_mapped_df = df[df[0].isin(multi_mapped_values)& (df[1] != '4')]
# filtering for just the TOP multi-mapped reads
topAlignmentsDF = multi_mapped_df[(multi_mapped_df[6] == 2*readSize) & ( multi_mapped_df[7] == 0)]
topReadCounts = topAlignmentsDF[0].value_counts()
# Get a list of alignments that appear 2 times or less
reads_to_keep = topReadCounts[topReadCounts <= 2].index
# Use the isin function to filter the DataFrame for top alignments
top2AlignmentsDF= topAlignmentsDF[topAlignmentsDF[0].isin(reads_to_keep)]
top2AlignmentsDF[5] = top2AlignmentsDF[5].astype(int)
# filteredList = top2AlignmentsDF[3].values.tolist()
# mapped coordinates of uniquely mapped reads
top2AlignsList = list(zip(top2AlignmentsDF[5], top2AlignmentsDF[4]))
top2AlignsList = sorted(top2AlignsList, key=lambda x: x[1])
# MULTIMAPPED READS
# Extract read IDs that appear more than once
# multi_mapped_values = value_counts[value_counts > 1].index.tolist()
# # Filter all reads to get the multi-mapped reads
# multi_mapped_df = df[df[0].isin(multi_mapped_values)]
# # is this necessary (yes it is)
# # reads that map to the correct chromosome
# sameChrDF = multi_mapped_df[multi_mapped_df[2].str.startswith(chr[:5])]
# # reads that map to a different chromosome
# diffChrDF = multi_mapped_df[~multi_mapped_df[2].str.startswith(chr[:5])] 

# # exclude reads that are in the diffChr dataframe
# sameChrDF = sameChrDF[~sameChrDF[0].isin(diffChrDF[0])]
# # filtering for just the top multi aligned reads
# completeAlignment = sameChrDF[(sameChrDF[3] == f'{readSize}M') & ( sameChrDF[4] == 2*readSize) & (sameChrDF[5] == 0)]

# read_counts = completeAlignment[0].value_counts()

# # Get a list of reads that appear 3 times or less
# readsTopAligns = read_counts[read_counts <= 2].index
# reads3Aligns = read_counts[read_counts == 3].index
# reads2Aligns = read_counts[read_counts == 2].index
# reads1Aligns = read_counts[read_counts == 1].index

# # Use the isin function to filter the DataFrame for top alignments
# readsTopAlignsDF= completeAlignment[completeAlignment[0].isin(readsTopAligns)]
# reads3AlignDF = completeAlignment[completeAlignment[0].isin(reads3Aligns)]
# reads2AlignDF = completeAlignment[completeAlignment[0].isin(reads2Aligns)]
# reads1AlignDF = completeAlignment[completeAlignment[0].isin(reads1Aligns)]

# # filter for just top 2 alignments to same hap
# newDF = reads2AlignDF[reads2AlignDF[2]==chr]
# newDFList = newDF[0].value_counts()
# reads2AlignsSameHap = newDFList[newDFList == 2].index
# reads2AlignsSameHapDF = reads2AlignDF[reads2AlignDF[0].isin(reads2AlignsSameHap)]
# reads2AlignsSameHapDF[0].nunique()

Purple= '#8B79A5'
Green ='#95BCA5'
darkGreen ='#325b38'
Blue ='#89B0D0'
Pink='#FFB6C1'

finalList = []

for read in uniqueList:
    read_start = read[0]
    read_end = read_start + read[1]
    block_start = read_start
    block_width = read[1]
    finalList.append([read_start, read_end, block_start, block_width,Purple,False])

for read in top2AlignsList:
    read_start = read[0]
    read_end = read_start + read[1]
    block_start = read_start
    block_width = read[1]
    finalList.append([read_start, read_end, block_start, block_width,darkGreen,False])

figureWidth=8
figureHeight=5
plt.figure(figsize=(figureWidth,figureHeight))
panelHeight=4
panelWidth=8
panel2 = plt.axes([0.1/figureWidth,1.7/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])

iBlue=(88/255,85/255,120/255)

finalList.sort(key=itemgetter(0))
countDict = {}
count=0
for ypos in range(1,len(finalList)):
    lastVal = 15900000  # HARDCODE
    for read in finalList:
        gstart,gend,blockstarts,blockwidths,readType,plotted=read[0],read[1],read[2],read[3],read[4],read[5]
        if plotted != True: # if not plotted yet, check if start coord is greater than curr end
            if gstart > lastVal: # if greater -> plot read
                # print(f'plotting rectangles {count}')
                rectangle = mplpatches.Rectangle((blockstarts, ypos),
                                                blockwidths,0.5,
                                                facecolor=readType,
                                                edgecolor='black',
                                                linewidth=0.05)
                panel2.add_patch(rectangle)
                lastVal = gend
                read[5] = True
                count+=1

panel2.set_xlim(15900000,17224298)  # HARDCODE
panel2.set_ylim(-1,350)
panel2.set_xlabel('Genomic Coordinate (Mb)')
# panel2.set_xticks(['57','58','59','60','61','62'])
# plt.title(f'{i} Read Coverage for {j}bp read size')
plt.savefig(output, dpi=2400)
plt.clf()