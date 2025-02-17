import pandas as pd
import os
from subprocess import Popen, PIPE
import gzip
from multiprocessing import Pool
import argparse
import pysam

# Need ReverseComplement function to make a strand-specific Bowtie2 index
def ReverseComplement(seq):
    seq = seq.upper()
    reversecomplement = ""
    for base in seq:
        if base == "A":
            reversecomplement = "T"+reversecomplement
        elif base == "T":
            reversecomplement = "A"+reversecomplement
        elif base == "G":
            reversecomplement = "C"+reversecomplement
        elif base == "C":
            reversecomplement = "G"+reversecomplement
        else:
            reversecomplement = base+reversecomplement
    return reversecomplement

def read_genome_fasta(fa):
    fa_dict = {}
    with open(fa) as f:
        for line in f:
            if line.startswith('>'):
                contig = line.split('>')[-1].strip()
                if contig not in fa_dict:
                    fa_dict[contig] = ''
                else:
                    raise KeyError("Redundant reference names in genome fasta: "+contig)
            else:
                fa_dict[contig] = fa_dict[contig]+line.strip()
    
    return fa_dict

# Create custom bowtie2 index
def build_custom_index_bowtie2(file_name, file_format='fasta', genome_fasta=None, flank_size=300):
    base = file_name.split('/')[-1].split('.')[0] + '_bowtie2'
    if 'index' not in os.listdir('./'):
        os.makedirs('index')
    
    # First option: provide fasta file (can be compressed)
    if file_format == 'fasta':
        if not file_name.endswith('gz'):
            args = ["bowtie2-build",file_name,'index/'+base]
            p = Popen(args)
        else:
            with open('fa.tmp','w') as fout:
                with gzip.open(file_name, 'rb') as f:
                    file_content = f.read()
                    for line in file_content:
                        fout.write(line)
            args = ["bowtie2-build","fa.tmp",'index/'+base]
            p = Popen(args)

    # Second option: provide a table of genome coordinates
    # Must have columns Chromosome, 5p boundary, 3p boundary, index should be CNAG IDs
    else:
        if file_format == 'tsv': df = pd.read_table(file_name, index_col=0)
        if file_format == 'csv': df = pd.read_csv(file_name, index_col=0)
        else:
            raise ValueError("Unrecognized file_format: "+file_format)

        fa = read_genome_fasta(genome_fasta)
            
        with open('fa.tmp','w') as fout:
            for ix, r in df.iterrows():
                chrom, start, end = r[['Chromosome','5p boundary','3p boundary']]
                start, end = (int(start),int(end))
                # Using reverse complement allows for stranded mapping with Bowtie2 downstream, important for minimizing incorrect mapping when gene flanks overlap
                if start < end:
                    lseq = ReverseComplement(fa[chrom][start-flank_size:start])
                    rseq = fa[chrom][end:end+flank_size]
                elif start > end:
                    lseq = fa[chrom][start:start+flank_size]
                    rseq = ReverseComplement(fa[chrom][end-flank_size:end])
                
                fout.write('>'+ix+'L\n')
                fout.write(lseq+'\n')
                fout.write('>'+ix+'R\n')
                fout.write(rseq+'\n')
        
        args = ["bowtie2-build","fa.tmp",'index/'+base]
        p = Popen(args)
    p.wait()
    if 'fa.tmp' in os.listdir('./'):
        os.remove('fa.tmp')
    return './index/'+base

#Function to add UMIs to read names
def add_UMI(UMI_fq, read1_fq, read2_fq):
    new_read1 = read1_fq.split('.f')[0]+'_umi.fastq.gz'
    new_read2 = read2_fq.split('.f')[0]+'_umi.fastq.gz'
    
    # Open new fastq for writing in zipped format
    with gzip.open(new_read1, 'wb') as fout:
        # Open UMI file
        with gzip.open(UMI_fq, 'rb') as umi:
            # Open read 1 file
            with gzip.open(read1_fq, 'rb') as f2:
                # Iterate through UMI and read 1 files together and count lines
                for n, (u, r1) in enumerate(zip(umi, f2)):
                    # 1st out of every 4 lines has read name - remember it from read 1
                    if n%4 == 0:
                        new_name = r1.strip().split()[0]
                    # 2nd out of ever 4 lines has read sequence - grab the UMI and add to read 1 name
                    elif n%4 == 1:
                        new_name = new_name+'_'.encode()+u[:8]+'\n'.encode() #[:8]+'\n'.encode() accounts for cases where there are extra bases from read into common adaptor seq
                        fout.write(new_name)
                        fout.write(r1)
                    else:
                        fout.write(r1)
                        
    with gzip.open(new_read2, 'wb') as fout:
        # Open UMI file
        with gzip.open(UMI_fq, 'rb') as umi:
            # Open read 1 file
            with gzip.open(read2_fq, 'rb') as f2:
                # Iterate through UMI and read 1 files together and count lines
                for n, (u, r2) in enumerate(zip(umi, f2)):
                    # 1st out of every 4 lines has read name - remember it from read 2
                    if n%4 == 0:
                        new_name = r2.strip().split()[0]
                    # 2nd out of ever 4 lines has read sequence - grab the UMI and add to read 2 name
                    elif n%4 == 1:
                        new_name = new_name+'_'.encode()+u[:8]+'\n'.encode() #[:8]+'\n'.encode() accounts for cases where there are extra bases from read into common adaptor seq
                        fout.write(new_name)
                        fout.write(r2)
                    else:
                        fout.write(r2)

    return (new_read1, new_read2)

#Function to trim reads with Cutadapt
def trim_filter(read1_fq, read2_fq):
    temp_out1 = read1_fq.split('.f')[0]+'_trim_temp.fastq'
    temp_out2 = read2_fq.split('.qf')[0]+'_trim_temp.fastq'
    out1 = read1_fq.split('.f')[0]+'_trim.fastq'
    out2 = read2_fq.split('.f')[0]+'_trim.fastq'
    log1_name = read1_fq.split('.f')[0]+'_trim1.log'
    log2_name = read1_fq.split('.f')[0]+'_trim2.log'
    
    # For downstream KO-seq reads, the 22 bp sequence CGGCCGCATCCCTGCATCCAAC from the 3' end of the NAT marker should be present
    # in all reads. The call -g "CGGCCGCATCCCTGCATCCAAC;required;o=22;e=0.2...AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;optional" looks
    # for all 22 bp of that sequence (e=0.2 tolerates 20% error rate) and also checks for (optional) Illumina adaptor sequence
    # at the 3' end of reads. The combination of --pair-filer=first and --discard-untrimmed tosses any reads that did not contain
    # the required 5' adaptor sequence; this is important for removing reads resulting from spurious PCR amplification in the
    # genome
    if 'down' in read1_fq.lower() or 'downstream' in read1_fq.lower():
        args1 = ["cutadapt", "-g", str("CGGCCGCATCCCTGCATCCAAC;required;o=22;e=0.2...AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;optional"), "-A", "GTTGGATGCAGGGATGCGGCCG", "-q", "20", "-m", "20", "--pair-filter=first", "--discard-untrimmed", "-o", temp_out1, "-p", temp_out2, read1_fq, read2_fq]
        p = Popen(args1, stdout=PIPE, stderr=PIPE)
        p.wait()
        with open(log1_name, 'w') as log:
            for line in p.stdout:
                log.write(line.decode())
            for line in p.stderr:
                log.write(line.decode())
                
    # Same description as for downstream junctions, except checks for presence of CGCGCCTAGCAGCGGATCCAAC 22 bp sequence at 5'
    # end of NAT marker that should be present in all reads from upstream junctions.
    elif 'up' in read1_fq.lower() or 'upstream' in read1_fq.lower():
        args1 = ["cutadapt", "-g", str("CGCGCCTAGCAGCGGATCCAAC;required;o=22;e=0.2...AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;optional"), "-A", "GTTGGATCCGCTGCTAGGCGCG", "-q", "20", "-m", "20", "--pair-filter=first", "--discard-untrimmed", "-o", temp_out1, "-p", temp_out2, read1_fq, read2_fq]
        p = Popen(args1, stdout=PIPE, stderr=PIPE)
        p.wait()
        with open(log1_name, 'w') as log:
            for line in p.stdout:
                log.write(line.decode())
            for line in p.stderr:
                log.write(line.decode())
    else:
        raise ValueError("Unclear if flank is downstream or upstream: "+read1_fq)
    
    #Extra round of cutadapt to remove any reads that are too short after initial trimming
    args2 = ["cutadapt", "-m", "20:20", "-o", out1, "-p", out2, temp_out1, temp_out2]
    p = Popen(args2, stdout=PIPE, stderr=PIPE)
    p.wait()
    with open(log2_name, 'w') as log:
        for line in p.stdout:
            log.write(line.decode())
        for line in p.stderr:
            log.write(line.decode())
    
    return (out1, out2)

#Function to map reads with Bowtie2
def call_bowtie2(mate, index, threads=1):
    sam_name = mate[0].split('_umi_trim')[0]+'.sam'
    log_name = mate[0].split('_umi_trim')[0]+'_bowtie2.log'
    
    # Combination of --fr and --norc is important to ensure that reads align to KO-junctions in the expected direction
    args = ['bowtie2', '--end-to-end', '--fr', '--norc', '--no-mixed', '--no-discordant', '-X500', '--score-min C,-22,0', '-p'+str(threads), '-x'+index, '-1', mate[0], '-2', mate[1], '-S', sam_name]
            
    p = Popen(args, stdout=PIPE, stderr=PIPE)
    p.wait()
    with open(log_name, 'w') as log:
        for line in p.stdout:
            log.write(line.decode())
        for line in p.stderr:
            log.write(line.decode())
    return sam_name

#Function to convert sam files to bam, sort and index
def sort_and_index_bam(sam):
    bam_name = sam.split('.sam')[0]+'.bam'
    sorted_bam_name = sam.split('.sam')[0]+'.sorted.bam'
    
    args1 = ["samtools", "view", "-bSo", bam_name, sam]
    p = Popen(args1)
    p.wait()
    
    args2 = ["samtools", "sort", bam_name, "-o", sorted_bam_name]
    p = Popen(args2)
    p.wait()
    
    args3 = ["samtools", "index", sorted_bam_name]
    p = Popen(args3)
    p.wait()
    
    return sorted_bam_name

#Function to deduplicate reads based on UMI and read 2 mapping coordinates with UMI-tools
def deduplicate_reads(sorted_bam):
    base = sorted_bam.split('.sorted.bam')[0]+'_dedup'
    deduplicated_bam = base+'.sorted.bam'
    
    args1 = ["umi_tools", "dedup", "--paired", "--output-stats="+base, "--log="+base+".log", "--stdin="+sorted_bam, "--stdout="+deduplicated_bam]
    p = Popen(args1)
    p.wait()
    
    args2 = ["samtools", "index", deduplicated_bam]
    p = Popen(args2)
    p.wait()
    
    return (sorted_bam, deduplicated_bam)

#Function to count reads, export both deduplicated and non-deduplicated output files
def count_flank_reads(nondedup_bam, dedup_bam):
    nondedup_bam_file = pysam.AlignmentFile(nondedup_bam)
    dedup_bam_file = pysam.AlignmentFile(dedup_bam)

    nondedup_flank_dict = {}
    dedup_flank_dict = {}
    
    for read in nondedup_bam_file.fetch():
        if read.is_read1:
            
            # Downstream reads should have 'R' in ref names, upstream should have 'L'. Ignore cases where that's not true.
            if '_down_' in nondedup_bam.lower() and 'L' in read.reference_name:
                continue
            if '_up_' in nondedup_bam.lower() and 'R' in read.reference_name:
                continue
                
            # Ignore cases where read starts 10 bp or more into the KO junction, as these are likely spurious PCR amplifications
            # from gDNA as opposed to true KO's
            if read.reference_start >= 10:
                continue
            if not read.has_tag('XS'):
                if read.reference_name.strip() not in nondedup_flank_dict.keys():
                    nondedup_flank_dict[read.reference_name.strip()] = 1
                else:
                    nondedup_flank_dict[read.reference_name.strip()] += 1
            if read.has_tag('XS'):
                if read.get_tag('AS') > read.get_tag('XS'):
                    if read.reference_name.strip() not in nondedup_flank_dict.keys():
                        nondedup_flank_dict[read.reference_name.strip()] = 1
                    else:
                        nondedup_flank_dict[read.reference_name.strip()] += 1
    
    nondedup_bam_file.close()   

    for read in dedup_bam_file.fetch():
        if read.is_read1:
            if '_down_' in dedup_bam.lower() and 'L' in read.reference_name:
                continue
            if '_up_' in dedup_bam.lower() and 'R' in read.reference_name:
                continue
            if read.reference_start >= 10:
                continue
            if not read.has_tag('XS'):
                if read.reference_name.strip() not in dedup_flank_dict.keys():
                    dedup_flank_dict[read.reference_name.strip()] = 1
                else:
                    dedup_flank_dict[read.reference_name.strip()] += 1
            if read.has_tag('XS'):
                if read.get_tag('AS') > read.get_tag('XS'):
                    if read.reference_name.strip() not in dedup_flank_dict.keys():
                        dedup_flank_dict[read.reference_name.strip()] = 1
                    else:
                        dedup_flank_dict[read.reference_name.strip()] += 1
    
    dedup_bam_file.close()
    
    nondedup_out_file = nondedup_bam.split('.')[0]+'_nondedup_counts.csv'
    dedup_out_file = dedup_bam.split('.')[0]+'_counts.csv'
    
    with open(nondedup_out_file, 'w') as fout:
        fout.write('Gene Flank,Gene,Count\n')
        
        for flank in nondedup_flank_dict.keys():
            fout.write(flank+','+flank[0:10]+','+str(nondedup_flank_dict[flank])+'\n')
    
    with open(dedup_out_file, 'w') as fout:
        fout.write('Gene Flank,Gene,Count\n')
        
        for flank in dedup_flank_dict.keys():
            fout.write(flank+','+flank[0:10]+','+str(dedup_flank_dict[flank])+'\n')
    
    return (nondedup_out_file, dedup_out_file)

#Function to remove KO's that failed QC
def remove_incorrect_KOs(nondedup_csv, dedup_csv):
    true_file = './KO_true_list.csv'
    fail_file = './KO_fail_list.csv'

    true_list = []
    fail_list = []

    with open(true_file, 'r', encoding = 'utf-8-sig') as true:
        for gene in true:
            true_list.append(gene.strip())

    with open(fail_file, 'r', encoding = 'utf-8-sig') as fail:
        for gene in fail:
            fail_list.append(gene.strip())

    nondedup_filtered_csv = nondedup_csv.split('.csv')[0]+'_filtered.csv'
    dedup_filtered_csv = dedup_csv.split('.csv')[0]+'_filtered.csv'
    
    with open(nondedup_csv, 'r') as in_file:
        with open(nondedup_filtered_csv, 'w') as out_file:
            out_file.write('Gene Flank,Gene,Count'+'\n')
            for line in in_file:
                gene = line.split(',')[1]
                if gene in true_list:
                    if not gene in fail_list:
                        out_file.write(line.strip()+'\n')
    
    with open(dedup_csv, 'r') as in_file:
        with open(dedup_filtered_csv, 'w') as out_file:
            out_file.write('Gene Flank,Gene,Count'+'\n')
            for line in in_file:
                gene = line.split(',')[1]
                if gene in true_list:
                    if not gene in fail_list:
                        out_file.write(line.strip()+'\n')
                                
                                
    return (nondedup_filtered_csv, dedup_filtered_csv)

#Function to remove unnecessary intermediate files
def remove_files(directory):
    for f in os.listdir(directory):
        if f.endswith('_umi.fastq.gz') or f.endswith('_umi.fq.gz'): os.remove(directory+f)
        if f.endswith('_trim_temp.fastq') or f.endswith('_trim_temp.fq'): os.remove(directory+f)
        if f.endswith('_umi_trim.fastq') or f.endswith('_umi_trim.fq'): os.remove(directory+f)
        if f.endswith('.sam'): os.remove(directory+f)
    return "Finished"

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(help='sub-command help', dest='subcommand')

    parser_index = subparsers.add_parser('build_index', help='Create Bowtie2 index')
    parser_index.add_argument("--input", help="Sequence file for creating Bowtie2 index")
    parser_index.add_argument("--format", default='fasta', help="Format of sequence file")
    parser_index.add_argument("--genome_fasta", default=None, help="Whole genome fasta for building index from coordinate table")
    parser_index.add_argument("--flank_size", type=int, default=300, help="Size of flanking sequences to use for index")
    
    parser_align = subparsers.add_parser('align', help='Align and count reads')
    parser_align.add_argument("--index", help="Bowtie2 index prefix (from build_index option)")
    parser_align.add_argument("--directory", default='./', help="Working directory containing fastq files")
    parser_align.add_argument("--threads", type=int, default=1, help="Number of processors")
    
    args = parser.parse_args()
    
    if args.subcommand == 'build_index':
        index_prefix = build_custom_index_bowtie2(args.input, file_format=args.format, 
                                          genome_fasta=args.genome_fasta, flank_size=args.flank_size)
        print("Index created at: "+index_prefix)

    elif args.subcommand == 'align':
        #Find all file sets in the directory
        fq_sets = []
        for file in os.listdir(args.directory):
            if file.endswith('fq.gz') or file.endswith('fastq.gz'):
                
                # Note: our Illumina outputs denote read 1 as R1, the UMI at the i7 index site as R2, and read 2 as R3.
                # This might differ between platforms or facilities.
                if '_R1' in file and 'undetermined_' not in file.lower():
                    r1 = file
                    r2 = r1.split('_R1')[0]+'_R3'+r1.split('_R1')[1]
                    umi = r1.split('_R1')[0]+'_R2'+r1.split('_R1')[1]
                    fq_sets.append((args.directory+umi,args.directory+r1,args.directory+r2)) 
        
        # Add UMI
        print("Adding UMI to read names...\n")
        p = Pool(args.threads)
        barcoded = p.starmap(add_UMI, fq_sets)
        
        # Trim with cutadapt
        print("Trimming and filtering with Cutadapt...\n")
        trimmed = p.starmap(trim_filter, barcoded)
        
        # Align with Bowtie2
        print("Aligning with Bowtie2...")
        sams = []
        for pair in trimmed:
            sam = call_bowtie2(pair, args.index, threads=args.threads)
            sams.append(sam)
        
        # Convert sam to bam, sort bam, index bam
        print("Generating, sorting, and indexing bam file...")
        p = Pool(args.threads)
        bams = p.map(sort_and_index_bam, sams)
        
        # Deduplicate reads with UMI-tools package
        print("Deduplicating reads...")
        p = Pool(args.threads)
        dedup_bams = p.map(deduplicate_reads, bams)
        
        # Count reads for each flank
        print("Counting reads...")
        p = Pool(args.threads)
        counts = p.starmap(count_flank_reads, dedup_bams)
        
        # Remove KOs not actually in collection or that failed QC
        print("Removing incorrect KOs...")
        p = Pool(args.threads)
        filtered = p.starmap(remove_incorrect_KOs, counts)
        
        remove_files(args.directory)

if __name__ == '__main__':
    main()