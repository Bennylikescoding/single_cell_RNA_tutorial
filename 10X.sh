### 从bcl->fq->reads_counts

# bcl2fasq (optional, sudo)

wget -c https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2-20-0-linux-x86-64.zip
unzip bcl2fastq2-v2-20-0-linux-x86-64.zip
rpm -i bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm 

# Cell ranger
wget -O cellranger-3.0.2.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.2.tar.gz?Expires=1550666525&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTMuMC4yLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU1MDY2NjUyNX19fV19&Signature=kR9wRykho-yCaWidpGAJe1dLuk9Is7CTu1aVxeul74joIEBBfH2ZFz2BJ8DEqhbepsOnt04uhS0c3YWqr6Uc8oSbTwkRWEjha1WNjvMDWuN8rZwcRj4ELxPYOZOzwwV4uGcmt34qE9nsm3utVtoT97b2wJ0H-t6PVBoBTpbk-hl3zIUC3vKBE-2YIhyiHbmawu2wiiTCP2QQXAw3Wl1Z0ovDTBcOvvmu2SCnXHfEqAa1nKTMNRWrthQgqeTuJx7d9ctd20bA2yp16Mfc8LyJSi5JYdh20g-rM~-vzUqYBWGfzllBeN65U8NceNGkBtJpT3CwHkp45yvz0sxFn~hiJg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

tar xvzf cellranger-3.0.2.tar.gz

ln -s `pwd`/cellranger-3.0.2/cellranger-cs/3.0.2/bin/cellranger ~/soft/bin

# 以Cellranger自带的fastq和genome为例
# 实际分析时把自己的测序文件放入fastq文件夹和基因组放入genome文件夹
mkdir -p fastq genome
cp `pwd`/cellranger-3.0.2/cellranger-tiny-fastq/3.0.0/*.fastq.gz fastq
cp `pwd`/cellranger-3.0.2/cellranger-tiny-ref/3.0.0/genes/genes.gtf genome
cp `pwd`/cellranger-3.0.2/cellranger-tiny-ref/3.0.0/fasta/genome.fa genome

# Test run
cellranger testrun --id=tiny


# 测试数据

wget -c http://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
wget -c http://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv
wget -c http://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-samplesheet-1.2.0.csv
mkdir -p bcl
mv cellranger-tiny-bcl* bcl/



# GTF过滤
# 从ENSEMBL或UCSC获得的GTF最好过滤下，去掉不关注的类型的转录本如pseudogene等。
# 采用Key-value对的形式进行过滤。
# 过滤是为了减少可能的基因重叠。
# GTF中必须包含 feature type (`exon`)，第三列
cd ~/scRNA
cellranger mkgtf genome/GRCh38.gtf genome/GRCh38.filtered.gtf \
                   --attribute=gene_biotype:protein_coding \
                   --attribute=gene_biotype:lincRNA \
                   --attribute=gene_biotype:antisense \
                   --attribute=gene_biotype:IG_LV_gene \
                   --attribute=gene_biotype:IG_V_gene \
                   --attribute=gene_biotype:IG_V_pseudogene \
                   --attribute=gene_biotype:IG_D_gene \
                   --attribute=gene_biotype:IG_J_gene \
                   --attribute=gene_biotype:IG_J_pseudogene \
                   --attribute=gene_biotype:IG_C_gene \
                   --attribute=gene_biotype:IG_C_pseudogene \
                   --attribute=gene_biotype:TR_V_gene \
                   --attribute=gene_biotype:TR_V_pseudogene \
                   --attribute=gene_biotype:TR_D_gene \
                   --attribute=gene_biotype:TR_J_gene \
                   --attribute=gene_biotype:TR_J_pseudogene \
                   --attribute=gene_biotype:TR_C_gene


#cellranger mkgtf (3.0.2)
#Copyright (c) 2019 10x Genomics, Inc.  All rights reserved.
#-------------------------------------------------------------------------------
#
#Writing new genes GTF file (may take 10 minutes for a 1GB input GTF file)...
#...done


# 建立索引##############------------更改路径-------------------------->>>>>>>>>>>.
# 结果输出在--genome指定的目录下
# Ref: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references
cd ~/scRNA
##############------------更改路径-------------------------->>>>>>>>>>>.
cellranger mkref --genome=GRCh8_10X_index --fasta=genome/GRCh38.fa --genes=genome/GRCh38.filtered.gtf


#?cellranger mkref (3.0.2)
#?Copyright (c) 2019 10x Genomics, Inc.  All rights reserved.
#?-------------------------------------------------------------------------------
#?
#?Creating new reference folder at /disk2/home/RNA_002/scRNA/GRCh8_10X_index
#?...done
#?
#?Writing genome FASTA file into reference folder...
#?...done
#?
#?Computing hash of genome FASTA file...
#?...done
#?
#?Indexing genome FASTA file...
#?...done
#?
#?Writing genes GTF file into reference folder...
#?...done
#?
#?Computing hash of genes GTF file...
#?...done
#?
#?Writing genes index file into reference folder (may take over 10 minutes for a 3Gb genome)...
#?...done
#?
#?Writing genome metadata JSON file into reference folder...
#?...done
#?
#?Generating STAR genome index (may take over 8 core hours for a 3Gb genome)...
#?Feb 20 10:55:33 ..... Started STAR run
#?Feb 20 10:55:33 ... Starting to generate Genome files
#?Feb 20 10:55:35 ... starting to sort  Suffix Array. This may take a long time...
#?Feb 20 10:55:36 ... sorting Suffix Array chunks and saving them to disk...
#?Feb 20 10:57:20 ... loading chunks from disk, packing SA...
#?Feb 20 10:57:22 ... Finished generating suffix array
#?Feb 20 10:57:22 ... Generating Suffix Array index
#?Feb 20 10:57:25 ... Completed Suffix Array index
#?Feb 20 10:57:25 ..... Processing annotations GTF
#?Feb 20 10:57:26 ..... Inserting junctions into the genome indices
#?Feb 20 10:57:32 ... writing Genome to disk ...
#?Feb 20 10:57:32 ... writing Suffix Array to disk ...
#?Feb 20 10:57:32 ... writing SAindex to disk
#?Feb 20 10:57:32 ..... Finished successfully
#?...done.
#?
#?>>> Reference successfully created! <<<
#?
#?You can now specify this reference on the command line:
#?cellranger --transcriptome=/disk2/home/RNA_002/scRNA/GRCh8_10X_index ...

# BCL demultiplexing

## --run: bcl directory
## --csv: simple csv
## --id: 输出文件的名字
## --qc: 质控结果在outs/qc_summary.json
##############------------更改路径-------------------------->>>>>>>>>>>.
cd ~/scRNA
/bin/rm -rf ysx-fastq
(cd bcl; tar xvzf cellranger-tiny-bcl-1.2.0.tar.gz)
cellranger mkfastq --id=ysx-fastq --run=bcl/cellranger-tiny-bcl-1.2.0 --csv=bcl/cellranger-tiny-bcl-simple-1.2.0.csv --delete-undetermined --qc


## cellranger mkfastq (3.0.2)
## Copyright (c) 2019 10x Genomics, Inc.  All rights reserved.
## -------------------------------------------------------------------------------
## 
## Martian Runtime - '3.0.2-v3.2.0'
## 2019-02-20 12:22:54 [jobmngr] WARNING: The current virtual address space size
##                               limit is too low.
## 	Limiting virtual address space size interferes with the operation of many
## 	common libraries and programs, and is not recommended.
## 	Contact your system administrator to remove this limit.
## 2019-02-20 12:22:54 [runtime] Reattaching in local mode.
## Serving UI at http://localhost.localdomain:33112?auth=jpYlAetrOJgg1vi7EoFcW-ob4fDg8v7vtnyhpcK86nY
## 
## 2019-02-20 12:22:54 [runtime] (reset-partial)   ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MAKE_FASTQS_PREFLIGHT.fork0.chnk0
## 2019-02-20 12:22:54 [runtime] (reset-partial)   ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MAKE_FASTQS_PREFLIGHT_LOCAL.fork0.chnk0
## 2019-02-20 12:22:54 [runtime] Found orphaned local stage: ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MAKE_FASTQS_PREFLIGHT
## 2019-02-20 12:22:54 [runtime] Found orphaned local stage: ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MAKE_FASTQS_PREFLIGHT_LOCAL
## Checking run folder...
## Checking RunInfo.xml...
## Checking system environment...
## Checking barcode whitelist...
## Emitting run information...
## Checking read specification...
## Checking samplesheet specs...
## Checking for dual index flowcell...
## 2019-02-20 12:22:54 [runtime] (ready)           ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.PREPARE_SAMPLESHEET
## 2019-02-20 12:22:54 [runtime] (run:local)       ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.PREPARE_SAMPLESHEET.fork0.chnk0.main
## 2019-02-20 12:22:55 [runtime] (chunks_complete) ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.PREPARE_SAMPLESHEET
## 2019-02-20 12:22:55 [runtime] (ready)           ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.BCL2FASTQ_WITH_SAMPLESHEET
## 2019-02-20 12:22:55 [runtime] (run:local)       ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.BCL2FASTQ_WITH_SAMPLESHEET.fork0.split
## 2019-02-20 12:22:55 [runtime] (split_complete)  ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.BCL2FASTQ_WITH_SAMPLESHEET
## 2019-02-20 12:22:55 [runtime] (run:local)       ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.BCL2FASTQ_WITH_SAMPLESHEET.fork0.chnk0.main
## 2019-02-20 12:23:04 [runtime] (chunks_complete) ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.BCL2FASTQ_WITH_SAMPLESHEET
## 2019-02-20 12:23:04 [runtime] (run:local)       ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.BCL2FASTQ_WITH_SAMPLESHEET.fork0.join
## 2019-02-20 12:23:05 [runtime] (join_complete)   ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.BCL2FASTQ_WITH_SAMPLESHEET
## 2019-02-20 12:23:05 [runtime] (ready)           ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MAKE_QC_SUMMARY
## 2019-02-20 12:23:05 [runtime] (run:local)       ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MAKE_QC_SUMMARY.fork0.split
## 2019-02-20 12:23:05 [runtime] (split_complete)  ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MAKE_QC_SUMMARY
## 2019-02-20 12:23:06 [runtime] (run:local)       ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MAKE_QC_SUMMARY.fork0.join
## 2019-02-20 12:23:06 [runtime] (join_complete)   ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MAKE_QC_SUMMARY
## 2019-02-20 12:23:06 [runtime] (ready)           ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MERGE_FASTQS_BY_LANE_SAMPLE
## 2019-02-20 12:23:06 [runtime] (run:local)       ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MERGE_FASTQS_BY_LANE_SAMPLE.fork0.split
## 2019-02-20 12:23:07 [runtime] (split_complete)  ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MERGE_FASTQS_BY_LANE_SAMPLE
## 2019-02-20 12:23:07 [runtime] (run:local)       ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MERGE_FASTQS_BY_LANE_SAMPLE.fork0.chnk0.main
## 2019-02-20 12:23:07 [runtime] (chunks_complete) ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MERGE_FASTQS_BY_LANE_SAMPLE
## 2019-02-20 12:23:07 [runtime] (run:local)       ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MERGE_FASTQS_BY_LANE_SAMPLE.fork0.join
## 2019-02-20 12:23:08 [runtime] (join_complete)   ID.ysx-fastq.MAKE_FASTQS_CS.MAKE_FASTQS.MERGE_FASTQS_BY_LANE_SAMPLE
## 
## Outputs:
## - Run QC metrics:        null
## - FASTQ output folder:   /disk2/home/RNA_002/scRNA/ysx-fastq/outs/fastq_path
## - Interop output folder: /disk2/home/RNA_002/scRNA/ysx-fastq/outs/interop_path
## - Input samplesheet:     /disk2/home/RNA_002/scRNA/ysx-fastq/outs/input_samplesheet.csv
## 
## Waiting 6 seconds for UI to do final refresh.
## Pipestance completed successfully!
## 
## 2019-02-20 12:23:14 Shutting down.
## Saving pipestance info to ysx-fastq/ysx-fastq.mri.tgz


# 10x count

# 拆分后的数据在 ysx-fastq/outs/fastq_path/H35KCBCXY/ysx/ 目录下，但这个数据量太少，只有88条reads，下面换一个测试数据在fastq目录下

cd ~/scRNA

# cell ranger默认使用所有的核，这里用--localcores做个限制
# --localmem 限制最大使用多少G内存
# 结果输出在 --id指定的目录下
# --expect-cells: default 3000
# 需要15分钟
##############------------更改路径-------------------------->>>>>>>>>>>.
cellranger count --id=ysx --transcriptome=GRCh8_10X_index --sample=tinygex --fastqs=fastq/ --localcores=2 --localmem=10

# Copyright (C) 2018 Genome Research Ltd.
# 
# 2019-02-20 13:52:32 [runtime] (ready)           ID.ysx.SC_RNA_COUNTER_CS.SC_RNA_COUNTER.DISABLE_FEATURE_STAGES
# 2019-02-20 13:52:32 [runtime] (run:local)       ID.ysx.SC_RNA_COUNTER_CS.SC_RNA_COUNTER.DISABLE_FEATURE_STAGES.fork0.chnk0.main
# Outputs:
# - Run summary HTML:                         /disk2/home/RNA_002/scRNA/ysx/outs/web_summary.html
# - Run summary CSV:                          /disk2/home/RNA_002/scRNA/ysx/outs/metrics_summary.csv
# - BAM:                                      /disk2/home/RNA_002/scRNA/ysx/outs/possorted_genome_bam.bam
# - BAM index:                                /disk2/home/RNA_002/scRNA/ysx/outs/possorted_genome_bam.bam.bai
# - Filtered feature-barcode matrices MEX:    /disk2/home/RNA_002/scRNA/ysx/outs/filtered_feature_bc_matrix
# - Filtered feature-barcode matrices HDF5:   /disk2/home/RNA_002/scRNA/ysx/outs/filtered_feature_bc_matrix.h5
# - Unfiltered feature-barcode matrices MEX:  /disk2/home/RNA_002/scRNA/ysx/outs/raw_feature_bc_matrix
# - Unfiltered feature-barcode matrices HDF5: /disk2/home/RNA_002/scRNA/ysx/outs/raw_feature_bc_matrix.h5
# - Secondary analysis output CSV:            /disk2/home/RNA_002/scRNA/ysx/outs/analysis
# - Per-molecule read information:            /disk2/home/RNA_002/scRNA/ysx/outs/molecule_info.h5
# - CRISPR-specific analysis:                 null
# - Loupe Cell Browser file:                  /disk2/home/RNA_002/scRNA/ysx/outs/cloupe.cloupe
# 
# Waiting 6 seconds for UI to do final refresh.
# Pipestance completed successfully!
# 
# 2019-02-20 14:07:28 Shutting down.
# Saving pipestance info to ysx/ysx.mri.tgz


# 查看count的输出结果

head ysx/outs/metrics_summary.csv

# 更好的查看格式
cat <<END | R --vanilla -q --slave
data <- read.table("ysx/outs/metrics_summary.csv", sep=",")
print(t(data))
END

#     [,1]                                             [,2]     
# V1  "Estimated Number of Cells"                      "1,107"  
# V2  "Mean Reads per Cell"                            "416"    
# V3  "Median Genes per Cell"                          "19"     
# V4  "Number of Reads"                                "461,083"
# V5  "Valid Barcodes"                                 "94.8%"  
# V6  "Sequencing Saturation"                          "69.4%"  
# V7  "Q30 Bases in Barcode"                           "94.5%"  
# V8  "Q30 Bases in RNA Read"                          "91.7%"  
# V9  "Q30 Bases in Sample Index"                      "91.8%"  
# V10 "Q30 Bases in UMI"                               "93.3%"  
# V11 "Reads Mapped to Genome"                         "99.9%"  
# V12 "Reads Mapped Confidently to Genome"             "84.8%"  
# V13 "Reads Mapped Confidently to Intergenic Regions" "5.8%"   
# V14 "Reads Mapped Confidently to Intronic Regions"   "43.7%"  
# V15 "Reads Mapped Confidently to Exonic Regions"     "35.3%"  
# V16 "Reads Mapped Confidently to Transcriptome"      "33.6%"  
# V17 "Reads Mapped Antisense to Gene"                 "0.6%"   
# V18 "Fraction Reads in Cells"                        "94.1%"  
# V19 "Total Genes Detected"                           "201"    
# V20 "Median UMI Counts per Cell"                     "29"


head ysx/outs/analysis/tsne/2_components/projection.csv

# Barcode,TSNE-1,TSNE-2
# AAACCCAAGGAGAGTA-1,-6.307645823966444,-21.61306482414778
# AAACGCTTCAGCCCAG-1,-10.020325222302018,9.689193966800287
# AAAGAACAGACGACTG-1,8.42357866010068,2.456971311396624
# AAAGAACCAATGGCAG-1,-2.7252532873977158,16.05445011680964
# AAAGAACGTCTGCAAT-1,7.945601833376288,11.63161709287395
# AAAGGATAGTAGACAT-1,-11.233212336507123,12.487366744008591
# AAAGGATCACCGGCTA-1,-0.65670826911923,15.061746386644229
# AAAGGATTCAGCTTGA-1,-15.937712546664537,-13.061572418865026
# AAAGGATTCCGTTTCG-1,-3.1954679996017807,-18.105088108147974


head ysx/outs/analysis/clustering/graphclust/clusters.csv

# Barcode,Cluster
# AAACCCAAGGAGAGTA-1,1
# AAACGCTTCAGCCCAG-1,4
# AAAGAACAGACGACTG-1,3
# AAAGAACCAATGGCAG-1,2
# AAAGAACGTCTGCAAT-1,2
# AAAGGATAGTAGACAT-1,4
# AAAGGATCACCGGCTA-1,2
# AAAGGATTCAGCTTGA-1,1
# AAAGGATTCCGTTTCG-1,1


# 过滤empty beads后的表达矩阵
ls ysx/outs/filtered_feature_bc_matrix

zcat ysx/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | head

# AAACCCAAGGAGAGTA-1
# AAACGCTTCAGCCCAG-1
# AAAGAACAGACGACTG-1
# AAAGAACCAATGGCAG-1
# AAAGAACGTCTGCAAT-1
# AAAGGATAGTAGACAT-1
# AAAGGATCACCGGCTA-1
# AAAGGATTCAGCTTGA-1
# AAAGGATTCCGTTTCG-1
# AAAGGGCTCATGCCCT-1

zcat ysx/outs/filtered_feature_bc_matrix/features.tsv.gz | head

# ENSG00000279493	CH507-9B2.2	Gene Expression
# ENSG00000277117	CH507-9B2.1	Gene Expression
# ENSG00000279687	CH507-9B2.8	Gene Expression
# ENSG00000280071	CH507-9B2.3	Gene Expression
# ENSG00000276612	CH507-9B2.4	Gene Expression
# ENSG00000275464	CH507-9B2.5	Gene Expression
# ENSG00000280433	CH507-9B2.9	Gene Expression
# ENSG00000279669	RP11-717F1.1	Gene Expression
# ENSG00000279094	CH507-24F1.1	Gene Expression
# ENSG00000274333	RP11-717F1.2	Gene Expression


zcat ysx/outs/filtered_feature_bc_matrix/matrix.mtx.gz | head

# %%MatrixMarket matrix coordinate integer general
# %metadata_json: {"format_version": 2, "software_version": "3.0.2"}
# 507 1107 23866
# 458 1 3
# 456 1 1
# 409 1 1
# 406 1 2
# 391 1 2
# 386 1 1
# 376 1 1

