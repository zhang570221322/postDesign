#!/bin/bash
### Set job name for example MyPBS
#PBS -N for build_db of my classfication

### set output files
#PBS -o MyPBS.stdout
#PBS -e MyPBS.stderr

### set queue name
#PBS -q batch

###set number of nodes,for example one node with 4 cpu.And you can specify the node: #PBS -l nodes=cu01:ppn=4

#PBS -l nodes=cu02:ppn=28
#the following is you own code
# loacl var
thread=4

jellyfish=/home/wlzhang/classfication/jellyfish-linux
viral_assembly_summary=/home/wlzhang/classfication/refseq/viral/assembly_summary_refseq.txt


k_number=$1
database_fna="/home/wlzhang/classfication/refseq/viral_database/kmer_count_k"$k_number
rm -f $database_fna
touch $database_fna
for fasta_gz_name  in  $(find /home/wlzhang/classfication/refseq/viral  -name '*.fna.gz')
    do
    # get filename
    temp=${fasta_gz_name##*/}
    # get sequence id
    array=(${temp//_/ })
    # joint GCF_ID
    GCF_id=${array[0]}"_"${array[1]}
    # find taxid
    assembly_summary_information=$(awk -F "\t"  -v target=$GCF_id '$1==target {print "taxid:|:"$6}'  $viral_assembly_summary)
    echo "find  $assembly_summary_information"
    # if  assembly_summary exists
    if [  -n "$assembly_summary_information" ]; then
        echo "[START] ***  $assembly_summary_information"
        echo $assembly_summary_information >> $database_fna
        zcat  $fasta_gz_name | $jellyfish count -m $k_number -s 100M -t $thread -o  /dev/fd/1   /dev/fd/0 | $jellyfish dump -c -t  /dev/fd/0  >>   $database_fna
        # zcat  $fasta_gz_name > ./temp_fa
        # $jellyfish count -m $k_number -s 100M -t $thread -o  ./temp ./temp_fa
        # $jellyfish dump -c -t  ./temp   >>   $database_fna
        echo "[DONE]  ***  $assembly_summary_information"
    fi
    done
echo "DONE kmer $k_number"


# ??
nohup  ./build_db.sh 20 > log20.out 2>&1 &
# test
jellyfish=/home/wlzhang/classfication/jellyfish-linux
fasta_gz_name=/home/wlzhang/classfication/refseq/viral/GCF_000869605.1/GCF_000869605.1_ViralProj18303_genomic.fna.gz
k_number=26
thread=10
database_fna="/home/wlzhang/classfication/refseq/viral_database/kmer_count_k"$k_number
zcat  $fasta_gz_name > ./temp_fa
$jellyfish count -m $k_number -s 100M -t $thread -o  ./temp ./temp_fa    
$jellyfish dump -c -t  ./temp   >>   $database_fna

while true;do a=1+1;done