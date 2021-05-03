wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxcat.tar.gz
tar zxvf taxcat.tar.gz
rm taxcat.tar.gz
wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar zxvf taxdump.tar.gz
rm taxdump.tar.gz
# 下载完成分类树
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{category}/assembly_summary.txt
# 下载assembly_summary_refse 所有的


jellyfish count -m 16 -s 100M -t 24 -o mer_counts   
#  k=16 , hash数组大小为100M  线程24 

# 多个taxid相同的
ls *.fna.gz | xargs -n 1 echo gunzip -c > generators 
sudo jellyfish count -g generators -m 16 -s 10M
# $fasta_gz_name"_k"$k_number
# 批处理
viral_assembly_summary=/mnt/e/Github/classification/refSeq/assembly_summary/viral/assembly_summary.txt
# $(awk -F "\t"  -v target=$GCF_id '$1==target {print "taxid:"$6",species_taxid:"$7",name:"$8}'  $viral_assembly_summary)
thread=8
basename="/mnt/e/Github/classification/refSeq/pipline/raw_data"
for k_number in 10 12 14
do 
database_fna="/mnt/e/Github/classification/refSeq/test/database_fna_k"$k_number
rm -f $database_fna 
touch $database_fna
for fasta_gz_name  in  $(find $basename  -name '*.fna.gz') 
    do
    # 取到文件名
    temp=${fasta_gz_name##*/}
    # 取到文件id号
    array=(${temp//_/ })
    # 拼接GCF_ID
    GCF_id=${array[0]}"_"${array[1]}
    # 寻找taxid
    assembly_summary_information=$(awk -F "\t"  -v target=$GCF_id '$1==target {print "taxid:"$6",species_taxid:"$7",id:"$1",name:"$8}'  $viral_assembly_summary) 
    # 如果找到assembly_summary
    echo $assembly_summary_information 
    if [  -n "$assembly_summary_information" ]; then 
        echo $assembly_summary_information >> $database_fna
        zcat  $fasta_gz_name | jellyfish count -m $k_number -s 100M -t $thread -o  /dev/fd/1   /dev/fd/0 | jellyfish dump -c -t  /dev/fd/0  >>   $database_fna
        echo "DONE $fasta_gz_name "
    fi  
    done
done 
a=$(find ../viral  -name '*.fna.gz' | wc -l) ;b=$(cat database_fna_k26 | grep taxid | wc -l); c=$(echo "scale=5;$b/$a*100" | bc ) ;echo "${c}%"