library=/home/cszhang/workSpace/classification/refSeq/pipline/for_kraken_db/library.removeNX.fna
grinder -reference_file $library -total_reads 100000 -diversity 100 -rd 150
ssh 49.235.23.28 -p 6062 -l wlzhang
nohup  python Union_Counter_save.py  > Union_Counter_save.py.log 2>Union_Counter_save.py.log &
/home/kyy/kraken/kraken2/kraken2 -db  /home/kyy/kraken/kraken2_viral_database3 ./grinder-reads.fa --report grinder-reads.fa.report --output grinder-reads.fa.output
find .  -name 'library.fna28_union_*' | xargs du -ck
C:/ProgramData/Anaconda3/Scripts/activate
jellyfish=/home/wlzhang/classfication/jellyfish-2.3.0/bin/jellyfish

# 显示https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi 的taxid
function s(){
for(i=0;i< $('#dummybodyid > table:nth-child(7)').getElementsByTagName('a').length;i++){
     a = $('#dummybodyid > table:nth-child(7)').getElementsByTagName('a')[i];
     taxid=a.getAttribute('href').split('id=')[1];
    old_innterText=a.innerText;
    a.innerText=old_innterText+":"+taxid;
}
}
# 下载SRA数据
# nohup  prefetch --option-file sra_list -O ./  >> prefetch.log  2>>prefetch.log &
# 查看多个文件总大小
ll library.removeNX.fna28_union_* | awk -F ' ' '{sum+=$5};END{print sum/1024/1024}'
# 解压
/opt/anaconda2/bin/fasterq-dump
# (wc -l library.removeNX.fna28_union_SPECIES$thread ; wc -l library.removeNX.fna28_union_LCA$thread )| awk -F ' ' '{sum+=$1};END{print sum/11545173}'

for((thread=0;thread<=23;thread++))
do 
echo ${thread}; 
(wc -l library.removeNX.fna28_union_SPECIES$thread ; wc -l library.removeNX.fna28_union_LCA$thread )| awk -F ' ' '{sum+=$1};END{print sum/11545173}'
done