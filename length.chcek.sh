#readlength.check in epsiolon 
for j in breast.GSE143956  breast.GSE167956 breast.GSE114937  breast.GSE162069  breast.GSE137752  breast.GSE157574  breast.GSE150099  breast.GSE162438  breast.GSE168706  breast.GSE173616 
do
	echo ${j}
	echo "################  Read length check ################"
	cd ${j}/fastq/
	ls |grep fastq |sed -e 's/.fastq//g' > ./f.list
	
	for i in $(cat ./f.list)
	do
		echo ${i}".fastq"
		cat ${i}.fastq |grep length |awk -F [=] '{print $2}' |sort |uniq -c > ./${i}_readL.txt
	done
	echo "################  Read Length check finished ################"

	cd ../../

done
