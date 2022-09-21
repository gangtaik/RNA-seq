#!/bin/bash
#SRR donwload to Base counting
#author G.T Lee
#executing server=MRC4(zeta) mapping >> MRC2(epsilon 가능)
#breast.{geo} and SraRunTable.txt or srr.list should be firstly made before execution
#bash Rseq.sh ${geo}

geo=${1}

dir=/data/Epsilon_data12/gangtl95/BRCA/
tool=/data/Epsilon_data7/gangtl95/BRCA/

#mkdir ./breast.${geo}

sra=${dir}breast.${geo}/sra
fastq=${dir}breast.${geo}/fastq
map=${dir}breast.${geo}/hisat2.GRCh38
cnt=${dir}breast.${geo}/htseq
results=${dir}breast.${geo}/results
gtf=${tool}Homo_sapiens.GRCh38.105.gtf
aligner=${tool}hisat2-2.2.1
genome=${tool}hg38_GRCh38


echo "breast."${geo}

#mkdir -p ${sra}
#mkdir -p ${fastq}
#mkdir -p ${map}
#mkdir -p ${cnt}
#mkdir -p ${results}
#
#echo "################ Downloading SRA files ################"
#
#Sratable=SraRunTable.txt
#
#if [ -f "${Sratable}" ]
#then
#       	echo "SraRunTable.txt exist"
#
#       	more SraRunTable.txt | awk -F "," '{OFS="\t"}{print $1}' |sed -e '1d' > srr.list	
#
#	cd ${sra}
#       	
#	for i in $(cat ../srr.list)
#       	do
#               	echo ${i}
#               	
#		${tool}sratoolkit.3.0.0-ubuntu64/bin/prefetch ${i} -X 100152919112
#
#		cd ${i}
#
#		mv ${i}.sra ../
#
#		cd ../
#
#		rm -r ${i}
#       	done
#
#       	#Check Run data download properly
#       	
#	mkdir -p ${sra}/log
#       	
#	echo "#Check Run data download properly"
#       	
#	ls -al |awk -F " " '{print $9}'| grep sra > ./log/dw.srr.list
#       	
#	for i in $(cat ./log/dw.srr.list)
#       	do
#       	        echo ${i}
#
#
#               	
#		${tool}sratoolkit.3.0.0-ubuntu64/bin/vdb-validate ${i} 2> ./log/${i}.validation_out
#               	
#		if grep -q 'err' ./log/${i}.validation_out;
#               	
#		then
#                       	
#			echo 'Errors found in '${i} > error.list
#               	else
#       	                
#			echo 'No Errors found in '${i}
#
#               	fi
#       	done
#else
#       	echo "SraRunTable.txt does not exist"
#	echo "make SRR.LIST manually"
#       	
#	cd ${sra}
#       	
#	for i in $(cat ../srr.list)
#       	
#	do
#       	        
#		echo ${i}
#
#               	${tool}sratoolkit.3.0.0-ubuntu64/bin/prefetch ${i} -X 100152919112
#
#
#		cd ${i}		
#		
#		mv ${i}.sra ../
#		
#		cd ../
#
#		rm -r ${i}/
#       	done
#
#       	#Check Run data download properly
#       	
#	mkdir -p ${sra}/log
#       	
#	echo "#Check Run data download properly"
#       	
#	ls -al |awk -F " " '{print $9}'| grep SRR > ./log/dw.srr.list
#       	
#	for i in $(cat ./log/dw.srr.list)
#       	do
#               	echo ${i}
#
#               	${tool}sratoolkit.3.0.0-ubuntu64/bin/vdb-validate ${i} 2> ./log/${i}.validation_out
#               	
#		if grep -q 'err' ./log/${i}.validation_out;
#               	
#		then
#                       	echo 'Errors found in '${i} > error.list
#               	
#		else
#                       	echo 'No Errors found in '${i}
#               	
#		fi
#
#       	done
#
#fi
#echo "################ Downloading Finished ################"
#
#
#echo "################ Transform SRA files to Fastq files ################"
#
#cd ${fastq}
#
#for j in $(cat ${sra}/log/dw.srr.list)
#do
#	
#	echo ${j}
#	
#	${tool}sratoolkit.3.0.0-ubuntu64/bin/fasterq-dump ../sra/${j}
#
#done
#
#cd ${sra}
#
#ls -lh |grep sra |awk -F " " '{OFS="\t"}{print $9,$5}' > ./size.txt
#
#rm *sra
#
#echo "################ Transforming Finished ################"
#
#
#echo "################ Combining Same Sample FASTQ files ################"
#
#cd ${fastq}
#
#ls |grep fastq |awk -F '_' '{print $1}' |uniq -c |awk -F " " '{OFS="\t"}{print $2,$1}' > ./cnt.txt
#
#cp ${dir}combine.py ${dir}breast.${geo}/
#		
#python3 ../combine.py -i ${fastq} -o ${fastq} -T ../SraRunTable.txt -C ./cnt.txt
#
#
#more ../SrrGsm.txt |sort -k 2 |uniq -f1 -D |awk -F "\t" '{print "rm ",$1"*"}' > rm_uncombined_srr.sh
#
#bash rm_uncombined_srr.sh
#
#echo "################ Combining Finished ################"
#
#
#
#
#echo "################  Base Counting ################"
#
#cd ${fastq}
#
#ls |grep fastq |awk -F '_' '{print $1}' |uniq > tmp.list
#
#cat tmp.list |sed -e 's/.fastq//g' > srr.list
#
#rm tmp.list
#
#for sname in $(cat ./srr.list)
#do
#	if [ $(ls | grep ${sname} | grep fastq | wc -l) == 2 ]
#	then
#		echo "Library layout is PAIRED"
#
#		echo -e ${sname}_1.fastq >> base.p.1.txt 2>&1 
#
#		awk 'BEGIN{sum=0;}{if(NR%4==2){sum+=length($0);}}END{print sum;}' ${sname}_1.fastq >> base.p.2.txt 2>&1
#	
#		echo -e ${sname}_2.fastq >> base.p.3.txt 2>&1
#
#		awk 'BEGIN{sum=0;}{if(NR%4==2){sum+=length($0);}}END{print sum;}' ${sname}_2.fastq >> base.p.4.txt 2>&1
#
#
#	else
#		echo "Library layouts is Single"
#
#		echo -e ${sname}.fastq >> base.s.1.txt 2>&1
#
#		awk 'BEGIN{sum=0;}{if(NR%4==2){sum+=length($0);}}END{print sum;}' ${sname}.fastq >> base.s.2.txt  2>&1
#
#
#	fi
#
#
#done
#
#if [ $(cat base.p.1.txt | wc -l) != 0 ]
#then
#	paste -d "\t" base.p.1.txt base.p.2.txt base.p.3.txt base.p.4.txt > base.p.5.txt
#	cat base.p.5.txt |awk -F "\t" '{OFS="\t"}{print $1, $2+$4}' > base.paired.txt
#
#
#	rm base.p.1.txt
#	rm base.p.2.txt
#	rm base.p.3.txt
#	rm base.p.4.txt
#	rm base.p.5.txt
#
#else
#	
#	echo "No Paired Data"
#
#fi
#
#if [ $(cat base.s.1.txt | wc -l) != 0 ]
#then
#	
#	paste -d "\t" base.s.1.txt base.s.2.txt > base.single.txt
#	
#	rm base.s.1.txt
#	
#	rm base.s.2.txt
#
#
#else
#	
#	echo "No Single Data"
#
#
#fi
#
#echo "################  Base Counting Finished ################"
#
#chmod 777 ${dir}breast.${geo}
#chmod 777 ${fastq}
#chmod 777 ${fastq}/*fastq
#
#echo "################  Alignment  ################"
#
#logf=${map}/log.txt
#
#
#cd ${fastq}
#
#for sname in $(cat ${fastq}/srr.list)
#do
#	echo -e '\n'${sname} >> ${logf} 2>&1
#	if [ $(ls| grep ${sname} |wc -l) == 2 ]
#	then
#		echo ${sname} is paired-end
#		sname_1=${sname}_1.fastq
#		sname_2=${sname}_2.fastq
#
#		${aligner}/hisat2 \
#                        -p 4 \
#                        -x ${genome}/genome -1 ${fastq}/${sname_1} \
#                        -2 ${fastq}/${sname_2} \
#                        -S ${map}/${sname}.sam >> ${logf} 2>&1
#                echo ${map}/${sname}.bam
#                samtools view -bS ${map}/${sname}.sam > ${map}/${sname}.bam
#                samtools view -F 268 -bS ${map}/${sname}.bam > ${map}/${sname}.uniqued.bam
#                samtools sort ${map}/${sname}.uniqued.bam -o ${map}/${sname}.sorted.bam
#                samtools index ${map}/${sname}.sorted.bam
#
#	else
#                echo ${sname} is single-end
#                sname_sgl=${sname}.fastq
#                
#                ${aligner}/hisat2 \
#                        -p 4 \
#                        -x ${genome}/genome \
#                        -U ${fastq}/${sname_sgl} \
#                        -S ${map}/${sname}.sam >> ${logf} 2>&1
#                samtools view -bS ${map}/${sname}.sam > ${map}/${sname}.bam
#                samtools view -F 268 -bS ${map}/${sname}.bam > ${map}/${sname}.uniqued.bam
#                samtools sort ${map}/${sname}.uniqued.bam -o ${map}/${sname}.sorted.bam
#                samtools index ${map}/${sname}.sorted.bam
#        fi
#
#	rm ${map}/${sname}.sam
#	rm ${map}/${sname}.uniqued.bam
#
#done
#
#echo "################  Alignment  Finished ################"
#
#
#
#echo "################  Flag check  ################"
#
#cd ${map}
#mkdir ./strandness
#
##positive strand feature gene NUMA1, ACTB (+ strand)
#
#ci_NUMA1=11:72013929-72015217
#fn_NUMA1=sorted.NUMA1.txt
#
#ci_ACTB=7:5528303-5528711
#fn_ACTB=sorted.ACTB.txt
#
##negative strand feature gene GAPDH (- strand)
#ci_GAPDH=12:6537601-6537982
#fn_GAPDH=sorted.GAPDH.txt
#
#
##Flag check from sorted.bam
#
#ls |grep bai |sed -e 's/.sorted.bam.bai//g' > sname.list
#
#for sname in $(cat ./sname.list)
#do
#
#        echo ${sname}
#
#        samtools view ${sname}.sorted.bam ${ci_NUMA1}|\
#                cut -f 2 |\
#                sort |\
#                uniq -c |\
#                sort -rn |\
#                sed -e 's/^ *//g' |\
#                sed -e 's/ /\t/g' > ./strandness/${sname}.${fn_NUMA1}
#        samtools view ${sname}.sorted.bam ${ci_ACTB}|\
#                cut -f 2 |\
#                sort |\
#                uniq -c |\
#                sort -rn |\
#                sed -e 's/^ *//g'|\
#                sed -e 's/ /\t/g' > ./strandness/${sname}.${fn_ACTB}
#        samtools view ${sname}.sorted.bam ${ci_GAPDH}|\
#                cut -f 2 |\
#                sort |\
#                uniq -c |\
#                sort -rn |\
#                sed -e 's/^ *//g' |\
#                sed -e 's/ /\t/g' > ./strandness/${sname}.${fn_GAPDH}
#
#        samtools view ${sname}.sorted.bam |\
#                cut -f 2 |\
#                sort |\
#                uniq -c |\
#                sort -rn |\
#                sed -e 's/^ *//g' |\
#                sed -e 's/ /\t/g' > ./strandness/${sname}_flag.txt
#
#
#
#        samtools view ${sname}.sorted.bam |\
#                cut -f 3 |\
#                sort |\
#                uniq -c |\
#                sort -rn |\
#                sed -e 's/^ *//g' |\
#                sed -e 's/ /\t/g' > ./strandness/${sname}_pos.txt
#
#	
#	
#done
#
##First=BRCA.dir
##Second=geo
#
#Rscript ${dir}Flag_check.R ${dir} ${geo}
#
#echo "################  Flag check Finished  ################"


echo "################  Read Counting  ################"

cd ${map}
ls |grep bai |sed -e 's/.sorted.bam.bai//g' > sortbam.list
cd ${cnt}

readcount () {
  echo "Input 1 if strand option is yes"
  echo "Input 2 if strand option is no"
  echo "Input 3 if strand option is reverse"
  echo "Input x to exit the script"
  read -n 1 -p "Input Strand Status:" strandness
  if [ "$strandness" == 1 ]; then
	  echo -e  "\nStrand exists"

	  for sname in $(cat ${map}/sortbam.list)
	  do
		  echo ${sname}.sorted.bam
                  file_samout=${sname}.out.sam
                  echo ${file_samout}
                  file_in=${map}/${sname}.sorted.bam
                  echo ${file_in}
                  file_out=${sname}.htseq.txt
                  echo ${file_out}
                  /home/tools/anaconda3/bin/htseq-count --format=bam \
			  --nprocesses=4 \
                          --order=pos \
                          --stranded=yes \
                          --type=exon \
                          --mode=union \
                          --samout=${file_samout} ${file_in} ${gtf} > ${cnt}/${file_out}
	  done
            
        elif [ "$strandness" == 2 ]; then
            echo -e "\nNo Strand"
	    for sname in $(cat ${map}/sortbam.list)
	    do
		    echo ${sname}.sorted.bam
		    file_samout=${sname}.out.sam		    
		    echo ${file_samout}
		    file_in=${map}/${sname}.sorted.bam                  
		    echo ${file_in}                 
		    file_out=${sname}.htseq.txt
		    echo ${file_out}                  
		    /home/tools/anaconda3/bin/htseq-count --format=bam \
			    --nprocesses=4 \
			    --order=pos \
			    --stranded=no \
			    --type=exon \
			    --mode=union \
			    --samout=${file_samout} ${file_in} ${gtf} > ${cnt}/${file_out}

	    done


        elif [ "$strandness" == 3 ]; then
            echo -e "\nStrand Reverse"

	    for sname in $(cat ${map}/sortbam.list)
	    do
		    echo ${sname}.sorted.bam
		    file_samout=${sname}.out.sam
		    echo ${file_samout}
		    file_in=${map}/${sname}.sorted.bam
		    echo ${file_in}
		    file_out=${sname}.htseq.txt
		    echo ${file_out}
		    /home/tools/anaconda3/bin/htseq-count --format=bam \
			    --nprocesses=4 \
			    --order=pos \
			    --stranded=reverse \
			    --type=exon \
			    --mode=union \
                            --samout=${file_samout} ${file_in} ${gtf} > ${cnt}/${file_out}



	    done

        elif [ "$strandness" == "x" ];then
            quitprogram

        elif [ "$strandness" == "X" ];then
            quitprogram
        else
            echo "You have entered an invallid selection!"
            echo "Please try again!"
            echo ""
            echo "Press any key to continue..."
            read -n 1
            clear
            readcount
        fi
}

readcount

echo "################  Read Counting Finished  ################"


echo "################  Merging  ################"

Rscript ${dir}Merge.R ${dir} ${geo}
cd ${cnt}
chmod 777 *_readcounts.txt

echo "################  Merging Finished  ################"

echo "################  Read length check ################"

cd ${fastq}
ls |grep fastq |sed -e 's/.fastq//g' > ${fastq}/f.list
for i in $(cat ${fastq}/f.list)
do
        echo ${i}".fastq"
        cat ${i}.fastq |grep length |awk -F [=] '{print $2}' |sort |uniq -c > ${fastq}/${i}_readL.txt

done

echo "################  Read Length check finished ################"

echo "breast."${geo}
cd ${dir}


