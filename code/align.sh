#!/bin/bash
usage="./align.sh [-h] [-o output dir -f full file name for input bam -b base name -i index] -- program to align ICEMaP data"

if [[ "$1" =~ ^((-{1,2})([Hh]$|[Hh][Ee][Ll][Pp])|)$ ]]; then
    echo $usage; exit 1

fi
while getopts :o:f:b:i: option; do
    case "$option" in
        o)
            dir=${OPTARG}
            echo "output directory is ${outdir}"
            ;;
        f)
            file=${OPTARG}
            echo "filename is ${file}"
            ;;
        b)
            base=${OPTARG}
            echo "base name is ${base}"
            ;;
        i)
            index=${OPTARG}
            echo "index is ${index}"
            ;;
        :)
            echo "Option -${OPTARG} requires an argument."
            exit 1
            ;;
        ?)
            echo "Invalid option: -${OPTARG}."
            exit 1
            ;;
    esac
done

bedtools bamtofastq -i ${file} -fq ${dir}/${base}.1.fastq -fq2 ${dir}/${base}.2.fastq

flash2 -O -m 40 -M 260 -o ${dir}/${base} ${dir}/${base}.1.fastq ${dir}/${base}.2.fastq

bowtie2 --very-sensitive-local --no-unal -x ${index} -U ${dir}/${base}.extendedFrags.fastq -S ${dir}/${base}.mapped.raw.sam 2> ${dir}/${base}.MapStats.txt

samtools view -h -b -S -o ${dir}/${base}.mapped.raw.bam ${dir}/${base}.mapped.raw.sam; done

samtools sort -o ${dir}/${base}.mapped.sorted.bam ${dir}/${base}.mapped.raw.bam; done
