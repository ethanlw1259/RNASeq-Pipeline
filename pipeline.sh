#!/usr/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=%u@tgen.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=120gb
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/%u/logFiles/pipeline%j.log

usage() { echo "Usage: $0 [-i <FASTQ1>] [-p <FASTQ2>] [-d <refGenome>] [-r <salmonIndex>] [-a <refGTF>] [-q <workingDir>] [-h identifier] [-f = Skip FastQC] [-t = download Salmon index] [-s = download STAR index] [-g = download hg38 GTF]" 1>&2; exit 1; }

skipFastQC='false'
downloadSalmonIndex='false'
downloadStarIndex='false'
downloadGTF='false'

while getopts ":i:p:d:r:a:q:h:ftsg" o; do
	case ${o} in
		i)
			fastqFile1=${OPTARG}
			;;
		p)
			fastqFile2=${OPTARG}
			;;
		d)
			refGenome=${OPTARG}
			;;
		r)
			salmonIndex=${OPTARG}
			;;
		a)
			referenceGTF=${OPTARG}
			;;
		q)
			directory=${OPTARG}
			;;
		h)
			identifier=${OPTARG}
			;;
		f)
			skipFastQC='true'
			;;
		t)
			downloadSalmonIndex='true'
			;;
		s)
			downloadStarIndex='true'
			;;
		g)
			downloadGTF='true'
			;;
		*)
			usage
			;;
	esac
done
shift $((OPTIND-1))

if [ -z "$fastqFile1" ] || [ -z "$fastqFile2" ] ||  [ -z "$directory" ] | [ -z "$identifier" ]; then
    echo "Failed to read one or more required input arguements"
	echo "
	fastqFile1 == ${fastqFile1}
	fastqFile2 == ${fastqFile2}
	refGenome == ${refGenome}
	salmonIndex == ${salmonIndex}
	referenceGTF == ${referenceGTF}
	directory == ${directory}
	identifier == ${identifier}
	"
	usage
fi

module purge
module load singularity/3.11.5
module load FastQC/0.12.1-Container
module load petagene

cd $directory

echo "Inputs args are as follows:
fastqFile1 == ${fastqFile1}
fastqFile2 == ${fastqFile2}
refGenome == ${refGenome}
salmonIndex == ${salmonIndex}
referenceGTF == ${referenceGTF}
directory == ${directory}
identifier == ${identifier}
skipFastQC == ${skipFastQC}
downloadSalmonIndex == ${downloadSalmonIndex}
downloadStarIndex == ${downloadStarIndex}
downloadGTF == ${downloadGTF}
"

#Concatenate all fastq read 1 files into one temporarily, works whether or not files are gzipped
echo 'Starting zcat of R1 FastQ files'
zcat -f $fastqFile1 > ${identifier}.input_R1.fq
status=$?
if [ $status != 0 ];
then
	echo "Failed to zcat together all R1 FastQ files";
	exit $status
fi
echo 'R1 FastQ files combined successfully'


#Concatenate all fastq read 2 files into one temporarily, works whether or not files are gzipped
echo 'Starting zcat of R2 FastQ files'
zcat -f $fastqFile2 > ${identifier}.input_R2.fq
status=$?
if [ $status != 0 ];
then
	echo "Failed to zcat together all R2 FastQ files";
	exit $status
fi
echo 'R2 FastQ files combined successfully'

if [ "$skipFastQC" = "false" ]; then
	echo "Starting pre-trimming FastQC"
	mkdir -p $identifier/fastQC
	fastqc --extract -t 16 -o $identifier/fastQC ${identiier}.input_R1.fq ${identiier}.input_R2.fq 
	status=$?
	if [ $status != 0 ];
	then
		echo "Failed to complete pre-trimming FastQC"
		exit $status
	fi
	echo 'Pre-trimming FastQC complete'
fi


tFile1=trimmed_${identifier}.input_R1.fq
tFile2=trimmed_${identifier}.input_R2.fq

singularity exec --bind $PWD --pwd $PWD --workdir /scratch/$USER --cleanenv --contain -B /home -B /scratch /scratch/eleiter-weintraub/rnaSeqProj/scripts/cutadapt_4.8--py39hff71179_1.sif \
cutadapt -m 15 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  \
-o ${tFile1} -p ${tFile2} ${identifier}.input_R1.fq ${identifier}.input_R2.fq
status=$?
if [ $status != 0 ];
then
	echo "Cutadapt failed to complete";
	exit $status
fi
echo 'Cutadapt completed'


rm ${identifier}.input_R1.fq
rm ${identifier}.input_R2.fq


if [ "$skipFastQC" = "false" ]; then
	echo "Starting post-trimming FastQC"
	fastqc -o ${identifier}/fastQC $tFile1 $tFile2
	status=$?
	if [ $status != 0 ];
	then
		echo "Failed to complete post-trimming FastQC"
		exit $status
	fi
	echo 'Post-trimming FastQC complete'
fi

echo 'Exporting .yaml file for refgenie initiation'
rm genome_config.yaml
export REFGENIE='genome_config.yaml'

echo 'Initiating refgenie'
singularity exec  -B /home,/scratch /home/clegendre/sif_images/refgenie_0.12.1--pyhdfd78af_0.sif refgenie init -c $REFGENIE
status=$?
if [ $status != 0 ];
then
	echo "Failed to initiate refgenie";
	exit $status
fi
echo 'Refgenie initiated successfully'

if [ "$downloadGTF" = "true" ]; then
	echo 'Downloading GTF from refgenie'
	rm -r alias/hg38/ensembl_gtf
	singularity exec -B /home,/scratch /home/clegendre/sif_images/refgenie_0.12.1--pyhdfd78af_0.sif refgenie pull --force-overwrite hg38/ensembl_gtf:default
	status=$?
	if [ $status != 0 ];
	then
		echo "Failed to download GTF from refgenie";
		exit $status
	fi
	echo 'GTF downloaded from refgenie successfully'
	echo 'Ungzipping gtf file'
	zcat alias/hg38/ensembl_gtf/default/hg38.gtf.gz > hg38.gtf
	referenceGTF='hg38.gtf'
fi

if [ "$downloadSalmonIndex" = "true" ]; then
	echo 'Downloading Salmon index from refgenie'
	rm -r alias/hg38/salmon_sa_index/default
	singularity exec  -B /home,/scratch /home/clegendre/sif_images/refgenie_0.12.1--pyhdfd78af_0.sif refgenie pull --force-overwrite --pull-large hg38/salmon_sa_index:default
		status=$?
	if [ $status != 0 ];
	then
		echo "Failed to download Salmon index from refgenie";
		exit $status
	fi
	echo 'Salmon index downloaded from refgenie successfully'
	salmonIndex='alias/hg38/salmon_sa_index/default'
fi


if [ "$downloadStarIndex" = "true" ]; then
	echo 'Downloading STAR index from refgenie'
	rm -r alias/hg38/star_index/default
	singularity exec  -B /home,/scratch /home/clegendre/sif_images/refgenie_0.12.1--pyhdfd78af_0.sif refgenie pull --force-overwrite --pull-large hg38/star_index:default
	if [ $status != 0 ];
	then
		echo "STAR index failed to download from refgenie";
		exit $status
	fi
	echo 'STAR index downloaded from refgenie successfully'
	refGenome='alias/hg38/star_index/default'
	
fi


#STILL NEED TO CHANGE INPUT FILES TO REFLECT TRIMMING AND UPDATE SOME PARAMETERS#
echo 'Starting Salmon run'
singularity exec --bind $PWD --pwd $PWD --workdir /scratch/$USER --cleanenv --contain -B /home -B /scratch /home/tgenref/containers/salmon_1.10.1--hecfa306_2.sif \
	salmon quant \
	--threads 12 \
	--index $salmonIndex \
	--libType A \
	-1 $tFile1 \
	-2 $tFile2 \
	--validateMappings \
	--numBootstraps 100 \
    --seqBias \
    --gcBias \
	--writeMappings=${identifier}/${identifier}.sam \
	--output ${identifier}

if [ $status != 0 ];
then
	echo "Failed to complete Salmon run, please note trimmed and catted files will still be in working directory";
	exit $status
fi
echo 'Salmon run complete'

module load STAR/2.7.11a-container


#Runs the STAR alignment in TranscriptomeSAM quantMode on 20 threads
echo 'Starting STAR run'
STAR \
	--twopassMode Basic \
	--limitBAMsortRAM 9600000000 \
	--runMode alignReads \
	--outSAMtype BAM Unsorted \
    --outSAMmode Full \
    --outFilterType BySJout \
    --outSAMunmapped Within \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --alignMatesGapMax 1000000 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignSJoverhangMin 18 \
    --alignSJDBoverhangMin 18 \
    --chimSegmentMin 18 \
    --chimJunctionOverhangMin 18 \
    --outSJfilterOverhangMin 18 18 18 18 \
    --alignTranscriptsPerReadNmax 50000 \
    --genomeLoad NoSharedMemory \
	--runThreadN 20 \
	--genomeDir $refGenome \
	--readFilesIn $tFile1 $tFile2 \
	--quantMode GeneCounts \
	--sjdbGTFfile $referenceGTF \
	--outFileNamePrefix ${identifier}/star.${identifier}
if [ $status != 0 ];
then
	echo "Failed to complete STAR run, please note trimmed and catted files will still be in working directory";
	exit $status
fi
echo 'STAR run complete'

rm $tFile1
rm $tFile2
