#!/bin/bash

#SBATCH --mail-user sonia.sarmientocabello@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition cpu
#SBATCH --array 1-2
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --job-name=Switch
#SBATCH -o %x_%A_%a.stdout
#SBATCH -e %x_%A_%a.stderr
#SBATCH --account jgoudet_barn_owl

# Load modules
module load gcc shapeit4 bcftools htslib/1.2

# Set your path
myPATH=/Users/soniasarmiento/Reference_panel

# Set path to Switch code (available in github)
Switch=${myPATH}/SHAPEIT5_switch_static

#Run in array for each scaffold (see above array option). 
#Here, we decide to run for the two biggest scaffolds in barn owl: Super-Scaffold_14 and Super-Scaffold_6
runLINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${myPATH}/SwitchOutput/BigScaffolds.txt)

SwitchComparison(){
	CHR=$1
	#Validation for code is a VCF file with all unphased parents of the offspring chosen (as well as the offspring themselves)
	VAL=${myPATH}/Unphased_Parents_and_Offspring

	#PED File with pedigree information of parents and offspring
	ped=${myPATH}/PEDfile_SwitchInput.list

	#Phasing output to compare with Mendelian inheritance from the Validation VCF file
	EST=${myPATH}/RefPanelShapeIt_OffspringONLY

	#Define output directory and file prefix
	OUTPUT_PREFIX=${myPATH}/SwitchOutput

	${Switch} --validation ${VAL}_${CHR}.vcf.gz --estimation ${EST}_${CHR}.vcf.gz --region ${CHR} --pedigree ${ped} --output ${OUTPUT_PREFIX}_${CHR}
}

#Export the function
export -f SwitchComparison

#Launch the function
SwitchComparison ${runLINE}
