#SBATCH -o slurm-%x-%A_%2a.out # Template for the std output of the job uses the job name, the job id and the array id
#SBATCH -e slurm-%x-%A_%2a.err # Template for the std error of the job
#SBATCH --nodes 1 # We always use 1 node
#SBATCH --ntasks 1 # In this script everything is sequencial
#SBATCH --mem 50G 
#SBATCH --cpus-per-task 24 # This allows to speed the indexing
#SBATCH --time 23:00:00 # This depends on the size of the fasta
#SBATCH --job-name runMicroC # Job name that appear in squeue as well as in output and error text files


sample=mc_120h_R1_unbio
pairtools sort --nproc ${SLURM_CPUS_PER_TASK} --tmpdir=$(mktemp -d) --output-stats "${sample}.parsed.pairsam" > "${sample}.sorted.pairsam"
