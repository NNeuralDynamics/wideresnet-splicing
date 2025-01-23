#$ -q gpu
#$ -cwd
#$ -N test_spliceai
#$ -e Logs/
#$ -o Logs/
#$ -l gpus=1
#$ -l h_vmem=500g

python -u test_model.py $1 $2 > Outputs/TEST_SpliceAI${1}_${2}_.txt
