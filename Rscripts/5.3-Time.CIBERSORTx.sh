## run this script with 
## ./5.3-Time.CIBERSORTx.sh 2>&1 | tee ../Result/Runtime/cibersortx/CIBERSORTx.log


#!/bin/bash
# Basic while loop
docker pull cibersortx/hires
counter=1
ncopies=10
while [ $counter -le $ncopies ]
do
echo working on copy $counter

echo start k_3.m_1000.n_500.t_$counter :  $(date -u)
docker run -v /Users/johnsonchen/Documents/Research/Unico/Unico2023/Result/Runtime/cibersortx/src:/src/data -v /Users/johnsonchen/Documents/Research/Unico/Unico2023/Result/Runtime/cibersortx/res:/src/outdir cibersortx/hires --username zchen05@g.ucla.edu --token 2dab6585dbf9c9e98b7771dd8d15f46a --mixture X.k_3.m_1000.n_500.t_$counter.txt  --cibresults W.k_3.m_1000.n_500.t_$counter.txt --sigmatrix  LM22.txt --threads 8 
echo finished $counter :  $(date -u)

echo start k_4.m_1000.n_500.t_$counter :  $(date -u)
docker run -v /Users/johnsonchen/Documents/Research/Unico/Unico2023/Result/Runtime/cibersortx/src:/src/data -v /Users/johnsonchen/Documents/Research/Unico/Unico2023/Result/Runtime/cibersortx/res:/src/outdir cibersortx/hires --username zchen05@g.ucla.edu --token 2dab6585dbf9c9e98b7771dd8d15f46a --mixture X.k_4.m_1000.n_500.t_$counter.txt  --cibresults W.k_4.m_1000.n_500.t_$counter.txt --sigmatrix  LM22.txt --threads 8 
echo finished $counter :  $(date -u)

echo start k_5.m_1000.n_500.t_$counter :  $(date -u)
docker run -v /Users/johnsonchen/Documents/Research/Unico/Unico2023/Result/Runtime/cibersortx/src:/src/data -v /Users/johnsonchen/Documents/Research/Unico/Unico2023/Result/Runtime/cibersortx/res:/src/outdir cibersortx/hires --username zchen05@g.ucla.edu --token 2dab6585dbf9c9e98b7771dd8d15f46a --mixture X.k_5.m_1000.n_500.t_$counter.txt  --cibresults W.k_5.m_1000.n_500.t_$counter.txt --sigmatrix  LM22.txt --threads 8 
echo finished $counter :  $(date -u)

echo start k_6.m_1000.n_500.t_$counter :  $(date -u)
docker run -v /Users/johnsonchen/Documents/Research/Unico/Unico2023/Result/Runtime/cibersortx/src:/src/data -v /Users/johnsonchen/Documents/Research/Unico/Unico2023/Result/Runtime/cibersortx/res:/src/outdir cibersortx/hires --username zchen05@g.ucla.edu --token 2dab6585dbf9c9e98b7771dd8d15f46a --mixture X.k_6.m_1000.n_500.t_$counter.txt  --cibresults W.k_6.m_1000.n_500.t_$counter.txt --sigmatrix  LM22.txt --threads 8 
echo finished $counter :  $(date -u)


echo start k_6.m_1000.n_100.t_$counter :  $(date -u)
docker run -v /Users/johnsonchen/Documents/Research/Unico/Unico2023/Result/Runtime/cibersortx/src:/src/data -v /Users/johnsonchen/Documents/Research/Unico/Unico2023/Result/Runtime/cibersortx/res:/src/outdir cibersortx/hires --username zchen05@g.ucla.edu --token 2dab6585dbf9c9e98b7771dd8d15f46a --mixture X.k_6.m_1000.n_100.t_$counter.txt  --cibresults W.k_6.m_1000.n_100.t_$counter.txt --sigmatrix  LM22.txt --threads 8 
echo finished $counter :  $(date -u)

echo start k_6.m_1000.n_250.t_$counter :  $(date -u)
docker run -v /Users/johnsonchen/Documents/Research/Unico/Unico2023/Result/Runtime/cibersortx/src:/src/data -v /Users/johnsonchen/Documents/Research/Unico/Unico2023/Result/Runtime/cibersortx/res:/src/outdir cibersortx/hires --username zchen05@g.ucla.edu --token 2dab6585dbf9c9e98b7771dd8d15f46a --mixture X.k_6.m_1000.n_250.t_$counter.txt  --cibresults W.k_6.m_1000.n_250.t_$counter.txt --sigmatrix  LM22.txt --threads 8 
echo finished $counter :  $(date -u)

echo finished on copy $counter
echo $(date -u)
((counter++))
done
echo All done


