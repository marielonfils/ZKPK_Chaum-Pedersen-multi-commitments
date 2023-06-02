#! /bin/bash
for l in 1 2 4
do
    for k in 1 2 4
    do
        python logarithmic_precomputation.py $k $l >> ./precomp_test2/pre_log_$k\_$l.txt

    done
done

for l in 8 16
do
    for k in 1 2 4 
    do
        python logarithmic_precomputation.py $k $l >> ./precomp_test2/pre_log_$k\_$l.txt

    done
done

for l in 1 2 4
do
    for k in 6 8 10
    do
        python logarithmic_precomputation.py $k $l >> ./precomp_test2/pre_log_$k\_$l.txt

    done
done

for l in 8 16 
do
    for k in 6 8 10
    do
        python logarithmic_precomputation.py $k $l >> ./precomp_test2/pre_log_$k\_$l.txt

    done
done


for l in 32 64 128
do
    for k in 1 2 4 
    do
        python logarithmic_precomputation.py $k $l >> ./precomp_test2/pre_log_$k\_$l.txt

    done
done

#for l in 32 64 128
#do
#    for k in 6 8 10
#    do
#        python logarithmic_precomputation.py $k $l >> ./precomp_test2/pre_log_$k\_$l.txt
#
#    done
#done