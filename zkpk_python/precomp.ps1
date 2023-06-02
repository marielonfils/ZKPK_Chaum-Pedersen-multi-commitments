for (($l=1); $l -lt 5; $l *=2){
    for (($k=1); $k -lt 5; $k *=2)
    {
        python logarithmic_precomputation.py $k $l >> .\precomp_test3\pre_log_$k`_$l.txt
    }
}
for (($l=8); $l -lt 17; $l *=2){
    for (($k=1); $k -lt 5; $k *=2)
    {
        python logarithmic_precomputation.py $k $l >> .\precomp_test3\pre_log_$k`_$l.txt
    }
}
for (($l=1); $l -lt 5; $l *=2){
    for (($k=1); $k -lt 5; $k *=2)
    {
        python logarithmic_precomputation.py $k $l >> .\precomp_test3\pre_log_$k`_$l.txt
    }
}
for (($l=1); $l -lt 5; $l *=2){
    for (($k=6); $k -lt 11; $k +=2)
    {
        python logarithmic_precomputation.py $k $l >> .\precomp_test3\pre_log_$k`_$l.txt
    }
}
for (($l=8); $l -lt 17; $l *=2){
    for (($k=6); $k -lt 11; $k +=2)
    {
        python logarithmic_precomputation.py $k $l >> .\precomp_test3\pre_log_$k`_$l.txt
    }
}
for (($l=32); $l -lt 129; $l *=2){
    for (($k=1); $k -lt 5; $k *=2)
    {
        python logarithmic_precomputation.py $k $l >> .\precomp_test3\pre_log_$k`_$l.txt
    }
}
for (($l=32); $l -lt 129; $l *=2){
    for (($k=6); $k -lt 11; $k +=2)
    {
        python logarithmic_precomputation.py $k $l >> .\precomp_test3\pre_log_$k`_$l.txt
    }
}