export OMP_NUM_THREADS=4
for ((L=1024; L<=$((64*1024*1024)); L*=2)) ; do
    ./space-time-soa $L
done | tee soa.txt
