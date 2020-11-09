export OMP_NUM_THREADS=4
for VL in 2 4 8 16 32 64 128 256; do
    make -B VL=${VL} B
    for ((L=1024; L<=$((64*1024*1024)); L*=2)) ; do
	./space-time-soa $L
    done | tee soa-vl${VL}.txt
done
