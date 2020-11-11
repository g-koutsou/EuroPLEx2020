export OMP_NUM_THREADS=4
for ((VL=4; VL<=64; VL*=2)); do
    make clean && VL=$VL make B
    for ((L=128; L<=$((4*1024)); L*=2)) ; do
	echo $(./laplv $L | grep lapl\() $L
    done | tee laplv-VL${VL}.txt
done
