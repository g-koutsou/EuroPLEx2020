export OMP_NUM_THREADS=4
for ((L=128; L<=$((4*1024)); L*=2)) ; do
	echo $(./lapl $L | grep lapl\() $L
done | tee lapl.txt
