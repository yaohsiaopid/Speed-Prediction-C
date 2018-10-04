for i in `seq 1 10;`
do
    echo $i
    diff "./matlabout/t$i.txt" "./couts/t$i.txt"
    echo "------------"
done