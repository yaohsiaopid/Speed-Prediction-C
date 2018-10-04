<<<<<<< HEAD
for i in `seq 1 100;`
do
    echo "-----------$i-------------"
=======
for i in `seq 1 10;`
do
    echo $i
>>>>>>> 69e31b155a183946695a0171f099b860e4c6564c
    diff "./matlabout/t$i.txt" "./couts/t$i.txt"
    echo "------------"
done