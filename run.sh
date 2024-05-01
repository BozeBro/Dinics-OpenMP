for i in 1 2 4 8; do 
    $1 -f testInput/$2.txt -m 1 -n $i
done
