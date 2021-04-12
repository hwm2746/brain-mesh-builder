## merge all w_dist*.txt files to view 

cd distance/
echo 'mol new' > w_dist_all.txt
for j in {16..49}; do                                                         
    cat w_dist$j.txt >> w_dist_all.txt
done

cd .. 
