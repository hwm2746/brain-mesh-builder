## merge all w_dist*.txt files to view in vmd
# In vmd console, type: source w_dist_all.txt

cd distance/
echo 'mol new' > w_dist_all.txt
for j in {18..68}; do                                                         
    cat w_dist$j.txt >> w_dist_all.txt
done                                                                            
cd ../
            
