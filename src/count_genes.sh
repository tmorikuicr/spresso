cd exprs_go
echo -e "go_id\tsize" > ../output/go2size.go.txt
find . -name "*.txt" -print | xargs wc -l | awk '{print $2,($1-1)}' | perl -pe 's/.\/exprs.log10.E1.//g' | perl -pe 's/.txt//g' | perl -pe 's/ /\t/g' | grep '^GO' >> ../output/go2size.go.txt
cd ../
cd exprs_peng
echo -e "go_id\tsize" > ../output/go2size.peng.txt
find . -name "*.txt" -print | xargs wc -l | awk '{print $2,($1-1)}' | perl -pe 's/.\/exprs.log10.E1.//g' | perl -pe 's/.txt//g' | perl -pe 's/ /\t/g' | grep '^[a-z]' | grep -v '^total' >> ../output/go2size.peng.txt
cd ../

