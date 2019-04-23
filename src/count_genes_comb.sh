cd exprs_go_comb2
echo -e "go_id\tsize" > ../output/go2size.go_comb2.txt
find . -name "*.txt" -print | xargs wc -l | awk '{print $2,($1-1)}' | perl -pe 's/.\/exprs.log10.E1.//g' | perl -pe 's/.txt//g' | perl -pe 's/ /\t/g' | grep '^GO' >> ../output/go2size.go_comb2.txt
cd ../

cd exprs_go_comb3
echo -e "go_id\tsize" > ../output/go2size.go_comb3.txt
find . -name "*.txt" -print | xargs wc -l | awk '{print $2,($1-1)}' | perl -pe 's/.\/exprs.log10.E1.//g' | perl -pe 's/.txt//g' | perl -pe 's/ /\t/g' | grep '^GO' >> ../output/go2size.go_comb3.txt
cd ../

cd exprs_go_comb2
echo -e "go_id\tsize" > ../output/go2size.go_comb2.txt
find . -name "*.txt" -print | xargs wc -l | awk '{print $2,($1-1)}' | perl -pe 's/.\/exprs.log10.E1.//g' | perl -pe 's/.txt//g' | perl -pe 's/ /\t/g' | grep '^GO' >> ../output/go2size.go_comb2.txt
cd ../

cd exprs_go_comb4
echo -e "go_id\tsize" > ../output/go2size.go_comb4.txt
find . -name "*.txt" -print | xargs wc -l | awk '{print $2,($1-1)}' | perl -pe 's/.\/exprs.log10.E1.//g' | perl -pe 's/.txt//g' | perl -pe 's/ /\t/g' | grep '^GO' >> ../output/go2size.go_comb4.txt
cd ../

cd exprs_go_comb5
echo -e "go_id\tsize" > ../output/go2size.go_comb5.txt
find . -name "*.txt" -print | xargs wc -l | awk '{print $2,($1-1)}' | perl -pe 's/.\/exprs.log10.E1.//g' | perl -pe 's/.txt//g' | perl -pe 's/ /\t/g' | grep '^GO' >> ../output/go2size.go_comb5.txt
cd ../

cd exprs_go_comb6
echo -e "go_id\tsize" > ../output/go2size.go_comb6.txt
find . -name "*.txt" -print | xargs wc -l | awk '{print $2,($1-1)}' | perl -pe 's/.\/exprs.log10.E1.//g' | perl -pe 's/.txt//g' | perl -pe 's/ /\t/g' | grep '^GO' >> ../output/go2size.go_comb6.txt
cd ../

cd exprs_go_comb5_del1
echo -e "go_id\tsize" > ../output/go2size.go_comb5_del1.txt
find . -name "*.txt" -print | xargs wc -l | awk '{print $2,($1-1)}' | perl -pe 's/.\/exprs.log10.E1.//g' | perl -pe 's/.txt//g' | perl -pe 's/ /\t/g' | grep '^comb_del-' >> ../output/go2size.go_comb5_del1.txt
cd ../

cd exprs_go_comb5_del2
echo -e "go_id\tsize" > ../output/go2size.go_comb5_del2.txt
find . -name "*.txt" -print | xargs wc -l | awk '{print $2,($1-1)}' | perl -pe 's/.\/exprs.log10.E1.//g' | perl -pe 's/.txt//g' | perl -pe 's/ /\t/g' | grep '^comb_del' >> ../output/go2size.go_comb5_del2.txt
cd ../

cd exprs_go_comb5_del3
echo -e "go_id\tsize" > ../output/go2size.go_comb5_del3.txt
find . -name "*.txt" -print | xargs wc -l | awk '{print $2,($1-1)}' | perl -pe 's/.\/exprs.log10.E1.//g' | perl -pe 's/.txt//g' | perl -pe 's/ /\t/g' | grep '^comb_del' >> ../output/go2size.go_comb5_del3.txt
cd ../
