#!/bin/sh

#path = ~/linxing2/workspace/DMY_sci_data/results/DEGs

#---------------D10 RNA-Seq DEMs extract----------------------------#

cat 10day_XXXYD_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 >=1' |awk '$7 <= 0.05' |sed 's/"//g' |awk '{print $1}' >10d_XXXYD_orth_XYD_high_degs_list1.csv
cat 10day_XXXYD_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 <=-1' |awk '$7 <= 0.05' |grep -v padj|sed 's/"//g' |awk '{print $1}' >10d_XXXYD_orth_XX_high_degs_list1.csv
cat 10day_XXXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 >=1' |awk '$7 <= 0.05' |sed 's/"//g' |awk '{print $1}' >10d_XXXY_orth_XY_high_degs_list1.csv
cat 10day_XXXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 <=-1' |awk '$7 <= 0.05' |grep -v padj|sed 's/"//g' |awk '{print $1}' >10d_XXXY_orth_XX_high_degs_list1.csv
cat 10day_XYDXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 >=1' |awk '$7 <= 0.05' |sed 's/"//g' |awk '{print $1}' >10d_XYDXY_orth_XYD_high_degs_list1.csv
cat 10day_XYDXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 <=-1' |awk '$7 <= 0.05' |grep -v padj|sed 's/"//g' |awk '{print $1}' >10d_XYDXY_orth_XY_high_degs_list1.csv

#---------------D30 RNA-Seq DEMs extract----------------------------#

cat 30day_XXXYD_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 >=1' |awk '$7 <= 0.05' |sed 's/"//g' |awk '{print $1}' >30d_XXXYD_orth_XYD_high_degs_list1.csv
cat 30day_XXXYD_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 <=-1' |awk '$7 <= 0.05' |grep -v padj|sed 's/"//g' |awk '{print $1}' >30d_XXXYD_orth_XX_high_degs_list1.csv
cat 30day_XXXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 >=1' |awk '$7 <= 0.05' |sed 's/"//g' |awk '{print $1}' >30d_XXXY_orth_XY_high_degs_list1.csv
cat 30day_XXXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 <=-1' |awk '$7 <= 0.05' |grep -v padj|sed 's/"//g' |awk '{print $1}' >30d_XXXY_orth_XX_high_degs_list1.csv
cat 30day_XYDXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 >=1' |awk '$7 <= 0.05' |sed 's/"//g' |awk '{print $1}' >30d_XYDXY_orth_XYD_high_degs_list1.csv
cat 30day_XYDXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 <=-1' |awk '$7 <= 0.05' |grep -v padj|sed 's/"//g' |awk '{print $1}' >30d_XYDXY_orth_XY_high_degs_list1.csv

#---------------D120 RNA-Seq DEMs extract----------------------------#

cat 120day_XXXYD_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 >=1' |awk '$7 <= 0.05' |sed 's/"//g' |awk '{print $1}' >120d_XXXYD_orth_XYD_high_degs_list1.csv
cat 120day_XXXYD_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 <=-1' |awk '$7 <= 0.05' |grep -v padj|sed 's/"//g' |awk '{print $1}' >120d_XXXYD_orth_XX_high_degs_list1.csv
cat 120day_XXXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 >=1' |awk '$7 <= 0.05' |sed 's/"//g' |awk '{print $1}' >120d_XXXY_orth_XY_high_degs_list1.csv
cat 120day_XXXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 <=-1' |awk '$7 <= 0.05' |grep -v padj|sed 's/"//g' |awk '{print $1}' >120d_XXXY_orth_XX_high_degs_list1.csv
cat 120day_XYDXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 >=1' |awk '$7 <= 0.05' |sed 's/"//g' |awk '{print $1}' >120d_XYDXY_orth_XYD_high_degs_list1.csv
cat 120day_XYDXY_deseq2_res.csv |sed 's/,/\t/g' |awk '$3 <=-1' |awk '$7 <= 0.05' |grep -v padj|sed 's/"//g' |awk '{print $1}' >120d_XYDXY_orth_XY_high_degs_list1.csv

