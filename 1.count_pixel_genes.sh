export DIR=/storage/peiweikeLab/guochenyu/Liuchang/project2/02.STtools/liver/02.all_liver/01.binsize300/allSolo.out/GeneFull/ordered
python 1.count_pixel_UMIs_genes.py $DIR/matrix.mtx.gz $DIR/spatialcoordinates.sorted.txt.gz $DIR/barcodes.tsv.gz $DIR/features.tsv.gz UMI_count gene_count 1>notincorrdi
