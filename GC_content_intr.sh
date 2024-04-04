for file in Hf02 Hf03 Hf04
do
	awk '{ if ($20 <= 0.01 && $21 == "-") {print $1, $2, $3, $4, $5, $6, $9, $19, $20, $21}}' /spectrum/GSCT/veronikap/DEEP/fastq/Tcells/BlEM/${file}/IRFinder-IR-dir.txt > /spectrum/GSCT/veronikap/DEEP/GCcount/BlEM/BlEM_${file}_nonretint.txt
	awk '{print $1,$2,$3,$4,$5,$6}' /spectrum/GSCT/veronikap/DEEP/GCcount/BlEM/BlEM_${file}_nonretint.txt | sed 's/  */\t/g' > /spectrum/GSCT/veronikap/DEEP/GCcount/BlEM/BlEM_${file}_nonretint.bed
	bedtools nuc -fi /spectrum/GSCT/REF/Human-GRCh38-release86/genome.fa -bed /spectrum/GSCT/veronikap/DEEP/GCcount/BlEM/BlEM_${file}_nonretint.bed -s > /spectrum/GSCT/veronikap/DEEP/GCcount/BlEM/GC_BlEM_${file}_nonretint.txt
done

for file in Hf02 Hf03 Hf04
do
	awk '{ if ($20 >= 0.1 && $21 == "-") {print $1, $2, $3, $4, $5, $6, $9, $19, $20, $21}}' /spectrum/GSCT/veronikap/DEEP/fastq/Tcells/BlEM/${file}/IRFinder-IR-dir.txt > /spectrum/GSCT/veronikap/DEEP/GCcount/BlEM/BlEM_${file}_retint.txt
	awk '{print $1,$2,$3,$4,$5,$6}' /spectrum/GSCT/veronikap/DEEP/GCcount/BlEM/BlEM_${file}_retint.txt | sed 's/  */\t/g' > /spectrum/GSCT/veronikap/DEEP/GCcount/BlEM/BlEM_${file}_retint.bed
	bedtools nuc -fi /spectrum/GSCT/REF/Human-GRCh38-release86/genome.fa -bed /spectrum/GSCT/veronikap/DEEP/GCcount/BlEM/BlEM_${file}_retint.bed -s > /spectrum/GSCT/veronikap/DEEP/GCcount/BlEM/GC_BlEM_${file}_retint.txt
done