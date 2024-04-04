for file in BlTN_Hf02 BlTN_Hf03 BlTN_Hf04 BlCM_Hf02 BlCM_Hf03 BlCM_Hf04 BlEM_Hf02 BlEM_Hf03 BlEM_Hf04
do
	awk '{gsub(/^chr/,""); print}' /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.nonretint.flankexons3ss.txt > /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.nonretint.flankexons3ss.nochr.txt #removing chr from the chromosome name
	sed 's/  */\t/g' /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.nonretint.flankexons3ss.nochr.txt > /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.nonretint.flankexons3ss.bed
	bedtools nuc -fi /spectrum/GSCT/REF/Human-GRCh38-release86/genome.fa -bed /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.nonretint.flankexons3ss.bed -s > /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/GC_${file}.nonretint.flankexons3ss.txt
done

for file in BlTN_Hf02 BlTN_Hf03 BlTN_Hf04 BlCM_Hf02 BlCM_Hf03 BlCM_Hf04 BlEM_Hf02 BlEM_Hf03 BlEM_Hf04
do
	awk '{gsub(/^chr/,""); print}' /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.nonretint.flankexons5ss.txt > /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.nonretint.flankexons5ss.nochr.txt #removing chr from the chromosome name
	sed 's/  */\t/g' /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.nonretint.flankexons5ss.nochr.txt > /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.nonretint.flankexons5ss.bed
	bedtools nuc -fi /spectrum/GSCT/REF/Human-GRCh38-release86/genome.fa -bed /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.nonretint.flankexons5ss.bed -s > /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/GC_${file}.nonretint.flankexons5ss.txt
done

for file in BlTN_Hf02 BlTN_Hf03 BlTN_Hf04 BlCM_Hf02 BlCM_Hf03 BlCM_Hf04 BlEM_Hf02 BlEM_Hf03 BlEM_Hf04
do
	awk '{gsub(/^chr/,""); print}' /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.retint.flankexons3ss.txt > /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.retint.flankexons3ss.nochr.txt #removing chr from the chromosome name
	sed 's/  */\t/g' /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.retint.flankexons3ss.nochr.txt > /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.retint.flankexons3ss.bed
	bedtools nuc -fi /spectrum/GSCT/REF/Human-GRCh38-release86/genome.fa -bed /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.retint.flankexons3ss.bed -s > /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/GC_${file}.retint.flankexons3ss.txt
done

for file in BlTN_Hf02 BlTN_Hf03 BlTN_Hf04 BlCM_Hf02 BlCM_Hf03 BlCM_Hf04 BlEM_Hf02 BlEM_Hf03 BlEM_Hf04
do
	awk '{gsub(/^chr/,""); print}' /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.retint.flankexons5ss.txt > /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.retint.flankexons5ss.nochr.txt #removing chr from the chromosome name
	sed 's/  */\t/g' /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.retint.flankexons5ss.nochr.txt > /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.retint.flankexons5ss.bed
	bedtools nuc -fi /spectrum/GSCT/REF/Human-GRCh38-release86/genome.fa -bed /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/${file}.retint.flankexons5ss.bed -s > /spectrum/GSCT/veronikap/DEEP/GCcount/Tcells/GC_${file}.retint.flankexons5ss.txt
done