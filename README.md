# Demultiplexing pipeline for Illumina data

This pipeline is based in [USEARCH](https://www.drive5.com/usearch/) that is a propetary software. Please, read the software policy before using it. 

This short pipeline is supposed to run on HPCC mounting the SLURM system.

You will need to:
* clone the repo `git clone git@github.com:Gian77/Demultiplex-illumina.git`;
* put all your files in the `rawdata` directory, including the barcode files.
* modify the `config.yml` file adding the full PATHs to your files;
* run the pipeline `sbatch demultiplex.sb`.

All generated files will be in the `output_dir/`.

__NOTE__
* To increased the alllocated time and CPUs you should modify the `demultiplex.sb` script.
* To generate a reverse complemented version of the `barcodes.fasta` you can use this online tool https://www.bioinformatics.org/sms2/rev_comp.html or write a custom script that does the job for you.





