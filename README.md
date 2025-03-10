# Demultiplexing pipeline for Illumina data

This pipeline is based on [USEARCH](https://github.com/rcedgar/usearch12) that is a propetary software. Please, read the software policy before using it. 

This short pipeline is supposed to run on HPCC mounting the SLURM system.

You will need to:
* clone the repo `git clone git@github.com:Gian77/Demultiplex-USEARCH.git`;
* put all your files in the `rawdata` directory, including the barcode files.
* modify the `config.yml` file adding the full PATHs to your files;
* run the pipeline `sbatch demultiplex-usearch.sb`.

All generated files will be in the `output_dir/`.

__NOTE__
* If you have a very large datasets, to increased the alllocated time, RAM, and CPU requirements you should modify the `demultiplex-usearch.sb` script. 

```
#SBATCH --time=00:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
```

The barcode file should be structured for [USEARCH](https://github.com/rcedgar/usearch12), as showed below, and saved it as a `.fasta`. The name of sample should be as simple as possible for convenience. 

```
>bc1
TACGACTCTG
>bc2
TTCTTAACGC
>bc3
ACCCAGTATG
>bc4
AACGGCTGGA
...
```

Make sure you change the USEARCH executable in the `demultiplex-usearch.sb` script before running or add the full path to the `config.yml` file.
