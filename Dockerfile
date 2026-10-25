# Reproducible environment for CNAppR: R 4.4 + Bioconductor 3.20 + the packages
# needed to run the full BAM pipeline (CBS and ichorCNA), GSEA and survival analysis.
#
#   docker build -t cnappr .
#   docker run --rm -it -v "$PWD":/data cnappr
#
# Inside the container, CNAppR is already installed: library(CNAppR)
FROM bioconductor/bioconductor_docker:RELEASE_3_20

RUN R -q -e 'BiocManager::install(c( \
      "Rsamtools", "GenomicRanges", "IRanges", "GenomeInfoDb", "BSgenome", \
      "BSgenome.Hsapiens.UCSC.hg38", "DNAcopy", "HMMcopy", "fgsea", "msigdbr", \
      "survival", "survminer", "randomForest", "caret", "dplyr", "valr"), \
      ask = FALSE, update = FALSE)'

# Installs the version published on GitHub (master branch).
RUN R -q -e 'remotes::install_github("Albertotorrejon/CNAppR", upgrade = "never")'

WORKDIR /data
CMD ["R"]
