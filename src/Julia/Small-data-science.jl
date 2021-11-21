# Load libraries
using CSV, DataFrames, CodecZlib, Mmap, DataFramesMeta, Chain, Plots
using Statistics: mean
using FloatingTableView: browse

# Download GENCODE v36
run(`curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz --output gencode.v36.annotation.gtf.gz`)

# Load GENCODE v36
@time gencode = CSV.File(transcode(GzipDecompressor, 
    Mmap.mmap("gencode.v36.annotation.gtf.gz")),
    delim = "\t", skipto = 6, 
    header = ["seqnames", "source", "feature", "start", "end", 
    "score", "strand", "phase", "info"]) |> DataFrame

#browse(gencode)
size(gencode)
names(gencode)
typeof(gencode)
gencode[:, ["start", "end"]]
gencode[:, r".e"]

# 32 exons with length of 1 bp
@chain gencode begin
    @transform(length = :end - :start .+ 1)
    @where(:length .== 1, :feature .== "exon")
    @combine(count = length(:score))
end

# Extract transcript support level (tsl)
gencode_tx = gencode[gencode.feature .== "transcript", :]
#gencode_tx = filter(row -> row.feature == "transcript", gencode)
#gencode_tx = filter(row -> row.feature == "transcript" && row.seqnames == "chr1", gencode)
match.(r"transcript_support_level \"(.*?)\";", gencode_tx.info[1]).captures
match.(r"gene_id \"(.*?)\";", gencode_tx.info[1]).captures
match.(r"gene_name \"(.*?)\";", gencode_tx.info[1]).captures
match.(r"transcript_id \"(.*?)\";", gencode_tx.info[1]).captures
match.(r"gene_type \"(.*?)\";", gencode_tx.info[1]).captures

