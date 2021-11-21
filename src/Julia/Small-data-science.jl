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

# Extract relevant information
for meta in ["gene_id", "gene_name", "gene_type"]
  gencode[:, "$(meta)"] = [getindex(m.captures, 1) for m in match.(Regex("$(meta) \"(.*?)\";"), gencode.info)]
end
#gencode_tx = gencode[gencode.feature .== "transcript", :]
#gencode_tx = filter(row -> row.feature == "transcript" && row.seqnames == "chr1", gencode)

# Extract transcript id's and transcript support level (tsl)
for meta in ["transcript_id", "transcript_support_level"]
  storage = Vector(undef, size(gencode)[1])
  for (n, m) in enumerate(match.(Regex("$(meta) \"(.*?)\";"), gencode.info))
    if m == nothing 
      storage[n] = missing
    else
      storage[n] = getindex(m.captures, 1)
    end
  end
  gencode[:, "$(meta)"] = storage
end

# 32 exons with length of 1 bp
@chain gencode begin
  @transform(:length = :end - :start .+ 1)
  @subset(:length .== 1, :feature .== "exon")
  @combine(:count = length(:score))
end

# Number of genes by different source
@chain gencode begin
  @subset(:feature .== "gene")
  groupby(:source)
  @combine(:count = length(:score))
end

# Number of genes by gene type
@chain gencode begin
  @subset(:feature .== "gene")
  groupby(:gene_type)
  @combine(:count = length(:score))
end

# Average gene length
@chain gencode begin
  @subset(:feature .== "gene")
  @transform(:length = :end - :start .+ 1)
  @combine(:avg = mean(:length))
end

@chain gencode begin
  @subset(:feature .== "gene")
  @transform(:length = :end - :start .+ 1)
  @orderby(-:length)
  @select(:gene_name, :length)
end

# Number of genes per chromosome
@chain gencode begin
  @subset(:feature .== "gene")
  groupby(:seqnames)
  @combine(:count = length(:score))
end

varinfo()
versioninfo()