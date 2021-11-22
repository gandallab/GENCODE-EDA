# Load libraries
using CSV, DataFrames, CodecZlib, Mmap, DataFramesMeta, Chain, Plots
using Statistics: mean
using FloatingTableView: browse
using GLM, StatsModels

# Download GENCODE v36
download("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz", "gencode.v36.annotation.gtf.gz")
#run(`curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz --output gencode.v36.annotation.gtf.gz`)

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
first(gencode, 5)
last(gencode, 5)
gencode[:, ["start", "end"]]
gencode[:, r".e"]

# Extract relevant information
for meta in ["gene_id", "gene_name", "gene_type"]
  gencode[:, "$(meta)"] = [getindex(m.captures, 1) for m in match.(Regex("$(meta) \"(.*?)\";"), gencode.info)]
end
#gencode_tx = gencode[gencode.feature .== "transcript", :]
#gencode_tx = filter(row -> row.feature == "transcript" && row.seqnames == "chr1", gencode)

# Extract transcript id's and transcript support lev  el (tsl)
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
df = @chain gencode begin
  @subset(:feature .== "gene")
  groupby(:seqnames)
  @combine(:count = length(:score))
  @subset(:seqnames .!= "chrX", :seqnames .!= "chrY", :seqnames .!= "chrM")
end

df.seqnames = parse.(Int64, replace.(df.seqnames, "chr" => ""))
lm = fit(LinearModel, @formula(count ~ seqnames), df)
summary(lm)
coef(lm)
f(x) = coef(lm)[2] * x + coef(lm)[1]
plot(df.seqnames, df.count, 
    smooth = true, 
    seriestype = :scatter, 
    title = "chromosome # vs # of genes", 
    linewidth = 8,
    linealpha = 0.5,
    line = (:black, 0.5, 6, :solid),
    size = (600, 300),
    label = nothing,
    xticks = (1:23),
    xlabel = "chromosome #",
    ylabel = "# of genes",
    background_color = :white,
    framestyle = :box, 
    widen = false,
    ylims = (0, 6000),
    grid = (0.3, :dot))
plot!(f, 0, 23, label = "linear model")

varinfo()
versioninfo()

df = DataFrame(transcript_id = unique(gencode.transcript_id)[2:4], test = 2)
@chain gencode begin
  @subset(.!ismissing.(:transcript_id))
  leftjoin(df, _; on = :transcript_id)
end

describe(gencode[:, ["start", "end"]], :mean, :std, :min, :q25, :median, 
  :q75, :max, :eltype, :nunique, :first, :last, :nmissing)

filter(endswith("id"), names(gencode))