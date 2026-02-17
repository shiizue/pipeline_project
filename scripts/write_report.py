with open("PipelineReport.txt", "w") as f:
    for input in [
        snakemake.input.cds,
        snakemake.input.sleuth,
    ]:
        with open(input) as r:
            f.write(r.read())
        f.write("\n")
