' row_bind_tables.R

Given a set of tab separated values (TSV) files with identical headers,
combines the rows of the tables into a single output TSV file.

Usage: row_bind_tables.R ( -p PATHS | -i INPUT | -G GLOB ) -o OUTPUT

Options:
    -p --paths PATHS            Comma separated list of paths to TSV files
    -i --input INPUT            Path to a file containing paths to tables, one per line
    -G --glob GLOB              Globstring matching paths
    -o --output OUTPUT          Path to output

Examples:
    Rscript basic-functions/row_bind_tables.R -p file1.tsv,file2.tsv,file3.tsv -o combined_file.tsv
    Rscript basic-functions/row_bind_tables.R -G file*.tsv -o combined_file.tsv
' -> doc

library(docopt)
library(plyr)
library(readr)
library(doParallel)
args <- docopt(doc)

if ( ! is.null(args[['paths']]) ) {
    paths = strsplit(args[['paths']], ',')[[1]]
} else if (! is.null(args[['input']])) {
    paths = readLines(args[['input']])
} else if ( ! is.null(args[['glob']]) ) {
    paths = Sys.glob(args[['glob']])
} else {
    stop('Must provide one of --paths, --input, or --glob.')
}

message('Merging files')

registerDoParallel()

options(readr.show_progress = FALSE)

output <- ddply(data.frame(paths), 'paths', function(z) {
      path = as.character(z$paths)
      message(paste0('Reading file: ', path))
      suppressMessages(read_tsv(path))
}, .parallel = TRUE)

message('Done merging')

write_tsv(output, args[['output']])

message(paste0('Output written to ', args[['output']]))
