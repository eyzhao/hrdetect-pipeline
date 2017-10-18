' 
Usage: flank_indels.R -v VCF -o OUTPUT

Options:
    -v --vcf VCF            Variant call file with indels
    -o --output OUTPUT      Path to output table of annotated indels
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)
library(BSgenome)
library(VariantAnnotation)

rev_char <- function(char) {
    sapply(char, function(a) {
           paste(rev(substring(a,1:nchar(a),1:nchar(a))),collapse="")
    })
}

get_match_length <- function(x, y) { 
    min_length <- min(nchar(x), nchar(y))
    if (min_length == 0) {
        return(0) 
    } else { 
        return(sum(sapply(1:min_length, function(i) {
            substr(x, 1, i) == substr(y, 1, i) 
        })))
    }
}

get_overlap_length <- function(seq1, seq2) {
  overlap_matrix <- sapply(1:nchar(seq1), function(i) {
    subseq1 <- substring(seq1, i, nchar(seq1))
    sapply(1:nchar(seq2), function(j) {
      subseq2 <- substring(seq2, j, nchar(seq2))
      get_match_length(subseq1, subseq2)
    })
  })
  
  return(max(overlap_matrix))
}

parse_indel <- function(indel_vcf_path) {
    vcf <- readVcf(indel_vcf_path)
    vr <- rowRanges(vcf)
    print(vr)
    rowsToKeep <- elementNROWS(vr$ALT) == 1
    vr$ALT <- unlist(CharacterList(vr$ALT))
    all_indels <- as_tibble(vr[rowsToKeep, ]) %>%
        mutate(
            seqnames = gsub('chr', '', seqnames),
            seqnames = gsub('23', 'X', seqnames),
            seqnames = gsub('24', 'Y', seqnames)
        )

    indels <- all_indels[paste0('chr', all_indels$seqnames) %in% seqlevels(hg19)
                         & nchar(all_indels$REF) > nchar(all_indels$ALT), ]
    return(indels)
}

get_flanking_sequence <- function(indel_vcf_path) {
    print(indel_vcf_path)
    indels <- parse_indel(as.character(indel_vcf_path))
    print(paste('Number of indels:', dim(indels)[1]))
    print('retrieving flanking sequences')
    five_prime_flank <- getSeq(hg19, paste0('chr', indels$seqnames), indels$start - 25, indels$start)
    three_prime_flank <- getSeq(hg19, paste0('chr', indels$seqnames), indels$end + 1, indels$end + 26)
    deleted <- as.character(sapply(indels$REF, function(z) {substring(z, 2, nchar(z))}))
    return(cbind(indels, tibble(deleted = as.character(deleted),
                                five_prime_flank = as.character(five_prime_flank), 
                                three_prime_flank = as.character(three_prime_flank))))
}

get_match_length <- function(x, y) { 
    min_length <- min(nchar(x), nchar(y))
    if (min_length == 0) {
        return(0) 
    } else { 
        return(sum(sapply(1:min_length, function(i) {
            substr(x, 1, i) == substr(y, 1, i) 
        })))
    }
}

microhomology_summary <- function(indel_table) {
    indel_table <- indel_table %>% filter(
        ! is.na(deleted),
        ! is.na(three_prime_flank),
        ! is.na(five_prime_flank)
    ) 

    print('computing match statistics')

    three_prime <- apply(indel_table, 1, function(z) {
        get_match_length(z['deleted'], z['three_prime_flank'])
    })

    five_prime <- apply(indel_table, 1, function(z) {
        get_match_length(rev_char(z['deleted']), rev_char(z['five_prime_flank']))
    })

    max_match <- apply(cbind(three_prime, five_prime), 1, max)

    indel_table$microhomology_length <- max_match
    indel_table$is_microhomology <- max_match > 2

    return(indel_table)
}

hg19 <- getBSgenome('BSgenome.Hsapiens.UCSC.hg19')

indel_flanked <- get_flanking_sequence(args[['vcf']]) %>% microhomology_summary()
write_tsv(indel_flanked, path = args[['output']])
print(paste('Wrote INDEL flanking data to', args[['output']]))
