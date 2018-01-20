' compute_hrdetect.R

Usage: compute_hrdetect.R -i INPUT -o OUTPUT

Options:
    -i --input INPUT       Path to input HRDetect table
    -o --output OUTPUT     Path to output table with HRDetect data
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

makeglm <- function(formula, ..., family, data=NULL) {
    dots <- list(...)
    out<-list()
    tt <- terms(formula, data=data)
    if(!is.null(data)) {
        mf <- model.frame(tt, data)
        vn <- sapply(attr(tt, "variables")[-1], deparse)
 
        if((yvar <- attr(tt, "response"))>0)
            vn <- vn[-yvar]
            xlvl <- lapply(data[vn], function(x) if (is.factor(x))
           levels(x)
        else if (is.character(x))
           levels(as.factor(x))
        else
            NULL)
        attr(out, "xlevels") <- xlvl[!vapply(xlvl,is.null,NA)]
        attr(tt, "dataClasses") <- sapply(data[vn], stats:::.MFclass)
    }
    out$terms <- tt
    coef <- numeric(0)
    stopifnot(length(dots)>1 & !is.null(names(dots)))
    for(i in seq_along(dots)) {
        if((n<-names(dots)[i]) != "") {
            v <- dots[[i]]
            if(!is.null(names(v))) {
                coef[paste0(n, names(v))] <- v
            } else {
                stopifnot(length(v)==1)
                coef[n] <- v
            }
        } else {
            coef["(Intercept)"] <- dots[[i]]
        }   
    }
    out$coefficients <- coef
    out$rank <- length(coef)
    if (!missing(family)) {
        out$family <- if (class(family) == "family") {
            family
        } else if (class(family) == "function") {
            family()
        } else if (class(family) == "character") {
            get(family)()
        } else {
            stop(paste("invalid family class:", class(family)))
        }
        out$qr <- list(pivot=seq_len(out$rank))
        out$deviance <- 1
        out$null.deviance <- 1
        out$aic <- 1
        class(out) <- c("glm","lm")
    } else {
        class(out) <- "lm"
        out$fitted.values <- predict(out, newdata=dd)
        out$residuals <- out$mf[attr(tt, "response")] - out$fitted.values
        out$df.residual <- nrow(data) - out$rank
        out$model <- data
        #QR doesn't work
    }
    out
}

get_hrdetect_score <- function(hrd_table) {
  if (! 'outcome' %in% colnames(hrd_table)) {
    hrd_table$outcome = rbernoulli(n = dim(hrd_table)[1], p = 0.5)
  }
  hrdetect_model <- makeglm(
    outcome ~ deletion_microhomology_proportion +
        snv_signature_3 + 
        rearrangement_signature_3 +
        rearrangement_signature_5 + 
        hrd_score + 
        snv_signature_8,
    family=binomial, data=hrd_table, 
    -3.364,
    deletion_microhomology_proportion = 2.398,
    snv_signature_3 = 1.611,
    rearrangement_signature_3 = 1.153,
    rearrangement_signature_5 = 0.847,
    hrd_score = 0.667,
    snv_signature_8 = 0.091
  )
  
  return(predict(hrdetect_model, newdata=hrd_table, type="response"))
}

input_table <- read_tsv(args[['input']]) %>%
    select(-data_type) %>%
    group_by(variable) %>%
    mutate(
        value = log(value + 1),
        value = (value - mean(value)) / sd(value)
    ) %>%
    ungroup() %>%
    spread(variable, value) %>%
    select(
        patient,
        sample,
        deletion_microhomology_proportion, 
        hrd_score = total, 
        snv_signature_3 = `Signature 3`, 
        snv_signature_8 = `Signature 8`, 
        rearrangement_signature_3 = `Rearrangement Signature 3`, 
        rearrangement_signature_5 = `Rearrangement Signature 5`
    )

missing_data_rows <- input_table %>% apply(1, function(row) { any(is.na(row)) })

if (any(missing_data_rows)) {
    message('Some cases are missing data and will not be included in the output.')
    print(input_table[missing_data_rows, ])
}

output_table <- input_table[ ! missing_data_rows, ]
output_table$HRDetect <- get_hrdetect_score(output_table)

output_table %>%
    write_tsv(args[['output']])
