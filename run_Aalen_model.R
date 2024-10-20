library(timereg)

train = read.csv('tmp/train.csv', row.names=1, check.names=F)
metadata = read.csv('tmp/metadata.csv', row.names=1, check.names=F)

res = list()
i = 0
for (gene in colnames(train)) {
    df = cbind(train[gene], metadata[rownames(train), ])
    df = df[complete.cases(df), ]
    colnames(df) = c('gene', colnames(metadata))
    g = 'const(gene)'
    for (covariate in tail(colnames(df), n=-3)) {
        g = paste0(g, '+const(', covariate, ')')
    }
    out = aalen(as.formula(paste('Surv(time_to_event,event)', '~', g)), df, n.sim=100)
    out = coef.aalen(out)[1,,drop=F]
    rownames(out) = c(gene)
    i = i+1
    res[[i]] = out
}

res = do.call(rbind, res)
write.table(res, 'tmp/results.txt', sep='\t', quote=F, col.names=NA)
